#!/usr/bin/env python3
"""Extract exact closed wall loops from the attributed SAAB grid, no fitting.

This extractor is original GPL-2.0-or-later code. Upstream data attribution and
rights are separate; see README.md. All source float values are retained.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import json
import sys
import zipfile

root = Path(__file__).resolve().parent
archive = root / 'UFR3-01_nlr7301grid.zip'
member = 't/nlrgap26.dat'
with zipfile.ZipFile(archive) as source:
    original = source.read(member)
tokens = iter(original.decode('ascii').split())
blocks = int(next(tokens))
assert blocks == 9
edges = defaultdict(list)
dimensions = []
for block in range(1, blocks + 1):
    ni, nj = int(next(tokens)), int(next(tokens))
    dimensions.append([ni, nj])
    grid = [[(float(next(tokens)), float(next(tokens))) for _ in range(ni)] for _ in range(nj)]
    for side, points in [('imin', [r[0] for r in grid]), ('imax', [r[-1] for r in grid]),
                         ('jmin', grid[0]), ('jmax', grid[-1])]:
        for k, (a, b) in enumerate(zip(points, points[1:]), 1):
            assert a != b
            edges[tuple(sorted((a, b)))].append([block, side, k])
assert next(tokens, None) is None
assert all(len(v) in (1, 2) for v in edges.values())
adjacency = defaultdict(list)
for (a, b), copies in edges.items():
    if len(copies) == 1:
        adjacency[a].append(b)
        adjacency[b].append(a)
assert all(len(v) == 2 for v in adjacency.values())
visited, loops = set(), []
for start in adjacency:
    if start in visited:
        continue
    points, previous, current = [start], None, start
    while True:
        nxt = next(v for v in adjacency[current] if v != previous)
        if nxt == start:
            break
        assert nxt not in visited
        visited.add(current)
        points.append(nxt)
        previous, current = current, nxt
    visited.update(points)
    loops.append(points)
assert len(loops) == 3
area = lambda p: .5 * sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(p, p[1:] + p[:1]))
loops.sort(key=lambda p: abs(area(p)))
flap, main, farfield = loops
assert [len(main), len(flap), len(farfield)] == [456, 248, 1216]
elements = []
for name, loop, upper_te, lower_te in [
    ('Main element', main, (.9436, .01499), (.9436, .0141)),
    ('Flap', flap, (1.201771, -.1037033), (1.201377, -.104784)),
]:
    if area(loop) < 0:
        loop.reverse()
    k = loop.index(upper_te)
    loop = loop[k:] + loop[:k]
    assert lower_te in loop
    spans = defaultdict(list)
    for a, b in zip(loop, loop[1:] + loop[:1]):
        block, side, index = edges[tuple(sorted((a, b)))][0]
        spans[(block, side)].append(index)
    elements.append({'name': name, 'points': [list(p) for p in loop + loop[:1]],
        'area': area(loop), 'upperTrailingEdge': upper_te, 'lowerTrailingEdge': lower_te,
        'trailingEdgeBaseLength': sum((a-b)**2 for a,b in zip(upper_te, lower_te))**.5,
        'sourceEdges': [{'block': b, 'side': s, 'firstSegment': min(v), 'lastSegment': max(v), 'count': len(v)}
                        for (b, s), v in sorted(spans.items())]})
report = {'sourceArchive': archive.name, 'sourceArchiveHash': sha256(archive.read_bytes()).hexdigest(),
    'sourceMember': member, 'sourceMemberHash': sha256(original).hexdigest(), 'blockDimensions': dimensions,
    'method': 'Exact cancellation of duplicate block-boundary segments; the two bounded wall loops are oriented CCW and cyclically start at their supplied upper TE corners. Repeated closure endpoints are appended; no existing point moves or is removed.',
    'edgeMultiplicity': dict(Counter(len(v) for v in edges.values())),
    'unpairedVertexDegrees': dict(Counter(len(v) for v in adjacency.values())),
    'loopPointCounts': {'main': len(main), 'flap': len(flap), 'farfield': len(farfield)},
    'referenceChord': 1, 'elements': elements}
outputs = {'geometry.json': json.dumps(report, indent=2) + '\n'}
lines = ['NLR 7301 SAAB FLOWNET gap 2.6% / flap 20deg; exact grid walls including finite TE bases', '-5 5 -5 5']
for k, element in enumerate(elements):
    if k:
        lines.append('999 999')
    lines += [f'{x!r} {y!r}' for x, y in element['points']]
outputs['nlr7301-gap26.dat'] = '\n'.join(lines) + '\n'
assert sys.argv[1:] in ([], ['--check'])
for name, content in outputs.items():
    if sys.argv[1:] == ['--check']:
        assert (root / name).read_text() == content, f'{name} differs from the exact source extraction'
    else:
        (root / name).write_text(content)
print(json.dumps({k: v for k, v in report.items() if k != 'elements'}, indent=2))
for element in elements:
    print(json.dumps({k: v for k, v in element.items() if k != 'points'}))
