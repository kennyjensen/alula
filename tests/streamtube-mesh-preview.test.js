import test from 'node:test';
import assert from 'node:assert/strict';
import { buildStreamtubeMeshPreview } from '../src/euler/streamtube-mesh-preview.js';
import { createInitialStreamtubeTopology } from '../src/geometry/streamtube-topology.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { slatMainFlap } from './fixtures/subcritical-assemblies.js';
import { solveStreamtubeAssembly } from '../src/euler/streamtube-result.js';
import { streamtubeAuditRegions, measureStreamtubeGridSpacing } from '../src/geometry/streamtube-grid-audit.js';

const flap = [{ points: naca4('2412', 160) }, { points: transform(naca4('0012', 160), { chord: .3, x: .94, y: -.08, angle: -15 }) }];

for (const [gridTubes, gridIntervals] of [[7, 16], [9, 16], [11, 8], [11, 16], [11, 32]]) {
  test(`default main/flap Build mesh: ${gridTubes} tubes and ${gridIntervals} surface intervals remain valid without Euler iteration`, () => {
    const events = [], r = solveStreamtubeAssembly({ elements: flap, alpha: 4, gridTubes, gridIntervals }, {
      meshOnly: true, onMesh: (mesh, stage) => events.push(stage),
      onIteration: () => assert.fail('Build mesh must not start an Euler iteration'),
    });
    assert.equal(r.status, 'mesh-ready'); assert.deepEqual(events, ['initial', 'initial']);
    const { mesh } = r;
    assert.equal(mesh.quality.valid, true); assert.equal(mesh.initialization.flowSolved, false);
    assert.equal(mesh.initialization.cutStationSpacing, 'physical-x');
    assert.equal(mesh.initialization.gridSmoothing, undefined);
    assert.equal(mesh.initialization.gridRepair, undefined);
    assert.deepEqual(mesh.initialization.tubes, [gridTubes + 1, gridTubes + 2, gridTubes + 1]);
    assert.equal(mesh.initialization.gridRefinement.flowInitialized, false);
    assert.equal(mesh.flow, undefined);
    assert.deepEqual(mesh.initialization.elementOrder, [1, 0]);
    assert.ok(mesh.initialization.maxStreamfunctionDrift < 1e-7);
    // Audit every corner independently of the production mesh-quality code.
    const edges = new Map();
    for (const cell of mesh.cells) {
      assert.equal(cell.length, 4); assert.equal(new Set(cell).size, 4);
      const p = cell.map(id => mesh.vertices[id]);
      for (let j = 0; j < 4; j++) {
        const [a, b, c] = [p[j], p[(j + 1) % 4], p[(j + 2) % 4]];
        assert.ok((b.x - a.x) * (c.y - b.y) > (b.y - a.y) * (c.x - b.x), 'nonpositive corner');
        const u = cell[j], v = cell[(j + 1) % 4], key = u < v ? `${u}:${v}` : `${v}:${u}`;
        edges.set(key, (edges.get(key) ?? 0) + 1);
      }
    }
    assert.ok([...edges.values()].every(count => count === 1 || count === 2));
    assert.equal(mesh.vertices.length - edges.size + mesh.cells.length, -1, 'two body holes must be retained');
    if (gridTubes === 7 && gridIntervals === 16) {
      const regions = streamtubeAuditRegions(mesh), spacing = measureStreamtubeGridSpacing(regions);
      // Regression screens for the old 61:1 flap jump and 0.429-chord
      // incoming skew. These are not physical-flow acceptance thresholds.
      assert.ok(spacing.maximumAdjacentRatio < 4.1, 'initial station spacing regressed');
      assert.ok(spacing.minimumCornerSine > .45, 'initial cell skew regressed');
      const blocks = mesh.initialization.potentialCrosslines;
      const firstEdge = blocks.x.indexOf(Math.min(...Object.values(blocks.panelBlocks.rank)));
      assert.ok(firstEdge > 0);
      for (const row of regions[1].nodes.slice(0, firstEdge + 1))
        assert.ok(Math.abs(row[0].x - row.at(-1).x) < .05, 'incoming cut banks lost their station correspondence');
    }
  });
}

test('one/two/three-element quadrilateral previews have positive cells, connected cuts and the correct number of body holes', () => {
  for (const elements of [[{ points: naca4('0012', 160) }], flap, slatMainFlap([120, 240, 120]).elements]) {
    const mesh = buildStreamtubeMeshPreview({ elements, alpha: 4 });
    assert.equal(mesh.topology, 'intrinsic-quadrilateral-streamtubes');
    assert.equal(mesh.quality.valid, true); assert.ok(mesh.quality.minArea > 0); assert.ok(mesh.quality.minCornerSine > .1);
    assert.equal(mesh.initialization.flowSolved, false); assert.equal(mesh.initialization.groups, elements.length + 1);
    assert.ok(mesh.initialization.maxStreamfunctionDrift < 1e-7);
    const edges = new Map(), neighbors = mesh.cells.map(() => []), used = new Set();
    mesh.cells.forEach((cell, id) => {
      assert.equal(cell.length, 4); assert.equal(new Set(cell).size, 4);
      for (let j = 0; j < 4; j++) {
        const a = cell[j], b = cell[(j + 1) % 4]; used.add(a);
        assert.ok(Number.isInteger(a) && a >= 0 && a < mesh.vertices.length);
        const key = a < b ? `${a}:${b}` : `${b}:${a}`, existing = edges.get(key);
        if (existing) { assert.equal(existing.length, 1); neighbors[id].push(existing[0]); neighbors[existing[0]].push(id); existing.push(id); }
        else edges.set(key, [id]);
      }
    });
    assert.equal(used.size, mesh.vertices.length);
    assert.equal(mesh.vertices.length - edges.size + mesh.cells.length, 1 - elements.length, 'wrong physical hole/cut topology');
    const seen = new Set([0]), stack = [0];
    while (stack.length) for (const next of neighbors[stack.pop()]) if (!seen.has(next)) { seen.add(next); stack.push(next); }
    assert.equal(seen.size, mesh.cells.length, 'fluid regions do not connect across cuts');
    for (const e of elements) assert.ok(mesh.vertices.some(p => Math.hypot(p.x - e.points[0].x, p.y - e.points[0].y) < 1e-12));
  }
});

test('topology controls and element permutations preserve the actual input geometry', () => {
  const input = { elements: flap, alpha: 4 }, before = structuredClone(input);
  const a = createInitialStreamtubeTopology(input), b = createInitialStreamtubeTopology({ ...input, elements: [...flap].reverse() });
  assert.deepEqual(input, before);
  for (let i = 0; i < a.bodies.length; i++) {
    assert.deepEqual(a.bodies[i].points, b.bodies[i].points);
    assert.deepEqual(a.bodies[i].surfaceFractions, b.bodies[i].surfaceFractions);
    assert.equal(a.bodies[i].element, 1 - b.bodies[i].element);
  }
  for (const controls of [{ tubes: 0 }, { surfaceIntervals: 5 }, { tubeGrowth: NaN }, { padding: .1 }])
    assert.throws(() => createInitialStreamtubeTopology(input, controls), /controls/);
  assert.throws(() => createInitialStreamtubeTopology({ elements: [] }), /controls/);
});
