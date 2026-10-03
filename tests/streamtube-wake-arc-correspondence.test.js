import test from 'node:test';
import assert from 'node:assert/strict';
import { initializeStreamtubeWakeArcCorrespondence as initialize } from '../src/euler/tests/streamtube-wake-arc-correspondence.js';
import { initializeStreamtubeWakeCorrespondence as legacy } from '../src/euler/streamtube-wake-correspondence.js';
import { streamtubeWakeGap } from '../src/euler/streamtube-wake-geometry.js';
const close = (a, b, tol = 3e-12) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);
const midpoint = (p, q) => ({ x: .5 * (p.x + q.x), y: .5 * (p.y + q.y) });
const center = (nodes, b, i) => midpoint(nodes[b][i].at(-1), nodes[b + 1][i][0]);
const distance = (p, q) => Math.hypot(p.x - q.x, p.y - q.y);
const normal = (p, q) => { const ds = distance(p, q); return { x: -(q.y - p.y) / ds, y: (q.x - p.x) / ds }; };
const gap = (f, nodes, b, i) => { const ids = [i - 1, i, Math.min(i + 1, f.layout.nx)];
  return streamtubeWakeGap(ids.map(k => nodes[b][k].at(-1)), ids.map(k => nodes[b + 1][k][0])).gap; };
function fixture({ elements = 1, curved = false, zeroGap = false } = {}) {
  const nx = 7, tubes = Array(elements + 1).fill(3), bodies = Array.from({ length: elements }, (_, b) => ({ leadingIndex: 1, trailingIndex: b + 2, element: elements - b - 1 }));
  const layout = { nx, tubes, bodies, elements, independentWakeBanks: true };
  const massFractions = tubes.map((_, g) => [0, .08 + .04 * g, .69 + .03 * g, 1]);
  const paths = [], banks = [], referenceCutPaths = [];
  for (let b = 0; b < elements; b++) {
    const te = bodies[b].trailingIndex, length = nx - te;
    const points = Array.from({ length: nx + 1 }, (_, i) => {
      const k = i - te, z = k < 0 ? k / length : (k / length) ** 2;
      return { x: 2 * z, y: b * 3 + (curved ? .3 * z + .22 * Math.sin(z * 2) : 0) };
    });
    paths.push(points);
    referenceCutPaths.push(points.map((p, i) => ({ x: (i - te) / length * 3, y: 0 })));
    banks.push(points.map((p, i) => {
      const n = normal(points[Math.max(0, i - 1)], points[Math.min(nx, i + 1)]), t = { x: n.y, y: -n.x };
      const thickness = zeroGap ? 0 : .0625 + Math.max(0, i - te) / 512;
      const tau = i <= te ? 1 / 64 : 0;
      const d = { x: thickness * n.x + tau * t.x, y: thickness * n.y + tau * t.y };
      return [{ x: p.x - d.x / 2, y: p.y - d.y / 2 }, { x: p.x + d.x / 2, y: p.y + d.y / 2 }];
    }));
  }
  const nodes = tubes.map((n, g) => Array.from({ length: nx + 1 }, (_, i) => {
    const lower = g === 0 ? { x: paths[0][i].x, y: -3 } : banks[g - 1][i][1];
    const upper = g === elements ? { x: paths.at(-1)[i].x, y: elements * 3 } : banks[g][i][0];
    return Array.from({ length: n + 1 }, (_, j) => ({ x: lower.x + j / n * (upper.x - lower.x), y: lower.y + j / n * (upper.y - lower.y) }));
  }));
  return { nodes, layout, massFractions, referenceCutPaths };
}
function invariants(f, r) {
  for (let b = 0; b < f.layout.elements; b++) {
    const te = f.layout.bodies[b].trailingIndex, nx = f.layout.nx;
    for (let i = 0; i <= te; i++) {
      assert.deepEqual(r.nodes[b][i].at(-1), f.nodes[b][i].at(-1));
      assert.deepEqual(r.nodes[b + 1][i][0], f.nodes[b + 1][i][0]);
    }
    close(center(r.nodes, b, nx).x, center(f.nodes, b, nx).x);
    close(center(r.nodes, b, nx).y, center(f.nodes, b, nx).y);
    for (let i = te + 1; i <= nx; i++) {
      close(gap(f, r.nodes, b, i), gap(f, f.nodes, b, i));
      const a = center(r.nodes, b, i - 1), c = center(r.nodes, b, i), ds = distance(a, c);
      for (const [g, j] of [[b, f.layout.tubes[b]], [b + 1, 0]]) {
        const p = r.nodes[g][i - 1][j], q = r.nodes[g][i][j];
        const forward = ((q.x - p.x) * (c.x - a.x) + (q.y - p.y) * (c.y - a.y)) / ds;
        assert.ok(forward > 0); close(forward, ds);
      }
    }
  }
  for (let i = 0; i <= f.layout.nx; i++) {
    assert.deepEqual(r.nodes[0][i][0], f.nodes[0][i][0]);
    assert.deepEqual(r.nodes.at(-1)[i].at(-1), f.nodes.at(-1)[i].at(-1));
  }
}

test('straight current wake follows analytic reference fractions and preserves original-index gaps', () => {
  const f = fixture(), r = initialize(f), te = f.layout.bodies[0].trailingIndex;
  invariants(f, r);
  for (let i = te; i <= f.layout.nx; i++) {
    close(center(r.nodes, 0, i).x, 2 * (i - te) / (f.layout.nx - te));
    close(center(r.nodes, 0, i).y, 0);
  }
  assert.ok(r.diagnostics.bodies[0].firstCenterAdvanceAfter > r.diagnostics.bodies[0].firstCenterAdvanceBefore);
});

test('curved polyline uses its cumulative arc, not x or straight endpoint interpolation', () => {
  const f = fixture({ curved: true }), r = initialize(f), b = r.diagnostics.bodies[0];
  invariants(f, r);
  const te = f.layout.bodies[0].trailingIndex, input = Array.from({ length: f.layout.nx - te + 1 }, (_, k) => center(f.nodes, 0, te + k));
  const arc = [0]; input.slice(1).forEach((p, i) => arc.push(arc.at(-1) + distance(input[i], p)));
  let differentFromChord = false;
  for (let k = 1; k < input.length - 1; k++) {
    const target = arc.at(-1) * k / (input.length - 1);
    let j = 0; while (arc[j + 1] < target) j++;
    const t = (target - arc[j]) / (arc[j + 1] - arc[j]);
    const c = center(r.nodes, 0, te + k);
    close(c.x, input[j].x * (1 - t) + input[j + 1].x * t);
    close(c.y, input[j].y * (1 - t) + input[j + 1].y * t);
    differentFromChord ||= Math.abs(c.y - (input[0].y + k / (input.length - 1) * (input.at(-1).y - input[0].y))) > 1e-3;
  }
  assert.ok(differentFromChord); close(b.originalCenterArcLength, arc.at(-1));
});

test('two bodies with different trailing stations blend both shared-passage corrections once by physical mass', () => {
  const f = fixture({ elements: 2, curved: true }), before = JSON.stringify(f), r = initialize(f);
  invariants(f, r); assert.equal(JSON.stringify(f), before);
  for (let g = 0; g < f.nodes.length; g++) for (let i = 0; i <= f.layout.nx; i++) for (let j = 1; j < f.layout.tubes[g]; j++) {
    const eta = f.massFractions[g][j];
    for (const key of ['x', 'y']) close(r.nodes[g][i][j][key], f.nodes[g][i][j][key]
      + (1 - eta) * (r.nodes[g][i][0][key] - f.nodes[g][i][0][key])
      + eta * (r.nodes[g][i].at(-1)[key] - f.nodes[g][i].at(-1)[key]));
  }
  assert.equal(r.diagnostics.bodies[0].trailingIndex, 2); assert.equal(r.diagnostics.bodies[1].trailingIndex, 3);
  const g = 1, i = 6, j = 1, geometricEta = 1 / 3;
  const wrong = f.nodes[g][i][j].x + (1 - geometricEta) * (r.nodes[g][i][0].x - f.nodes[g][i][0].x)
    + geometricEta * (r.nodes[g][i].at(-1).x - f.nodes[g][i].at(-1).x);
  assert.ok(Math.abs(r.nodes[g][i][j].x - wrong) > 1e-5, 'Index interpolation must not masquerade as mass interpolation.');
});

test('rigid rotation, translation and positive dimensional scaling commute with the construction', () => {
  const f = fixture({ elements: 2, curved: true }), a = initialize(f), angle = .48, scale = 3.7;
  const transform = p => ({ x: 13 + scale * (Math.cos(angle) * p.x - Math.sin(angle) * p.y),
    y: -8 + scale * (Math.sin(angle) * p.x + Math.cos(angle) * p.y) });
  const moved = structuredClone(f); moved.nodes = moved.nodes.map(g => g.map(row => row.map(transform)));
  moved.referenceCutPaths = moved.referenceCutPaths.map(row => row.map(transform));
  const b = initialize(moved);
  a.nodes.forEach((g, gi) => g.forEach((row, i) => row.forEach((p, j) => {
    const q = transform(p); close(q.x, b.nodes[gi][i][j].x, 2e-11); close(q.y, b.nodes[gi][i][j].y, 2e-11);
  })));
  close(b.diagnostics.minimumBankForwardAdvance, scale * a.diagnostics.minimumBankForwardAdvance);
});

test('reference length and frame do not change its normalized station fractions', () => {
  const f = fixture({ curved: true }), a = initialize(f), other = structuredClone(f);
  other.referenceCutPaths = other.referenceCutPaths.map(row => row.map(p => ({ x: 4 - 9 * p.y, y: 7 + 9 * p.x })));
  const b = initialize(other);
  a.nodes.forEach((g, gi) => g.forEach((row, i) => row.forEach((p, j) => {
    close(p.x, b.nodes[gi][i][j].x); close(p.y, b.nodes[gi][i][j].y);
  })));
});

test('zero-gap point aliases remain detached and input geometry is not mutated', () => {
  const f = fixture({ zeroGap: true });
  for (let i = f.layout.bodies[0].trailingIndex + 1; i <= f.layout.nx; i++) f.nodes[1][i][0] = f.nodes[0][i].at(-1);
  const old = JSON.stringify(f), r = initialize(f);
  assert.equal(JSON.stringify(f), old); invariants(f, r);
  assert.notEqual(r.nodes[0][4].at(-1), r.nodes[1][4][0]);
  r.nodes[0][0][0].x += 100; assert.equal(JSON.stringify(f), old);
});

test('malformed explicit reference, mass and geometry inputs reject without silent fallback', () => {
  for (const mutate of [
    f => delete f.referenceCutPaths,
    f => f.referenceCutPaths.pop(),
    f => f.referenceCutPaths[0].pop(),
    f => f.referenceCutPaths[0][4].x = NaN,
    f => f.referenceCutPaths[0][4] = { ...f.referenceCutPaths[0][3] },
    f => f.massFractions[0][1] = f.massFractions[0][2],
    f => f.massFractions[0][0] = .01,
    f => f.nodes[0][4][0].x = Infinity,
    f => f.layout.bodies[0].trailingIndex = f.layout.nx,
  ]) { const f = fixture(); mutate(f); assert.throws(() => initialize(f)); }
  const f = fixture();
  for (const [g, j] of [[0, 3], [1, 0]]) f.nodes[g][4][j] = { ...f.nodes[g][3][j] };
  assert.throws(() => initialize(f), /current arc/);
});

test('legacy correspondence entry point stays exact and independently available when reference arc is omitted', () => {
  const f = fixture({ curved: true });
  const { referenceCutPaths, ...oldInput } = f;
  const a = legacy(oldInput), b = legacy({ ...oldInput, referenceCutPaths });
  assert.deepEqual(a, b);
  assert.notDeepEqual(initialize(f).nodes, a.nodes);
});
