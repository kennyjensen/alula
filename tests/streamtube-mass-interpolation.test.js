import test from 'node:test';
import assert from 'node:assert/strict';
import { refineStreamtubeMassCoordinates as refine } from '../src/geometry/streamtube-mass-interpolation.js';

const levels = [0, .02, .07, .15], masses = levels.slice(1).map((v, i) => v - levels[i]);
const counts = [4, 2, 3];
const targets = [0, ...masses.flatMap((m, j) => Array.from({ length: counts[j] }, (_, k) => levels[j] + (k + 1) * m / counts[j]))];
const uniform = () => Array.from({ length: 5 }, (_, i) => levels.map((psi, j) => {
  const x = .23 * i + .035 * j * j + .015 * i * j; return { x, y: psi + .3 * x };
}));

test('mass refinement exactly preserves uniform flow on skew curved crosslines, all parent points and input arrays', () => {
  const nodes = uniform(), before = structuredClone(nodes), r = refine(nodes, masses, counts);
  assert.deepEqual(nodes, before); assert.equal(r.diagnostics.affineFallbacks, 0);
  assert.equal(r.diagnostics.limitedEndpointSlopes, 0);
  for (let i = 0; i < nodes.length; i++) {
    let index = 0; assert.deepEqual(r.nodes[i][index], nodes[i][0]);
    for (let j = 0; j < masses.length; j++) {
      for (let k = 1; k <= counts[j]; k++) {
        const p = r.nodes[i][++index], t = k / counts[j], a = nodes[i][j], b = nodes[i][j + 1];
        assert.ok(Math.abs(p.x - (a.x + t * (b.x - a.x))) < 2e-12);
        assert.ok(Math.abs(p.y - .3 * p.x - targets[index]) < 2e-13);
      }
      assert.deepEqual(r.nodes[i][index], nodes[i][j + 1]);
    }
  }
});

test('physical quadratic stagnation streamfunction is preserved on every parent station, with square-root wall-tube subdivision', () => {
  const phi = [-.1, -.025, 0, .025, .1];
  const grid = phi.map(f => levels.map(psi => { const radius = Math.hypot(f, psi); return { x: Math.sqrt(radius + f), y: Math.sqrt(radius - f) }; }));
  const r = refine(grid, masses, counts, { lowerStagnation: 2 });
  assert.equal(r.diagnostics.affineFallbacks, 0); assert.equal(r.diagnostics.limitedEndpointSlopes, 0);
  r.nodes.forEach(row => row.forEach((p, j) => assert.ok(Math.abs(p.x * p.y - targets[j]) < 2e-13)));
  assert.ok(Math.abs(r.nodes[2][1].x / grid[2][1].x - .5) < 2e-13);
  assert.ok(Math.abs(r.nodes[2][1].y / grid[2][1].y - .5) < 2e-13);
  // The opposite bank has the same square-root law with reversed mass order.
  const reverse = refine(grid.map(row => row.slice().reverse()), masses.slice().reverse(), counts.slice().reverse(), { upperStagnation: 2 });
  r.nodes.forEach((row, i) => row.forEach((p, j) => {
    const q = reverse.nodes[i].at(-1 - j); assert.ok(Math.hypot(q.x - p.x, q.y - p.y) < 2e-12);
  }));
});

test('mass refinement is covariant under rotation, translation and independent geometry/mass scaling', () => {
  const grid = uniform(), c = Math.cos(.73), s = Math.sin(.73);
  const transform = p => ({ x: 3 * (c * p.x - s * p.y) + 7, y: 3 * (s * p.x + c * p.y) - 2 });
  const a = refine(grid, masses, counts), b = refine(grid.map(row => row.map(transform)), masses.map(m => 13 * m), counts);
  a.nodes.forEach((row, i) => row.forEach((p, j) => {
    const q = transform(p); assert.ok(Math.hypot(q.x - b.nodes[i][j].x, q.y - b.nodes[i][j].y) < 3e-12);
  }));
});

test('mass interpolation rejects bad inputs and unresolved stencils without moving parent points', () => {
  const grid = uniform(), before = structuredClone(grid);
  for (const [m, n, options] of [[[], counts, {}], [masses, [0, 1, 1], {}], [[0, .05, .08], counts, {}], [masses, counts, { lowerStagnation: -1 }]])
    assert.throws(() => refine(grid, m, n, options), /Invalid/);
  assert.deepEqual(grid, before);
  const line = grid.map(() => levels.map(psi => ({ x: 0, y: psi })));
  assert.throws(() => refine(line, masses, counts), /Unresolved/);
});
