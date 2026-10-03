import test from 'node:test';
import assert from 'node:assert/strict';
import { createPotentialPlaneGrid, smoothPotentialPlaneGrid } from '../src/geometry/potential-plane-grid.js';
import { solveLinear } from '../src/numerics/linear.js';

const plane = (nx, nt) => Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
  const u = i / nx, v = Math.expm1(j / nt) / Math.expm1(1);
  return { x: u + .2 * v + .06 * Math.sin(Math.PI * u) * v, y: v };
}));
const boundary = (nodes, expected) => nodes.forEach((row, i) => row.forEach((p, j) => {
  if (!i || i === nodes.length - 1 || !j || j === row.length - 1) assert.deepEqual(p, expected[i][j]);
}));
const affine = w => ({ x: 2 + 1.4 * w.x - .3 * w.y, y: -1 + .3 * w.x + 1.4 * w.y });

test('potential-plane SLOR exactly recovers an affine conformal map on a skewed nonuniform computational grid', () => {
  const coordinates = plane(12, 7), exact = coordinates.map(row => row.map(affine));
  const nodes = exact.map((row, i) => row.map((p, j) => !i || i === 12 || !j || j === 7 ? { ...p }
    : { x: p.x + .02 * Math.sin(i + j), y: p.y - .015 * Math.cos(i - j) }));
  const before = structuredClone(nodes), input = structuredClone(coordinates);
  const system = createPotentialPlaneGrid({ nodes, coordinates });
  assert.ok(system.residuals(exact).residual < 1e-14);
  const r = smoothPotentialPlaneGrid(system, { tolerance: 1e-13 });
  assert.equal(r.converged, true, r.reason); assert.equal(r.physicsValidated, false); boundary(r.nodes, nodes);
  r.nodes.forEach((row, i) => row.forEach((p, j) => assert.ok(Math.hypot(p.x - exact[i][j].x, p.y - exact[i][j].y) < 2e-11)));
  assert.deepEqual(nodes, before); assert.deepEqual(coordinates, input);
  const exhausted = smoothPotentialPlaneGrid(system, { maxSweeps: 0 });
  assert.equal(exhausted.converged, false); assert.equal(exhausted.reason, 'sweep limit'); assert.deepEqual(exhausted.nodes, nodes);
  const moved = structuredClone(nodes); moved[0][0].x += .01;
  assert.throws(() => system.sweep(moved), /boundary/);
});

test('a factored potential-plane SLOR sweep matches an independent dense block SOR calculation', () => {
  const nx = 6, nt = 4, coordinates = plane(nx, nt), nodes = coordinates.map((row, i) => row.map((p, j) => {
    const a = affine(p), f = .02 * Math.sin(i * .5 + j); return { x: a.x + f, y: a.y - .3 * f };
  })), system = createPotentialPlaneGrid({ nodes, coordinates }), omega = 1.3;
  const locations = [];
  for (let i = 1; i < nx; i++) for (let j = 1; j < nt; j++) locations.push({ i, j });
  const flat = grid => system.residuals(grid).rows.map(r => r.x), base = flat(nodes), n = base.length, h = 1e-4;
  const matrix = new Float64Array(n * n), delta = new Float64Array(n);
  locations.forEach(({ i, j }, col) => {
    const p = structuredClone(nodes), m = structuredClone(nodes); p[i][j].x += h; m[i][j].x -= h;
    const a = flat(p), b = flat(m);
    for (let row = 0; row < n; row++) matrix[row * n + col] = (a[row] - b[row]) / (2 * h);
  });
  for (let j = 1; j < nt; j++) {
    const indices = locations.flatMap((p, i) => p.j === j ? [i] : []);
    const a = indices.flatMap(row => indices.map(col => matrix[row * n + col]));
    const b = indices.map(row => -base[row] - delta.reduce((s, v, col) => s + matrix[row * n + col] * v, 0));
    const update = solveLinear(a, b); indices.forEach((id, k) => { delta[id] += omega * update[k]; });
  }
  const result = system.sweep(nodes, omega);
  locations.forEach(({ i, j }, k) => assert.ok(Math.abs(result.nodes[i][j].x - nodes[i][j].x - delta[k]) < 2e-11));
  boundary(result.nodes, nodes);
});

const sqrt = (x, y) => {
  const r = Math.hypot(x, y); return { x: Math.sqrt(Math.max(0, (r + x) / 2)), y: Math.sqrt(Math.max(0, (r - x) / 2)) };
};
// Analytic inverse w=z+1/z on the upper half-plane, outside a unit circle.
// Multiplying the two principal square roots fixes the branch at both noses.
function cylinder(w) {
  const a = sqrt(w.x - 2, w.y), b = sqrt(w.x + 2, w.y);
  return { x: .5 * (w.x + a.x * b.x - a.y * b.y), y: .5 * (w.y + a.x * b.y + a.y * b.x) };
}
test('inverse-potential refinement approaches the analytic cylinder flow including boundary stagnation points', t => {
  let previous = Infinity; const evidence = [];
  for (const [nx, nt] of [[12, 6], [24, 12], [48, 24]]) {
    const coordinates = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => ({ x: -3 + 6 * i / nx, y: 1.4 * j / nt })));
    const exact = coordinates.map(row => row.map(cylinder)), system = createPotentialPlaneGrid({ nodes: exact, coordinates });
    const result = smoothPotentialPlaneGrid(system, { tolerance: 1e-12, maxSweeps: 2000 });
    assert.equal(result.converged, true, result.reason); boundary(result.nodes, exact);
    let error = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      // Check the independent analytic forward flow map, not the grid's
      // own discrete Laplace equation or stiffness matrix.
      const r2 = p.x * p.x + p.y * p.y, w = coordinates[i][j];
      assert.ok(r2 >= 1 - 2e-12, 'The inverse cylinder grid must remain outside the body');
      error = Math.max(error, Math.hypot(p.x + p.x / r2 - w.x, p.y - p.y / r2 - w.y));
    }));
    // Boundary square-root singularities preclude an assumed global h^2 rate.
    assert.ok(error < .85 * previous, `Forward-map error ${error} after ${previous}`); previous = error;
    evidence.push({ nx, nt, error, sweeps: result.history.length - 1, minimumCorner: result.quality.minCornerSine });
  }
  t.diagnostic(JSON.stringify(evidence));
});

test('inverse-potential grid equations preserve rigid transforms and potential units and reject folded computational cells', () => {
  const coordinates = plane(8, 5), nodes = coordinates.map(row => row.map(affine)), c = Math.cos(.7), s = Math.sin(.7);
  const physical = p => ({ x: 4 + 3 * (c * p.x - s * p.y), y: -2 + 3 * (s * p.x + c * p.y) });
  const transformed = createPotentialPlaneGrid({ nodes: nodes.map(row => row.map(physical)),
    coordinates: coordinates.map(row => row.map(p => ({ x: 10 + 7 * p.x, y: -3 + 7 * p.y }))) });
  assert.ok(transformed.residuals(transformed.initial).residual < 2e-14);
  const folded = structuredClone(coordinates); folded[3][2].x += 5;
  assert.throws(() => createPotentialPlaneGrid({ nodes, coordinates: folded }), /Invalid quadrilateral/);
  assert.throws(() => createPotentialPlaneGrid({ nodes, coordinates: [] }), /Invalid potential/);
});
