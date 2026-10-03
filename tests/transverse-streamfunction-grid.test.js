import test from 'node:test';
import assert from 'node:assert/strict';
import { createTransverseStreamfunctionGrid, smoothTransverseStreamfunctionGrid } from '../src/geometry/transverse-streamfunction-grid.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { potentialGridQuality } from '../src/geometry/potential-plane-grid.js';

const b = .12;
const clusteredX = u => u + .6 * u * (1 - u);
// Independent analytic harmonic field, Delta psi = -2b + 2b = 0.
const psi = p => p.y + b * (p.x * p.x - p.y * p.y);
// Rationalization avoids cancellation in the inverse near psi=b*x^2.
const inverseY = (x, eta) => 2 * (eta - b * x * x) / (1 + Math.sqrt(1 - 4 * b * (eta - b * x * x)));
const boundary = (i, j, nx, nt) => !i || !j || i === nx || j === nt;
const vertical = nodes => nodes.map(row => row.map(() => ({ x: 0, y: 1 })));
const masses = eta => eta.slice(1).map((v, j) => v - eta[j]);
const fixed = (actual, initial) => actual.forEach((row, i) => row.forEach((p, j) => {
  assert.equal(p.x, initial[i][j].x, 'Transverse moves must preserve every prescribed x station.');
  if (boundary(i, j, actual.length - 1, row.length - 1)) assert.deepEqual(p, initial[i][j]);
}));

function harmonicFixture(nx, nt = nx) {
  const eta = Array.from({ length: nt + 1 }, (_, j) => j / nt);
  const exact = Array.from({ length: nx + 1 }, (_, i) => eta.map(v => {
    const x = clusteredX(i / nx); return { x, y: inverseY(x, v) };
  }));
  const nodes = exact.map((row, i) => row.map((p, j) => boundary(i, j, nx, nt) ? { ...p } : {
    x: p.x, y: (1 - eta[j]) * row[0].y + eta[j] * row[nt].y
      + .01 * Math.sin(Math.PI * i / nx) * Math.sin(Math.PI * eta[j]),
  }));
  return { eta, exact, nodes, massFlows: masses(eta), directions: vertical(nodes) };
}

test('transverse smoothing preserves exact straight streamlines with clustered x and nonuniform mass spacing', () => {
  const eta = [0, .03, .11, .29, .58, .84, 1], nx = 12;
  const nodes = Array.from({ length: nx + 1 }, (_, i) => eta.map(y => ({ x: clusteredX(i / nx), y })));
  const before = structuredClone(nodes);
  const system = createTransverseStreamfunctionGrid({ nodes, massFlows: masses(eta), directions: vertical(nodes) });
  const value = system.evaluate(system.initial, { linearize: true });
  assert.ok(Math.max(...value.gradient.map(Math.abs)) < 1e-11);
  const result = smoothTransverseStreamfunctionGrid(system);
  assert.equal(result.converged, true, result.reason);
  assert.deepEqual(result.nodes, nodes, 'Already harmonic geometry must remain exactly unchanged.');
  assert.deepEqual(nodes, before); assert.equal(system.quality(result.nodes).valid, true);
});

test('transverse streamfunction energy has the analytic gradient and a symmetric positive matrix', () => {
  const fixture = harmonicFixture(7, 5);
  const system = createTransverseStreamfunctionGrid(fixture);
  const { energy, gradient, matrix } = system.evaluate(system.initial, { linearize: true });
  assert.ok(Number.isFinite(energy)); assert.equal(gradient.length, system.n);
  const direction = Float64Array.from(gradient, (_, k) => Math.sin(1.7 * (k + 1)));
  const exactDerivative = gradient.reduce((sum, value, k) => sum + value * direction[k], 0);
  let previous = Infinity;
  for (const h of [1e-4, 5e-5, 2.5e-5]) {
    const plus = system.evaluate(system.move(system.initial, direction, h)).energy;
    const minus = system.evaluate(system.move(system.initial, direction, -h)).energy;
    const error = Math.abs((plus - minus) / (2 * h) - exactDerivative);
    assert.ok(error < .35 * previous, `Energy-gradient error ${error} after ${previous}`); previous = error;
  }
  assert.ok(previous < 1e-6);
  const dense = sparseDense(matrix);
  for (let i = 0; i < system.n; i++) for (let j = 0; j < system.n; j++)
    assert.ok(Math.abs(dense[i * system.n + j] - dense[j * system.n + i]) < 1e-11);
  assert.ok(sparseProduct(matrix, direction).reduce((sum, value, k) => sum + value * direction[k], 0) > 0);
  const expected = sparseProduct(matrix, direction); let previousHessianError = Infinity;
  for (const h of [1e-4, 5e-5, 2.5e-5]) {
    const plus = system.evaluate(system.move(system.initial, direction, h), { linearize: true }).gradient;
    const minus = system.evaluate(system.move(system.initial, direction, -h), { linearize: true }).gradient;
    const error = Math.max(...expected.map((v, k) => Math.abs((plus[k] - minus[k]) / (2 * h) - v)));
    assert.ok(error < .35 * previousHessianError, `Exact-Hessian error ${error} after ${previousHessianError}`);
    previousHessianError = error;
  }
  assert.ok(previousHessianError < 1e-5);
});

test('the exact transverse Hessian equals independently integrated Laplace stiffness on an affine flow', () => {
  const nx = 4, eta = [0, .12, .35, .7, 1], nt = eta.length - 1;
  const nodes = Array.from({ length: nx + 1 }, (_, i) => eta.map(y => ({ x: clusteredX(i / nx), y })));
  const system = createTransverseStreamfunctionGrid({ nodes, massFlows: masses(eta), directions: vertical(nodes) });
  const { matrix } = system.evaluate(system.initial, { linearize: true }), n = system.n;
  // Discover the public move() indexing instead of duplicating an internal
  // flattening convention in this independent element assembly.
  const index = new Map();
  for (let k = 0; k < n; k++) {
    const unit = new Float64Array(n); unit[k] = 1;
    const moved = system.move(system.initial, unit, 1e-4); let changed = 0;
    moved.forEach((row, i) => row.forEach((p, j) => {
      if (Math.abs(p.y - nodes[i][j].y) > 1e-8) { index.set(`${i},${j}`, k); changed++; }
    }));
    assert.equal(changed, 1);
  }
  // Exact 1D integrals of the linear basis: stiffness D and mass M.
  // Tensor products give K=(hy/hx) D⊗M + (hx/hy) M⊗D on each rectangle.
  // This is independent of the source quadrature and energy differentiation.
  const D = [[1, -1], [-1, 1]], M = [[1 / 3, 1 / 6], [1 / 6, 1 / 3]];
  const expected = new Float64Array(n * n), corners = [[0, 0], [1, 0], [1, 1], [0, 1]];
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const hx = nodes[i + 1][j].x - nodes[i][j].x, hy = eta[j + 1] - eta[j];
    for (const [a, b] of corners) for (const [c, d] of corners) {
      const row = index.get(`${i + a},${j + b}`), col = index.get(`${i + c},${j + d}`);
      if (row !== undefined && col !== undefined)
        expected[row * n + col] += hy / hx * D[a][c] * M[b][d] + hx / hy * M[a][c] * D[b][d];
    }
  }
  const actual = sparseDense(matrix);
  actual.forEach((v, k) => assert.ok(Math.abs(v - expected[k]) < 5e-12, `Stiffness entry ${k}: ${v} != ${expected[k]}`));
});

test('transverse relaxation from perturbed seeds refines to a nonconformal exact harmonic streamfunction', t => {
  let previous = Infinity; const evidence = [];
  for (const n of [4, 8, 16]) {
    const fixture = harmonicFixture(n), before = structuredClone(fixture.nodes);
    const system = createTransverseStreamfunctionGrid(fixture);
    const result = smoothTransverseStreamfunctionGrid(system);
    assert.equal(result.converged, true, JSON.stringify({ n, reason: result.reason, last: result.history.at(-1) }));
    assert.equal(system.quality(result.nodes).valid, true); fixed(result.nodes, fixture.nodes);
    assert.deepEqual(fixture.nodes, before);
    let error = 0, initialError = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      error = Math.max(error, Math.abs(psi(p) - fixture.eta[j]));
      initialError = Math.max(initialError, Math.abs(psi(fixture.nodes[i][j]) - fixture.eta[j]));
    }));
    assert.ok(error < initialError / 10, `Solved psi error ${error} from seed error ${initialError}`);
    assert.ok(error < .35 * previous, `Refinement psi error ${error} after ${previous}`); previous = error;
    evidence.push({ n, error, initialError, iterations: result.history.length - 1,
      minCornerSine: system.quality(result.nodes).minCornerSine });
  }
  assert.ok(previous < 2e-4); t.diagnostic(JSON.stringify(evidence));
});

test('transverse moves keep physical boundaries and do not mutate the original grid', () => {
  const fixture = harmonicFixture(5, 4), before = structuredClone(fixture.nodes);
  const system = createTransverseStreamfunctionGrid(fixture);
  const delta = Float64Array.from({ length: system.n }, (_, k) => .002 * Math.cos(k));
  const moved = system.move(system.initial, delta, .5);
  fixed(moved, fixture.nodes); assert.deepEqual(system.initial, before); assert.deepEqual(fixture.nodes, before);
  assert.notDeepEqual(moved, before);
});

test('varying transverse directions retain the exact Hessian and coordinate invariance', () => {
  const fixture = harmonicFixture(6, 5);
  fixture.directions = fixture.nodes.map((row, i) => row.map((p, j) => {
    const angle = .55 * Math.sin(1.3 * i + .7 * j); return { x: Math.sin(angle), y: Math.cos(angle) };
  }));
  const system = createTransverseStreamfunctionGrid(fixture), e = system.evaluate(system.initial, { linearize: true });
  const direction = Float64Array.from(e.gradient, (_, k) => Math.cos(.8 * k));
  const expected = sparseProduct(e.matrix, direction); let previous = Infinity;
  for (const h of [6e-5, 3e-5, 1.5e-5]) {
    const plus = system.evaluate(system.move(system.initial, direction, h), { linearize: true }).gradient;
    const minus = system.evaluate(system.move(system.initial, direction, -h), { linearize: true }).gradient;
    const error = Math.max(...expected.map((v, k) => Math.abs((plus[k] - minus[k]) / (2 * h) - v)));
    assert.ok(error < .35 * previous, `${error} after ${previous}`); previous = error;
  }
  assert.ok(previous < 1e-5);
  const angle = .71, c = Math.cos(angle), s = Math.sin(angle), scale = 2.3;
  const rotate = p => ({ x: c * p.x - s * p.y, y: s * p.x + c * p.y });
  const mapped = createTransverseStreamfunctionGrid({
    nodes: fixture.nodes.map(row => row.map(p => { const q = rotate(p); return { x: 3 + scale * q.x, y: -4 + scale * q.y }; })),
    directions: fixture.directions.map(row => row.map(rotate)), massFlows: fixture.massFlows.map(m => 7 * m),
  });
  const m = mapped.evaluate(mapped.initial, { linearize: true });
  assert.ok(Math.abs(e.energy - m.energy) < 1e-12);
  e.gradient.forEach((v, k) => assert.ok(Math.abs(v - scale * m.gradient[k]) < 1e-11));
  e.matrix.values.forEach((v, k) => assert.ok(Math.abs(v - scale ** 2 * m.matrix.values[k]) < 1e-10));
});

test('scalar SLOR and a separately certified sparse reference find the same transverse grid', () => {
  const fixture = harmonicFixture(8), system = createTransverseStreamfunctionGrid(fixture);
  const slor = smoothTransverseStreamfunctionGrid(system);
  const reference = smoothTransverseStreamfunctionGrid(system, { linearSolve: solveSparseDirect, linearSolverName: 'KLU reference' });
  assert.equal(slor.converged, true); assert.equal(reference.converged, true);
  assert.equal(slor.linearBackend, 'scalar streamwise SLOR'); assert.equal(reference.linearBackend, 'KLU reference');
  slor.nodes.forEach((row, i) => row.forEach((p, j) => assert.ok(Math.hypot(p.x - reference.nodes[i][j].x, p.y - reference.nodes[i][j].y) < 1e-9)));
});

test('transverse chart violations and failed or uncertified solves cannot pass', () => {
  const fixture = harmonicFixture(5), system = createTransverseStreamfunctionGrid(fixture);
  const stopped = smoothTransverseStreamfunctionGrid(system, { maxIterations: 0 });
  assert.equal(stopped.converged, false); assert.equal(stopped.reason, 'iteration limit');
  const uncertified = smoothTransverseStreamfunctionGrid(system, {
    linearSolve: () => ({ x: new Float64Array(system.n), converged: true, relativeResidual: 0 }),
  });
  assert.equal(uncertified.converged, false); assert.equal(uncertified.reason, 'linear solve failed');
  const changed = structuredClone(system.initial); changed[2][2].x += .001;
  assert.throws(() => system.evaluate(changed), /guide line/);
  const movedBoundary = structuredClone(system.initial); movedBoundary[0][2].y += .001;
  assert.throws(() => system.evaluate(movedBoundary), /boundaries/);
  const straight = Array.from({ length: 4 }, (_, i) => Array.from({ length: 4 }, (_, j) => ({ x: i / 3, y: j / 3 })));
  assert.throws(() => createTransverseStreamfunctionGrid({ nodes: straight, massFlows: [1, 1, 1],
    directions: straight.map(row => row.map(() => ({ x: 1, y: 0 }))) }), /transverse guides/);
  const folded = structuredClone(straight); folded[1][1].y = 1.5;
  assert.throws(() => createTransverseStreamfunctionGrid({ nodes: folded, massFlows: [1, 1, 1], directions: vertical(straight) }), /positive cells/);
  assert.throws(() => system.move(system.initial, [0], 1), /direction/);
  assert.throws(() => smoothTransverseStreamfunctionGrid(system, { maxIterations: -1 }), /controls/);
});

test('positive cells and transverse nodal directions do not hide a tangent guide inside a cell', () => {
  const a = [{ x: 1, y: 0 }, { x: .5, y: Math.sqrt(3) / 2 }, { x: .5, y: Math.sqrt(3) / 2 }];
  const nodes = Array.from({ length: 3 }, (_, i) => a.map((t, j) => ({ x: i / 2 * t.x, y: 2 * j + i / 2 * t.y })));
  const dirs = [10, 230, 230].map(degrees => ({ x: Math.cos(degrees * Math.PI / 180), y: Math.sin(degrees * Math.PI / 180) }));
  assert.equal(potentialGridQuality(nodes).valid, true);
  a.forEach((t, j) => assert.ok(t.x * dirs[j].y - t.y * dirs[j].x > .17));
  // Inside the lower row of cells, the tangent and guide turn in opposite
  // directions. At t=1/2 their cross product is negative despite all six
  // corner tests being positive. The quadratic-in-t minimum must catch it.
  const t = { x: (a[0].x + a[1].x) / 2, y: (a[0].y + a[1].y) / 2 };
  const d = { x: (dirs[0].x + dirs[1].x) / 2, y: (dirs[0].y + dirs[1].y) / 2 };
  assert.ok(t.x * d.y - t.y * d.x < 0);
  assert.throws(() => createTransverseStreamfunctionGrid({ nodes, massFlows: [1, 1],
    directions: nodes.map(() => structuredClone(dirs)) }), /transverse guides/);
});
