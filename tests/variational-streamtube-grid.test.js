import test from 'node:test';
import assert from 'node:assert/strict';
import { createVariationalStreamtubeGrid, smoothVariationalStreamtubeGrid } from '../src/geometry/tests/variational-streamtube-grid.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { solveLinear } from '../src/numerics/linear.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { factorStreamwiseBlockLines, solveStreamwiseBlockLines } from '../src/numerics/tests/streamwise-block-lines.js';

const rectangle = (nx, nt) => Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => ({ x: -1 + 2 * i / nx, y: j / nt })));
const fixed = (a, b) => a.forEach((row, i) => row.forEach((p, j) => {
  if (!i || !j || i === a.length - 1 || j === row.length - 1) assert.deepEqual(p, b[i][j]);
}));
// Inverse of the nonconformal harmonic coordinates u=x, v=y+.12(x^2-y^2).
// Harmonic x(u,v), y(u,v) would be the wrong continuous problem here.
const inverse = p => ({ x: p.x, y: 2 * (p.y - .12 * p.x * p.x) / (1 + Math.sqrt(1 - .48 * (p.y - .12 * p.x * p.x))) });
const harmonic = p => ({ x: p.x, y: p.y + .12 * (p.x * p.x - p.y * p.y) });

function specimen(nx = 7, nt = 5) {
  const coordinates = rectangle(nx, nt), nodes = coordinates.map(row => row.map(inverse));
  for (let i = 1; i < nx; i++) for (let j = 1; j < nt; j++) {
    nodes[i][j].x += .006 * Math.sin(i + j); nodes[i][j].y -= .004 * Math.cos(i - j);
  }
  return createVariationalStreamtubeGrid({ nodes, coordinates });
}

test('Winslow functional has the analytic energy gradient and a positive symmetric Gauss-Newton matrix', () => {
  const system = specimen(), { energy, gradient, matrix } = system.evaluate(system.initial, { linearize: true });
  assert.ok(energy > 0);
  const direction = Float64Array.from(gradient, (_, k) => Math.sin(1.7 * k)), exact = gradient.reduce((s, v, k) => s + v * direction[k], 0);
  let previous = Infinity;
  for (const h of [1e-4, 5e-5, 2.5e-5]) {
    const fd = (system.evaluate(system.move(system.initial, direction, h)).energy - system.evaluate(system.move(system.initial, direction, -h)).energy) / (2 * h);
    const error = Math.abs(fd - exact); assert.ok(error < .3 * previous, `${error} after ${previous}`); previous = error;
  }
  const dense = sparseDense(matrix), n = system.n;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) assert.ok(Math.abs(dense[i * n + j] - dense[j * n + i]) < 1e-12);
  assert.ok(sparseProduct(matrix, direction).reduce((s, v, k) => s + v * direction[k], 0) > 0);
});

test('two-coordinate block Thomas SLOR matches dense line solves and the certified global solution', () => {
  const system = specimen(6, 4), { gradient, matrix } = system.evaluate(system.initial, { linearize: true });
  const rhs = Float64Array.from(gradient, v => -v), lines = factorStreamwiseBlockLines(matrix, system), omega = 1.25;
  const dense = sparseDense(matrix), n = system.n, reference = new Float64Array(n);
  for (let j = 1; j < system.nt; j++) {
    const ids = [];
    for (let i = 1; i < system.nx; i++) ids.push(2 * ((i - 1) * (system.nt - 1) + j - 1), 2 * ((i - 1) * (system.nt - 1) + j - 1) + 1);
    const a = ids.flatMap(r => ids.map(c => dense[r * n + c]));
    const b = ids.map(r => rhs[r] - reference.reduce((s, v, c) => s + dense[r * n + c] * v, 0));
    const d = solveLinear(a, b); ids.forEach((id, k) => { reference[id] += omega * d[k]; });
  }
  const result = lines.sweep(new Float64Array(n), rhs, omega);
  result.forEach((v, k) => assert.ok(Math.abs(v - reference[k]) < 1e-13));
  const direct = solveSparseDirect(matrix, rhs), slor = solveStreamwiseBlockLines(matrix, rhs, system);
  assert.equal(slor.converged, true); assert.ok(slor.relativeResidual < 1e-10);
  slor.x.forEach((v, k) => assert.ok(Math.abs(v - direct.x[k]) < 1e-10));
  assert.equal(solveStreamwiseBlockLines(matrix, rhs, system, { maxSweeps: 0 }).converged, false);
});

test('variational grid refinement converges to nonconformal harmonic forward coordinates', t => {
  let previous = Infinity; const evidence = [];
  for (const n of [4, 8, 16]) {
    const coordinates = rectangle(n, n), nodes = coordinates.map(row => row.map(inverse));
    const system = createVariationalStreamtubeGrid({ nodes, coordinates }), result = smoothVariationalStreamtubeGrid(system);
    assert.equal(result.converged, true, JSON.stringify({ reason: result.reason, last: result.history.at(-1) }));
    assert.equal(result.linearBackend, 'streamwise block SLOR');
    assert.equal(result.physicsValidated, false); assert.equal(result.quality.valid, true); fixed(result.nodes, nodes);
    let error = 0;
    result.nodes.forEach((row, i) => row.forEach((p, j) => {
      const w = harmonic(p); error = Math.max(error, Math.hypot(w.x - coordinates[i][j].x, w.y - coordinates[i][j].y));
    }));
    assert.ok(error < .3 * previous, `${error} after ${previous}`); previous = error;
    evidence.push({ n, error, iterations: result.history.length - 1, minCornerSine: result.quality.minCornerSine });
  }
  t.diagnostic(JSON.stringify(evidence));
});

test('nonconformal affine grids are stationary, transform invariant and keep nonzero conformal defect', () => {
  const coordinates = rectangle(6, 5), affine = p => ({ x: 2 + p.x + .4 * p.y, y: -1 + .2 * p.x + 1.8 * p.y });
  const nodes = coordinates.map(row => row.map(affine)), before = structuredClone(nodes), r = smoothVariationalStreamtubeGrid(createVariationalStreamtubeGrid({ nodes, coordinates }));
  assert.equal(r.converged, true); assert.equal(r.history.length, 1); assert.ok(r.history[0].energy > .1); assert.deepEqual(r.nodes, nodes);
  const c = Math.cos(.4), s = Math.sin(.4), physical = p => ({ x: 3 + 2 * (c * p.x - s * p.y), y: -4 + 2 * (s * p.x + c * p.y) });
  const transformed = createVariationalStreamtubeGrid({ nodes: nodes.map(row => row.map(physical)), coordinates: coordinates.map(row => row.map(p => ({ x: 8 + 3 * p.x, y: -2 + 3 * p.y }))) });
  const q = smoothVariationalStreamtubeGrid(transformed); assert.equal(q.converged, true);
  assert.ok(Math.abs(q.history[0].energy / r.history[0].energy - 9) < 1e-11); assert.deepEqual(nodes, before);
});

test('variational relaxation rejects infeasible grids and distinguishes a limit from convergence', () => {
  const system = specimen(), r = smoothVariationalStreamtubeGrid(system, { maxIterations: 0 });
  assert.equal(r.converged, false); assert.equal(r.reason, 'iteration limit'); assert.deepEqual(r.nodes, system.initial);
  const wrong = structuredClone(system.initial); wrong[0][1].x += .1;
  assert.throws(() => system.evaluate(wrong), /boundaries/);
  const folded = structuredClone(system.initial); folded[3][2].x += 10;
  assert.throws(() => createVariationalStreamtubeGrid({ nodes: folded, coordinates: system.coordinates }), /positive/);
  assert.throws(() => createVariationalStreamtubeGrid({ nodes: system.initial, coordinates: [] }), /dimensions/);
  const failed = smoothVariationalStreamtubeGrid(system, { linearSolve: () => ({ converged: false, relativeResidual: 1 }) });
  assert.equal(failed.converged, false); assert.equal(failed.reason, 'linear solve failed'); assert.deepEqual(failed.nodes, system.initial);
  const uncertified = smoothVariationalStreamtubeGrid(system, { linearSolve: () => ({ x: new Float64Array(system.n), converged: true, relativeResidual: 0 }) });
  assert.equal(uncertified.converged, false); assert.equal(uncertified.reason, 'linear solve failed');
  assert.equal(uncertified.history[0].linearResidual, 1);
});

test('variational refinement approaches analytic cylinder flow and preserves both boundary stagnation points', t => {
  const sqrt = (x, y) => ({ x: Math.sqrt(Math.max(0, (Math.hypot(x, y) + x) / 2)), y: Math.sqrt(Math.max(0, (Math.hypot(x, y) - x) / 2)) });
  const cylinder = w => { const a = sqrt(w.x - 2, w.y), b = sqrt(w.x + 2, w.y);
    return { x: .5 * (w.x + a.x * b.x - a.y * b.y), y: .5 * (w.y + a.x * b.y + a.y * b.x) }; };
  let previous = Infinity; const evidence = [];
  for (const [nx, nt] of [[12, 6], [24, 12], [48, 24]]) {
    const coordinates = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => ({ x: -3 + 6 * i / nx, y: 1.4 * j / nt })));
    const nodes = coordinates.map(row => row.map(cylinder)), r = smoothVariationalStreamtubeGrid(createVariationalStreamtubeGrid({ nodes, coordinates }));
    assert.equal(r.converged, true, JSON.stringify({ reason: r.reason, last: r.history.at(-1) })); fixed(r.nodes, nodes);
    let error = 0;
    r.nodes.forEach((row, i) => row.forEach((p, j) => {
      const r2 = p.x * p.x + p.y * p.y, w = coordinates[i][j];
      assert.ok(r2 >= 1 - 2e-12, 'Interior grid nodes cannot enter the cylinder');
      error = Math.max(error, Math.hypot(p.x + p.x / r2 - w.x, p.y - p.y / r2 - w.y));
    }));
    assert.ok(error < .85 * previous, `${error} after ${previous}`); previous = error;
    evidence.push({ nx, nt, error, iterations: r.history.length - 1, minCornerSine: r.quality.minCornerSine });
  }
  t.diagnostic(JSON.stringify(evidence));
});
