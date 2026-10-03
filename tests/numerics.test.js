import test from 'node:test';
import assert from 'node:assert/strict';
import { solveLinear, linearResidual, normInf } from '../src/numerics/linear.js';
import { solveNewton } from '../src/numerics/newton.js';

test('pivoted dense solve respects row scaling and preserves inputs', () => {
  const a = [0, 2e-14, 3e-14, 4e12, 5e12, 6e12, 7, 8, 10];
  const b = [13e-14, 32e12, 53];
  const x = solveLinear(a, b);
  assert.ok(normInf(x.map((v, i) => v - (i + 1))) < 1e-12);
  assert.equal(a[0], 0); assert.equal(b[0], 13e-14);
  assert.ok(normInf(linearResidual([2, 1, 1, 3], solveLinear([2, 1, 1, 3], [1, 2]), [1, 2])) < 1e-14);
});
test('singular and nonfinite linear systems fail explicitly', () => {
  assert.throws(() => solveLinear([1, 2, 2, 4], [1, 2]), /Singular/);
  assert.throws(() => solveLinear([NaN], [1]), /nonfinite/);
});
test('global Newton solves a coupled nonlinear system with damping', () => {
  const result = solveNewton({ initial: [0.2, 0.2], residual: ([x, y]) => [x * x + y - 5, x + y * y - 7],
    admissible: x => x.every(v => v > 0) });
  assert.equal(result.converged, true);
  assert.ok(result.history.at(-1).residual < 1e-10);
  assert.ok(result.history.some(h => h.step > 0 && h.step < 1));
  for (let i = 1; i < result.history.length; i++) assert.ok(result.history[i].residual < result.history[i - 1].residual);
});
test('Newton does not misreport iteration limits or singular Jacobians as convergence', () => {
  assert.equal(solveNewton({ initial: [1], residual: ([x]) => [x * x - 2], maxIterations: 0 }).converged, false);
  assert.equal(solveNewton({ initial: [0], residual: () => [1] }).converged, false);
});

test('Newton observers receive isolated accepted states, including a retained iteration-limit state', () => {
  const options = { initial: [.1], residual: ([x]) => [x * x - 2], maxIterations: 3 };
  const baseline = solveNewton(options); const frames = [];
  const observed = solveNewton({ ...options, onState: state => {
    frames.push(structuredClone(state));
    state.x.fill(NaN); state.iteration.residual = NaN;
  } });
  assert.deepEqual(observed, baseline);
  assert.equal(observed.reason, 'iteration limit');
  assert.equal(frames.length, observed.history.length);
  assert.ok(observed.history.some(h => h.step > 0 && h.step < 1));
  frames.forEach(({ x, iteration }, i) => {
    assert.deepEqual(iteration, observed.history[i]);
    assert.equal(iteration.residual, Math.abs(options.residual(x)[0]));
  });
  assert.deepEqual(frames.at(-1).x, observed.x);
});
