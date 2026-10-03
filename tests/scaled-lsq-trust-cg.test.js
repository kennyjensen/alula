import test from 'node:test';
import assert from 'node:assert/strict';
import { scaledLeastSquaresTrustCG } from '../scripts/validation/scaled-lsq-trust-cg.js';
import { solveLinear } from '../src/numerics/linear.js';
import { sparseMatrix, sparseAdd } from '../src/numerics/sparse.js';

const dense = Float64Array.from([2, 1, 0, 0, 3, 1, 1, 0, 2]);
const make = (units = [1, 1, 1]) => {
  const a = sparseMatrix([[0, 1, 2], [0, 1, 2], [0, 1, 2]]);
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) sparseAdd(a, i, j, dense[3 * i + j] / units[j]);
  return a;
};
const residual = Float64Array.from([1, -.4, 2]);

test('interior trust CG matches an independent dense Newton solve and certifies the actual model gradient', () => {
  const r = scaledLeastSquaresTrustCG(make(), residual, { radius: 100, gradientTolerance: 1e-12, maxIterations: 10 });
  const expected = solveLinear(dense, residual.map(v => -v));
  assert.ok(r.gradientConverged && r.relativeGradientResidual < 1e-12 && r.linearRelativeResidual < 1e-12);
  r.direction.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 1e-12));
});

test('boundary trust CG agrees with the analytic scaled-gradient intersection and is invariant to variable units', () => {
  const radius = .1, a = make(), r = scaledLeastSquaresTrustCG(a, residual, { radius });
  const gradient = Float64Array.from([4, -.2, 3.6], (v, i) => v / r.scales[i]);
  const gn = Math.hypot(...gradient);
  assert.equal(r.reason, 'trust boundary'); assert.equal(r.iterations, 1);
  r.direction.forEach((v, i) => assert.ok(Math.abs(v + radius * gradient[i] / gn / r.scales[i]) < 1e-12));
  const units = [1e3, 1e-3, 5], q = scaledLeastSquaresTrustCG(make(units), residual, { radius });
  q.direction.forEach((v, i) => assert.ok(Math.abs(v / units[i] - r.direction[i]) < 1e-12));
  assert.ok(Math.abs(r.scaledStepNorm - radius) < 1e-12);
});

test('a bounded Krylov proposal records unresolved Newton error while improving on its first Cauchy iterate', () => {
  const first = scaledLeastSquaresTrustCG(make(), residual, { radius: 100, maxIterations: 1, gradientTolerance: 1e-12 });
  const second = scaledLeastSquaresTrustCG(make(), residual, { radius: 100, maxIterations: 2, gradientTolerance: 1e-12 });
  assert.equal(second.reason, 'iteration limit'); assert.equal(second.gradientConverged, false);
  assert.ok(second.linearRelativeResidual > 1e-10);
  assert.ok(second.predictedReduction > first.predictedReduction);
});
