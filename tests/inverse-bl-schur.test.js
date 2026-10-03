import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { prepareCoupledInverseBL } from '../scripts/validation/coupled-inverse-bl.js';
import { reduceInverseBL } from '../scripts/validation/inverse-bl-schur.js';
import { solveLinear } from '../src/numerics/linear.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';

test('inverse BL Schur solution satisfies the original simultaneous two-element Euler/BL matrix', t => {
  const s = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    { transitionMode: 'automatic', edgeMatching: 'section-velocity', reynolds: 1e6, ncrit: 9 });
  const p = prepareCoupledInverseBL(s, s.initial); s.bl.restoreActive(p.phase);
  const matrix = s.jacobian(p.x), rhs = s.residual(p.x).map(v => -v), reduced = reduceInverseBL(s, matrix, rhs);
  const reference = solveLinear(sparseDense(matrix), rhs), y = solveLinear(sparseDense(reduced.matrix), reduced.rhs), x = reduced.lift(y);
  const product = sparseProduct(matrix, x); let error = 0, residual = 0;
  x.forEach((v, i) => { error = Math.max(error, Math.abs(v - reference[i]) / Math.max(1, Math.abs(v), Math.abs(reference[i]))); });
  product.forEach((v, i) => { residual = Math.max(residual, Math.abs(v - rhs[i])); });
  assert.ok(error < 1e-10 && residual < 1e-10);
  assert.equal(reduced.matrix.n, s.ne + s.bl.stations.length);
  assert.equal(reduced.rows.length, reduced.columns.length);
  // A separate arbitrary reduced vector verifies the complete operator,
  // including inhomogeneous BL rows, without relying on either LU result.
  const trial = Float64Array.from(y, (_, i) => .01 * Math.sin(i * .37)), lifted = reduced.lift(trial);
  const full = sparseProduct(matrix, lifted).map((v, i) => v - rhs[i]);
  const small = sparseProduct(reduced.matrix, trial).map((v, i) => v - reduced.rhs[i]);
  let operatorError = 0;
  reduced.rows.forEach((row, i) => { operatorError = Math.max(operatorError, Math.abs(full[row] - small[i])); });
  assert.ok(operatorError < 1e-10);
  t.diagnostic(JSON.stringify({ original: s.n, reduced: reduced.matrix.n, nonzeros: reduced.matrix.values.length, error, residual, operatorError }));
});
