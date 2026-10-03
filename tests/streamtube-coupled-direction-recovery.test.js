import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { coupledNewtonRayStalled, coupledRegularizedDirection, coupledStepStagnation } from '../src/euler/streamtube-coupled-direction-recovery.js';
import { sparseMatrix, sparseAdd, sparseProduct } from '../src/numerics/sparse.js';

test('regularization retains Euler and matching rows and never mutates the governing Jacobian', () => {
  const matrix = sparseMatrix(Array.from({ length: 5 }, (_, i) => [i]));
  for (let i = 0; i < 5; i++) sparseAdd(matrix, i, i, i + 1);
  const before = structuredClone(matrix), residual = Float64Array.of(1, 2, 3, 4, 5);
  const step = coupledRegularizedDirection(matrix, residual, 1, .5);
  assert.deepEqual(matrix, before);
  assert.equal(step.diagnostics.changedDiagonals, 3);
  assert.equal(step.linear.x[0], -1);
  assert.equal(step.linear.x[4], -1);
  assert.ok(step.diagnostics.linearizedMeritSlope < 0);
  assert.ok(step.diagnostics.relativeLinearResidual <= 1e-10);
  assert.equal(step.diagnostics.equationsChanged, false);
});

test('alternate directions retain the Euler/BL cross coupling and matching equations', () => {
  const matrix = sparseMatrix(Array.from({ length: 5 }, () => [0, 1, 2, 3, 4]));
  for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++)
    sparseAdd(matrix, i, j, i === j ? i + 2 : .1 * (i + 1));
  const residual = Float64Array.of(1, 2, 3, 4, 5);
  const candidate = coupledRegularizedDirection(matrix, residual, 1, .5);
  assert.ok(candidate);
  const product = sparseProduct(matrix, candidate.linear.x);
  for (const row of [0, 4]) assert.ok(Math.abs(product[row] + residual[row]) < 1e-12);
  assert.ok(Math.abs(product[2] + residual[2]) > .1);
});

test('direction recovery only interrupts repeated nearly stationary Newton steps', () => {
  const decrease = { beforeSquaredNorm: 1, afterSquaredNorm: .9999 };
  const h = { step: 1e-5, residualDecrease: decrease };
  assert.equal(coupledNewtonRayStalled([h, h], h, { changed: false }, decrease), true);
  assert.equal(coupledNewtonRayStalled([h], h, { changed: false }, decrease), false);
  for (const change of [{ activeChange: true }, { directionRecovery: {} }, { step: 1 },
    { residualDecrease: { beforeSquaredNorm: 1, afterSquaredNorm: .9 } }])
    assert.equal(coupledNewtonRayStalled([h, { ...h, ...change }], h, { changed: false }, decrease), false);
  assert.equal(coupledNewtonRayStalled([h, h], h, { changed: true }, decrease), false);
});

test('the reported stalled RAE state reaches the same equations with guarded alternate directions', () => {
  const cp = JSON.parse(gunzipSync(fs.readFileSync(new URL('./fixtures/rae16x7-stalled-coupled.json.gz', import.meta.url))));
  const before = structuredClone(cp), controls = cp.continuation;
  const legacy = solveCoupledStreamtubeIses(undefined, { ...controls, resume: cp, maxIterations: 1, tolerance: 1e-10 });
  assert.equal(legacy.converged, false);
  assert.equal(legacy.lastRejectedStep.code, 'COUPLED_RESIDUAL_DECREASE');
  const result = solveCoupledStreamtubeIses(undefined, { ...controls, resume: cp,
    directionRecovery: true, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(result.converged, true, result.reason);
  for (const value of Object.values(result.families)) assert.ok(value <= 1e-10);
  assert.ok(result.history.some(h => h.directionRecovery));
  for (const h of result.history.filter(h => h.directionRecovery)) {
    assert.ok(h.directionRecovery.linearizedMeritSlope < 0);
    assert.ok(h.residualDecrease.afterSquaredNorm <= h.residualDecrease.allowedSquaredNorm);
  }
  assert.equal(result.mesh.quality.valid, true);
  assert.deepEqual(result.solverInput.upwind, cp.restart.input.upwind);
  assert.equal(result.solverInput.mach, .74);
  assert.equal(result.conditions.reynolds, cp.restart.options.reynolds);
  assert.equal(result.conditions.ncrit, 9);
  assert.equal(result.checkpoint.continuation.directionRecovery, true);
  assert.deepEqual(cp, before);
});


test('roundoff-sized coupled updates stop without declaring convergence or crossing transition phases', () => {
  const rows = Array.from({ length: 3 }, (_, iteration) => ({ iteration, step: 1e-12, residual: .1,
    residualDecrease: { beforeSquaredNorm: 1, afterSquaredNorm: .9999 } }));
  assert.ok(coupledStepStagnation(rows, 1e-10));
  assert.equal(coupledStepStagnation(rows.slice(1), 1e-10), null);
  for (const patch of [{ step: .1 }, { activeChange: true }, { accepted: false }, { residual: 1e-12 },
    { residualDecrease: { beforeSquaredNorm: 1, afterSquaredNorm: .1 } }]) {
    const trial = structuredClone(rows); Object.assign(trial[1], patch);
    assert.equal(coupledStepStagnation(trial, 1e-10), null);
  }
});
