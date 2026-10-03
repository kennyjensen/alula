// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { coupledConstrainedDirection } from '../src/euler/streamtube-coupled-direction-recovery.js';
import { sparseMatrix, sparseAdd, sparseProduct } from '../src/numerics/sparse.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

function model(units = [1, 1]) {
  const matrix = sparseMatrix([[0, 1], [0, 1]]);
  for (const [i, j, value] of [[0, 0, 1], [0, 1, 2], [1, 1, 1]]) sparseAdd(matrix, i, j, value / units[j]);
  return { matrix, system: { stepConstraints: () => [{ gradient: new Map([[1, 1 / units[1]]]), lower: 0 }] },
    state: Float64Array.of(0, 0), residual: Float64Array.of(1, 1), newton: Float64Array.of(units[0], -units[1]) };
}

test('constrained fallback finds descent when the Newton ray points outside the domain', () => {
  for (const units of [[1, 1], [1e-6, 1e4]]) {
    const { system, state, matrix, residual, newton } = model(units), before = structuredClone(matrix);
    const p = coupledConstrainedDirection(system, state, matrix, residual);
    assert.ok(p);
    assert.ok(p.direction[0] < 0);
    assert.ok(p.direction[1] >= -1e-12 * units[1]);
    const jd = sparseProduct(matrix, p.direction);
    assert.ok(jd.reduce((sum, v, i) => sum + v * residual[i], 0) < 0);
    assert.ok(p.diagnostics.projection.converged);
    assert.ok(p.diagnostics.predictedReduction > 0);
    assert.equal(p.diagnostics.equationsChanged, false);
    assert.deepEqual(matrix, before);
    assert.deepEqual(state, Float64Array.of(0, 0));
  }
});

test('constrained fallback uses the same weighted residual model without mutating rows', () => {
  const { system, state, matrix, residual, newton } = model(), weights = Float64Array.of(4, .5);
  const before = structuredClone({ matrix, residual, weights });
  const p = coupledConstrainedDirection(system, state, matrix, residual, { weights });
  const jd = sparseProduct(matrix, p.direction);
  const reduction = residual.reduce((sum, r, i) => sum - weights[i] ** 2 * (r * jd[i] + .5 * jd[i] ** 2), 0);
  assert.ok(Math.abs(reduction - p.diagnostics.predictedReduction) < 1e-12);
  assert.deepEqual({ matrix, residual, weights }, before);
});

test('a stationary constrained model is not an acceptable fallback', () => {
  const { state, matrix, residual, newton } = model();
  const system = { stepConstraints: () => [0, 1].map(i => ({ gradient: new Map([[i, 1]]), lower: 0 })) };
  assert.equal(coupledConstrainedDirection(system, state, matrix, residual), null);
});

test('the stalled flap escapes its flattened-cell Newton ray without relaxing its equations', () => {
  const path = new URL('./fixtures/coupled-flattened-passage.json.gz', import.meta.url);
  const cp = JSON.parse(gunzipSync(fs.readFileSync(path))), before = structuredClone(cp);
  const r = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 3, tolerance: 1e-10 });
  const h = r.history.at(-1);
  assert.equal(h.directionRecovery.method, 'constrained-residual-gradient');
  assert.equal(h.step, 1);
  assert.ok(h.residualDecrease.afterSquaredNorm < .7 * h.residualDecrease.beforeSquaredNorm);
  assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.converged, false); // Escaping a stall is not a root.
  assert.deepEqual(r.checkpoint.restart.input, cp.restart.input);
  assert.deepEqual(cp, before);
});

test('an uncertified Newton solve can recover using a direction that needs no factorization', () => {
  const path = new URL('./fixtures/coupled-uncertified-newton.json.gz', import.meta.url);
  const cp = JSON.parse(gunzipSync(fs.readFileSync(path))), before = structuredClone(cp);
  const controls = { ...cp.continuation, resume: cp, maxIterations: 1, tolerance: 1e-10 };
  const stopped = solveCoupledStreamtubeIses(undefined, { ...controls, directionRecovery: false });
  assert.equal(stopped.lastRejectedStep.code, 'KLU_RESIDUAL_LIMIT');
  const r = solveCoupledStreamtubeIses(undefined, controls), h = r.history.at(-1);
  assert.equal(h.directionRecovery.cause, 'uncertified-newton-solve');
  assert.ok(h.directionRecovery.linearizedMeritSlope < 0);
  assert.ok(h.residualDecrease.afterSquaredNorm < h.residualDecrease.beforeSquaredNorm);
  assert.equal(r.linearDiagnostics.iterations[0].rejectedLinearSolve.code, 'KLU_RESIDUAL_LIMIT');
  assert.equal(r.linearDiagnostics.solves, 0); // No inaccurate LU is relabeled as certified.
  assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.converged, false);
  assert.deepEqual(cp, before);
});
