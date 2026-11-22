// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { recoverCoupledResolution } from '../src/euler/streamtube-coupled-resolution-recovery.js';
import { coupledChangeDecision } from '../src/euler/streamtube-coupled-convergence.js';

// Exercise the recovery acceptance contract without a numerical solve.
const input = { gridIntervals: 8 };
function source() {
  const families = { euler: 2e-7, boundaryLayer: 1e-7, edgeMatching: 3e-7 };
  const groups = { sample: { rms: 1e-6, maximum: 2e-6, count: 2 } };
  return { converged: true, checkpoint: {}, families,
    convergence: coupledChangeDecision({ accepted: groups, newton: groups, residual: groups, families, stable: true }),
    mesh: { quality: { valid: true }, cells: [0] }, history: [{ iteration: 0 }, { iteration: 1 }], initialization: {} };
}
function recover(parent) {
  return recoverCoupledResolution({ history: [{ iteration: 0 }], families: { euler: .1 } }, input,
    { sourceInput: { gridIntervals: 16 }, tolerance: 1e-10, maxIterations: 40, solve: () => parent });
}

test('resolution recovery accepts verified MSES changes without demanding an equation root', () => {
  const parent = source();
  const recovered = recover(parent);
  assert.equal(recovered.diagnostics.accepted, true);
  assert.equal(recovered.result.converged, true);
  assert.deepEqual(recovered.result.families, parent.families);
  assert.equal(recovered.diagnostics.finalNominalGridIntervals, 16);
});

test('recovery rejects damped, stale, incomplete, invalid or unconverged sources', () => {
  for (const mutate of [
    p => { delete p.convergence.newton; },
    p => { p.convergence.stable = false; },
    p => { p.convergence.newton = { sample: { rms: .01, maximum: .1, count: 2 } }; },
    p => { p.families.euler *= 2; },
    p => { p.families.euler = NaN; },
    p => { p.mesh.quality.valid = false; },
    p => { p.converged = false; },
    p => { delete p.checkpoint; },
  ]) {
    const parent = source(); mutate(parent);
    const recovered = recover(parent);
    assert.equal(recovered.diagnostics.accepted, false);
    assert.equal(recovered.result, undefined);
  }
});
