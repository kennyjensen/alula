// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { coupledResolutionRecoveryCase, recoverCoupledResolution } from '../src/euler/streamtube-coupled-resolution-recovery.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

const input = { gridIntervals: 32, gridTubes: 7, gridInletIntervals: 16, gridOutletIntervals: 16,
  mach: .74, alpha: 2.68, reynolds: 1e6, ncrit: 9, transitionMode: 'automatic', materialTrips: [[1, 1]] };
const stalled = () => ({ converged: false, families: { euler: .1, boundaryLayer: .7, edgeMatching: .1 },
  checkpoint: { continuation: { stepAcceptance: 'event-armijo' } }, mesh: { quality: { valid: true } },
  history: Array.from({ length: 9 }, (_, iteration) => ({ iteration, step: .001 })) });
const load = tubes => JSON.parse(gunzipSync(fs.readFileSync(new URL(`./fixtures/rae16x${tubes}-resolution-source.json.gz`, import.meta.url))));
function source(tubes) {
  const f = load(tubes), cp = f.checkpoint;
  const r = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 0, tolerance: 1e-10 });
  assert.equal(r.converged, true);
  return { ...r, history: f.history, linearDiagnostics: f.linearDiagnostics,
    referenceChord: 1, solverLength: r.conditions.referenceChord, referenceReynolds: 1e6, kernelReynolds: r.conditions.reynolds,
    initialization: { boundaryLayer: {}, euler: {} } };
}

test('resolution recovery preserves every input except nominal surface resolution and remains bounded', () => {
  const before = structuredClone(input), r = stalled();
  const candidate = coupledResolutionRecoveryCase(input, r, { maxIterations: 40 });
  assert.deepEqual(candidate, { ...input, gridIntervals: 16 });
  candidate.materialTrips[0][0] = .3;
  assert.deepEqual(input, before);
  for (const change of [{ mach: .2 }, { gridIntervals: 16 }, { gridIntervals: 128 },
    { transitionMode: 'fixed-trip' }, { materialTrips: [[.2, 1]] }])
    assert.equal(coupledResolutionRecoveryCase({ ...input, ...change }, r, { maxIterations: 40 }), null);
  for (const change of [{ converged: true }, { resolutionRecovery: {} }, { mesh: { quality: { valid: false } } }])
    assert.equal(coupledResolutionRecoveryCase(input, { ...r, ...change }, { maxIterations: 40 }), null);
  assert.equal(coupledResolutionRecoveryCase(input, r, { maxIterations: 0 }), null);
  assert.equal(coupledResolutionRecoveryCase(input, r, { maxIterations: 40, enabled: false }), null);
});

for (const tubes of [7, 9]) test(`a converged ${16}x${tubes} source must solve again on its refined surface grid`, () => {
  const parent = source(tubes), original = structuredClone(parent.checkpoint), failed = stalled();
  const conditions = { ...input, gridTubes: tubes };
  let initializations = 0;
  const options = { sourceInput: { ...conditions, gridIntervals: 16 }, maxIterations: 40, tolerance: 1e-10,
    solve: (c, controls) => { initializations++; assert.equal(controls.resolutionRecovery, false);
      assert.deepEqual(c, { ...conditions, gridIntervals: 16 }); return parent; } };
  const recovered = recoverCoupledResolution(failed, conditions, options);
  assert.equal(initializations, 1);
  const r = recovered.result;
  assert.equal(r.converged, true, recovered.diagnostics.reason);
  assert.ok(r.mesh.cells.length > parent.mesh.cells.length);
  for (const residual of Object.values(r.families)) assert.ok(residual <= 1e-10);
  assert.equal(r.stepAcceptance, parent.stepAcceptance);
  assert.equal(r.mesh.quality.valid, true);
  assert.deepEqual(r.solverInput.upwind, parent.solverInput.upwind);
  assert.equal(r.solverInput.mach, .74);
  assert.equal(r.conditions.reynolds, parent.conditions.reynolds); assert.equal(r.conditions.ncrit, 9);
  assert.equal(r.history.length - 1, failed.history.length - 1 + parent.history.length - 1 + recovered.diagnostics.refinement.iterations);
  assert.deepEqual(parent.checkpoint, original);
  const b = parent.solverInput.bodies[0], counts = recovered.diagnostics.refinement.plan.streamwiseSubdivisions;
  assert.ok(counts.slice(0, b.leadingIndex).every(n => n === 1));
  assert.ok(counts.slice(b.trailingIndex).every(n => n === 1));
});

test('a coarse root is never accepted for a fine request whose refined solve has not converged', () => {
  const parent = source(7), failed = stalled(), before = structuredClone(failed);
  const r = recoverCoupledResolution(failed, input, { sourceInput: { ...input, gridIntervals: 16 },
    maxIterations: 0, tolerance: 1e-10, solve: () => parent });
  assert.equal(r.result, undefined);
  assert.equal(r.diagnostics.accepted, false);
  assert.equal(r.diagnostics.refinement.converged, false);
  assert.deepEqual(failed, before);
});

test('an adaptive increase on the coarsest request is explicit and observer failures propagate', () => {
  const parent = source(7), failed = stalled(), coarse = { ...input, gridIntervals: 8 };
  const args = { sourceInput: { ...input, gridIntervals: 16 }, maxIterations: 40, tolerance: 1e-10, solve: () => parent };
  const r = recoverCoupledResolution(failed, coarse, args);
  assert.equal(r.result.converged, true);
  assert.equal(r.diagnostics.requestedGridIntervals, 8);
  assert.equal(r.diagnostics.finalNominalGridIntervals, 16);
  assert.equal(r.diagnostics.requestedGridRetained, false);
  const error = new Error('observer failure');
  assert.throws(() => recoverCoupledResolution(failed, coarse, { ...args,
    solve: (_, hooks) => { hooks.onStage({ stage: 'coupled' }); return parent; },
    onStage: () => { throw error; } }), e => e === error);
  const bad = recoverCoupledResolution(failed, coarse, { ...args,
    solve: () => ({ ...parent, families: { ...parent.families, boundaryLayer: .1 } }) });
  assert.equal(bad.result, undefined); assert.equal(bad.diagnostics.accepted, false);
});
