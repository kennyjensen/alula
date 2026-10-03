// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { coupledMachRecoveryCase, recoverCoupledMach } from '../src/euler/streamtube-coupled-mach-recovery.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
const input = { mach: .76, alpha: 2.68, gridIntervals: 16, gridTubes: 7, reynolds: 1e6,
  ncrit: 9, transitionMode: 'automatic', materialTrips: [[1, 1]] };
const failed = () => ({ converged: false, checkpoint: { restart: { options: { edgeMatching: 'section-velocity' } } },
  mesh: { quality: { valid: true } }, history: [{ iteration: 0 }] });

test('Mach recovery changes only the startup Mach and cannot recurse or replace a root', () => {
  const source = coupledMachRecoveryCase(input, failed(), { maxIterations: 40 });
  assert.deepEqual(source, { ...input, mach: .74 });
  source.materialTrips[0][0] = .2;
  assert.equal(input.materialTrips[0][0], 1);
  for (const patch of [{ mach: .2 }, { mach: .9 }, { transitionMode: 'fixed-trip' }, { materialTrips: [[.2, 1]] }])
    assert.equal(coupledMachRecoveryCase({ ...input, ...patch }, failed(), { maxIterations: 40 }), null);
  for (const patch of [{ converged: true }, { machRecovery: {} }, { mesh: { quality: { valid: false } } }])
    assert.equal(coupledMachRecoveryCase(input, { ...failed(), ...patch }, { maxIterations: 40 }), null);
  assert.equal(coupledMachRecoveryCase(input, failed(), { enabled: false, maxIterations: 40 }), null);
  assert.equal(coupledMachRecoveryCase(input, failed(), { maxIterations: 0 }), null);
});

function source() {
  const f = JSON.parse(gunzipSync(fs.readFileSync(new URL('./fixtures/rae16x7-resolution-source.json.gz', import.meta.url))));
  const cp = f.checkpoint;
  return solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 0, tolerance: 1e-10 });
}

test('Mach recovery requires a complete root at the requested Mach after both transfers', () => {
  const parent = source(), before = structuredClone(parent.checkpoint), events = [];
  let cold = 0;
  const recovered = recoverCoupledMach(failed(), input, { sourceInput: { ...input, mach: .74 },
    maxIterations: 40, tolerance: 1e-10, solve: (c, controls) => {
      cold++; assert.equal(c.mach, .74); assert.equal(controls.machRecovery, false);
      assert.equal(controls.resolutionRecovery, false); return parent;
    }, onStage: e => events.push(e) });
  assert.equal(cold, 1);
  assert.equal(recovered.diagnostics.accepted, true, recovered.diagnostics.reason);
  assert.deepEqual(recovered.diagnostics.attempts.map(a => a.mach), [.75, .76]);
  const r = recovered.result;
  assert.equal(r.conditions.mach, .76); assert.equal(r.solverInput.alpha, 2.68);
  assert.equal(r.conditions.ncrit, 9); assert.equal(r.mesh.quality.valid, true);
  for (const value of Object.values(r.families)) assert.ok(value <= 1e-10);
  assert.deepEqual(events.map(e => [e.actualMach, e.targetMach]), [[.75, .76], [.76, .76]]);
  assert.deepEqual(parent.checkpoint, before);
});

test('a failed intermediate solve is not success and observer failures propagate', () => {
  const parent = source(), original = failed(), before = structuredClone(original);
  const options = { sourceInput: { ...input, mach: .74 }, maxIterations: 0, tolerance: 1e-10, solve: () => parent };
  const r = recoverCoupledMach(original, input, options);
  assert.equal(r.result, undefined); assert.equal(r.diagnostics.accepted, false);
  assert.equal(r.diagnostics.attempts.length, 1);
  assert.deepEqual(original, before);
  const error = new Error('observer interrupted');
  assert.throws(() => recoverCoupledMach(original, input, { ...options, onStage: () => { throw error; } }), e => e === error);
});
