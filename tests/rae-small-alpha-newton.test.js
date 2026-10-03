// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { initializeCoupledStreamtubeFromFlow } from '../src/euler/tests/streamtube-coupled-flow-restart.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

test('RAE Mach .74, alpha .84 to .85 converges directly with unchanged topology and source', () => {
  const {checkpoint}=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-small-alpha-audit/seed.json'));
  const before=structuredClone(checkpoint);
  assert.equal(checkpoint.restart.input.mach,.74);
  assert.equal(checkpoint.restart.input.alpha,.84);
  assert.ok(Object.values(checkpoint.families).every(v=>v<=1e-10));
  const seed=initializeCoupledStreamtubeFromFlow(.74,checkpoint,{targetAlpha:.85});
  const r=solveCoupledStreamtubeIses(undefined,{resume:seed.checkpoint,...seed.checkpoint.continuation,
    tolerance:1e-10,maxIterations:6,iterationRecovery:false});
  assert.equal(r.converged,true,r.reason);
  assert.equal(r.checkpoint.restart.input.alpha,.85);
  assert.equal(r.checkpoint.restart.input.mach,.74);
  assert.ok(r.history.length<=5, 'nearby root should recover rapid terminal Newton convergence');
  assert.ok(r.history.slice(1).every(h=>h.step===1));
  assert.ok(Object.values(r.families).every(v=>v<=1e-10));
  assert.equal(r.checkpoint.restart.initialEuler.x.length,checkpoint.restart.initialEuler.x.length);
  assert.deepEqual(checkpoint,before);
});

test('RAE 32/11 to 64 profile preparation repairs one-ulp geometry replay without changing source', async () => {
  const {refineCoupledStreamtubeBody}=await import('../src/euler/streamtube-coupled-refinement.js');
  const {prepareRefinedCoupledGridProfile}=await import('../src/euler/streamtube-coupled-grid-profile.js');
  const {checkpoint}=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-small-alpha-audit/seed11.json'));
  const before=structuredClone(checkpoint),source=initializeCoupledStreamtubeFromFlow(.74,checkpoint);
  const refined=refineCoupledStreamtubeBody(checkpoint.restart.input,source.system,
    {streamwiseFactor:2,normalFactor:1,streamwiseInterpolation:'surface-pchip'});
  const saved=structuredClone({input:refined.input,options:refined.options,initialEuler:refined.initialEuler,initialBL:refined.initialBL});
  const r=prepareRefinedCoupledGridProfile(refined),d=r.transfer.profilePreparation;
  assert.equal(d.canonicalReplayExact,true);
  assert.equal(d.geometryReplayRecovery.accepted,true);
  assert.equal(d.geometryReplayRecovery.originalGatesChanged,false);
  assert.equal(d.geometryReplayRecovery.sourceFailure.diagnostics.geometryReplayDifference.changedCoordinates,1);
  assert.ok(d.geometryReplayRecovery.sourceFailure.diagnostics.geometryReplayDifference.maximumAbsoluteDifference<1e-17);
  assert.deepEqual(checkpoint,before);
  assert.deepEqual({input:refined.input,options:refined.options,initialEuler:refined.initialEuler,initialBL:refined.initialBL},saved);
});
