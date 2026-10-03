// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

test('a proven DSLIM descent failure is recovered before slow gradient fallback masks it', () => {
  const { checkpoint: cp } = JSON.parse(gunzipSync(fs.readFileSync(
    new URL('./fixtures/nlr-projected-newton-stall.json.gz', import.meta.url))));
  const before = structuredClone(cp), frames = [], checkpoints = [];
  const result = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp,
    maxIterations: 20, tolerance: 1e-10, onIteration: h => frames.push(h),
    onCheckpoint: (checkpoint, details) => checkpoints.push({ checkpoint, details }) });
  const recovery = frames.find(h => h.hkProjectionRecovery);
  assert.ok(recovery);
  assert.equal(recovery.iteration, 1);
  assert.equal(recovery.accepted, false);
  assert.equal(recovery.step, 0);
  const proof = recovery.hkProjectionRecovery;
  assert.ok(proof.rawMerit < proof.beforeMerit);
  assert.ok(proof.projectedMerit >= proof.beforeMerit);
  assert.equal(proof.equationsChanged, false);
  const switched = checkpoints.find(p => p.checkpoint.restart.options.hkFloorLinearization === 'native').checkpoint;
  assert.deepEqual(switched.families, cp.families);
  assert.deepEqual(Array.from(switched.restart.initialEuler.x), cp.restart.initialEuler.x);
  assert.deepEqual(Array.from(switched.restart.initialBL), cp.restart.initialBL);
  assert.deepEqual(switched.restart.input, cp.restart.input);
  assert.equal(result.converged, true, result.reason);
  assert.ok(Object.values(result.families).every(v => v <= 1e-10));
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(result.history.at(-1).iteration <= 15);
  assert.deepEqual(cp, before);
});
