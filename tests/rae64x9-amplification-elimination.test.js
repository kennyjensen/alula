// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeAssembly } from '../src/euler/streamtube-coupled-assembly.js';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';

test('RAE 64x9 cold handoff closes the original equations after eliminating laminar amplification', () => {
  const { checkpoint } = JSON.parse(gunzipSync(fs.readFileSync(new URL(
    './fixtures/rae64x9-amplification-start.json.gz', import.meta.url))));
  const before = structuredClone(checkpoint);
  // Select the new numerical update on the archived cold handoff. Old
  // checkpoints otherwise retain their original simultaneous-N policy.
  const seed = structuredClone(checkpoint);
  seed.continuation.eliminateAmplification = true;
  const result = solveCoupledStreamtubeIses(undefined, { ...seed.continuation,
    resume: seed, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.history.length <= 41, true);
  assert.equal(result.checkpoint.continuation.eliminateAmplification, true);
  assert.equal(result.conditions.mach, .74);
  assert.equal(result.conditions.ncrit, 4);
  assert.deepEqual(result.solverInput.bodies, checkpoint.restart.input.bodies);
  assert.equal(result.mesh.cells.length, 4032);
  assert.deepEqual(checkpoint, before);

  // Independently reconstruct every original governing row. A projected
  // amplification profile is not an alternative convergence criterion.
  const f = result.checkpoint.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: f.initialBL });
  const value = system.evaluate(system.initial);
  assert.equal(system.euler.layout.nx, 168);
  assert.deepEqual(system.euler.layout.tubes, [12, 12]);
  assert.deepEqual(value.families, result.families);
  assert.ok(value.residual.every(v => Number.isFinite(v) && Math.abs(v) <= 1e-10));
  assert.ok(system.bl.activeTargets(system.initial.subarray(0, system.ne),
    system.initial.subarray(system.ne)).every(t => t.from === t.to));
  const replay = solveCoupledStreamtubeIses(undefined, { ...result.checkpoint.continuation,
    resume: result.checkpoint, maxIterations: 0, tolerance: 1e-10 });
  assert.deepEqual(replay.residual, result.residual);
  assert.equal(replay.linearDiagnostics.solves, 0);
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { ...seed.continuation,
    resume: seed, eliminateAmplification: false, maxIterations: 0 }), /amplification-update controls/);
});

test('amplification elimination also converges from a cold coarse app input', () => {
  const { caseData } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-bl', changes: {
    alpha: 2.68, mach: .74, reynolds: 1e6, ncrit: 9, gridIntervals: 16, gridTubes: 7,
    gridInletIntervals: 16, gridOutletIntervals: 16,
  } });
  const result = solveCoupledStreamtubeAssembly(caseData, { direct: true, maxIterations: 40,
    eulerMaxIterations: 40, resolutionRecovery: false, machRecovery: false });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.mesh.cells.length, 1340);
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(result.history.length <= 41);
  assert.ok(Object.values(result.families).every(v => v <= 1e-10));
  assert.equal(result.checkpoint.continuation.eliminateAmplification, true);
});
