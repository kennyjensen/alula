// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { recoverStreamtubeStagnationMotion } from '../src/euler/streamtube-stagnation-startup.js';

test('stagnation chart recovery preserves the physical state, closes the equations and strictly replays', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    stagnationMotion: 'walls-only' };
  const controls = { stepAcceptance: 'armijo', iterationGeometry: 'ises-sampled', tolerance: 1e-10,
    retainCheckpoint: true };
  const partial = solveStreamtubeIses(input, { ...controls, maxIterations: 1 });
  assert.equal(partial.converged, false);
  // Inject only the trigger; geometry, residuals, flow and subsequent solves
  // are real. This tests the recovery independently of a particular plateau.
  partial.lastRejectedStep = { code: 'EULER_GRID_STAGNATION', diagnostics: { cell: { group: 1, i: 1, tube: 0 } } };
  const before = structuredClone(partial), frames = [];
  const result = recoverStreamtubeStagnationMotion(partial, { maxIterations: 20, tolerance: 1e-10,
    onIteration: h => frames.push(h) });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.solverInput.stagnationMotion, 'interpolated');
  assert(result.stagnationRecovery.maximumCoordinateChange < 1e-14);
  assert(result.history.at(-1).iteration <= 20);
  assert.equal(frames[0].iteration, 1);
  assert.deepEqual(result.history.map(h => h.iteration), Array.from({ length: result.history.length }, (_, i) => i));
  assert.deepEqual(result.solverInput.bodies, partial.solverInput.bodies);
  assert.equal(result.solverInput.mach, input.mach);
  assert.deepEqual(partial, before);
  const replay = solveStreamtubeIses(undefined, { ...result.checkpoint.continuation,
    resume: result.checkpoint, maxIterations: 0, tolerance: 1e-10 });
  assert.equal(replay.converged, true);
  assert.equal(replay.diagnostics.residual, result.diagnostics.residual);
  assert.equal(recoverStreamtubeStagnationMotion(partial, { maxIterations: 1 }), partial);
  assert.equal(recoverStreamtubeStagnationMotion({ ...partial, lastRejectedStep: null }, { maxIterations: 20 }).checkpoint,
    partial.checkpoint);
});
