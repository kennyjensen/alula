// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const input = () => intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
const controls = { maxIterations: 2, tolerance: 1e-10, iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible' };

test('best Euler retention preserves the full terminal trajectory and exactly replays the minimum admitted state', () => {
  const plain = solveStreamtubeIses(input(), controls), captures = [];
  const retained = solveStreamtubeIses(input(), { ...controls, retainBestCheckpoint: true,
    onCheckpoint: (checkpoint, details) => captures.push({ checkpoint, iteration: details.history.length - 1,
      residual: details.history.at(-1).residual }) });
  assert(captures.length > 1, 'Exercise actual admitted Newton updates.');
  const expected = captures.reduce((best, value) => value.residual < best.residual ? value : best);
  assert.deepEqual(retained.bestCheckpoint, expected);
  const { bestCheckpoint, checkpoint, ...rest } = retained;
  assert.deepEqual(rest, plain, 'Retention must not change the steps, packed state, geometry, phase or reported terminal residual.');
  const saved = JSON.parse(JSON.stringify(bestCheckpoint.checkpoint));
  const replay = solveStreamtubeIses(undefined, { ...controls, maxIterations: 0, resume: saved });
  assert.equal(replay.initialRedistribution.resumed, true);
  assert.equal(replay.diagnostics.residual, bestCheckpoint.residual);
  assert.deepEqual(replay.nodes, saved.initialEuler.nodes);
  assert.deepEqual(Array.from(replay.residual), saved.residual);
  assert(replay.finalQuality.valid);
  // A detached retained seed cannot mutate the terminal result or input.
  const firstX = retained.nodes[0][0][0].x;
  retained.bestCheckpoint.checkpoint.initialEuler.nodes[0][0][0].x += 1;
  assert.equal(retained.nodes[0][0][0].x, firstX);
});

test('an admitted zero-update initialization is recorded honestly and controls reject unchecked geometry modes', () => {
  const result = solveStreamtubeIses(input(), { ...controls, maxIterations: 0, retainBestCheckpoint: true });
  assert.equal(result.bestCheckpoint.iteration, 0);
  assert.equal(result.history.length, 1);
  assert.equal(result.converged, false);
  assert.equal(result.reason, 'iteration limit');
  assert.equal(result.linearDiagnostics.solves, 0);
  for (const options of [{ retainBestCheckpoint: 'yes' }, { retainBestCheckpoint: true, stepAcceptance: 'listing' }])
    assert.throws(() => solveStreamtubeIses(input(), { ...controls, ...options }), /Invalid ISES research iteration controls/);
});
