// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { selectCoupledEulerStartup } from '../src/euler/streamtube-coupled-startup.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

function fixture() {
  const flow = solveStreamtubeIses(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), {
    maxIterations: 0, iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', retainBestCheckpoint: true });
  assert.equal(flow.initialRedistribution.accepted, true);
  const expected = structuredClone(flow);
  // Only manufacture the terminal comparison envelope. The retained state
  // and its residual/geometry are real and must survive exact kernel replay.
  flow.history.push({ iteration: 20, residual: 2 * flow.diagnostics.residual });
  flow.diagnostics.residual *= 2;
  flow.finalQuality.minCornerSine *= .5;
  return { precursor: { status: 'unconverged', flow }, expected };
}

test('earlier inviscid seed replays exactly without iterations and preserves terminal evidence', () => {
  const { precursor, expected } = fixture(), before = structuredClone(precursor);
  const selected = selectCoupledEulerStartup(precursor);
  assert.ok(selected);
  assert.deepEqual(selected.flow.x, expected.x);
  assert.deepEqual(selected.flow.nodes, expected.nodes);
  assert.deepEqual(selected.flow.residual, expected.residual);
  assert.equal(selected.flow.linearDiagnostics.solves, 0);
  assert.equal(selected.flow.initialRedistribution.resumed, true);
  assert.equal(selected.selection.selectedIteration, 0);
  assert.equal(selected.selection.terminalIteration, 20);
  assert.equal(selected.selection.stateConverged, false);
  assert.deepEqual(precursor, before);
});

test('converged, rejected, unstarted or nondominating terminal states keep their original endpoint', () => {
  for (const change of [
    p => { p.status = 'research-converged'; },
    p => { p.flow.converged = true; },
    p => { p.flow.reason = 'admissibility failure'; },
    p => { p.flow.lastRejectedStep = { stage: 'Newton' }; },
    p => { p.flow.initialRedistribution.accepted = false; },
    p => { p.flow.finalQuality.valid = false; },
    p => { p.flow.finalQuality.minCornerSine = NaN; },
    p => { p.flow.history.pop(); },
    p => { p.flow.diagnostics.residual = p.flow.bestCheckpoint.residual; },
    p => { p.flow.finalQuality.minCornerSine = 1; },
  ]) {
    const { precursor } = fixture(); change(precursor);
    const before = structuredClone(precursor);
    assert.equal(selectCoupledEulerStartup(precursor), null);
    assert.deepEqual(precursor, before);
  }
});

test('corrupted earlier evidence fails its actual gas/residual replay', () => {
  for (const corrupt of [
    p => { p.flow.bestCheckpoint.checkpoint.residual[0] += .01; },
    p => { p.flow.bestCheckpoint.checkpoint.initialEuler.x[0] = 1000; },
    p => { p.flow.bestCheckpoint.residual *= .5; },
  ]) {
    const { precursor } = fixture(); corrupt(precursor);
    assert.throws(() => selectCoupledEulerStartup(precursor));
  }
});

test('a different complete input or update method cannot substitute for the retained precursor', () => {
  for (const change of [
    p => { p.flow.bestCheckpoint.checkpoint.input.mach += .01; },
    p => { p.flow.bestCheckpoint.checkpoint.input.bodies[0].points[0].y += .001; },
    p => { p.flow.bestCheckpoint.checkpoint.continuation.stagnationLimiter = 'prose'; },
  ]) {
    const { precursor } = fixture(); change(precursor);
    assert.throws(() => selectCoupledEulerStartup(precursor), /different geometry, equations or update controls/);
  }
});

test('adaptive precursor selection replays matching laws and does not rank different MCRIT residuals', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'momentum', upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const flow = solveStreamtubeIses(input, { maxIterations: 0, iterationGeometry: 'ises-sampled',
    stepAcceptance: 'armijo', adaptiveMcrit: true, retainBestCheckpoint: true });
  flow.history.push({ iteration: 20, residual: 2 * flow.diagnostics.residual });
  flow.diagnostics.residual *= 2;
  flow.finalQuality.minCornerSine *= .5;
  const precursor = { status: 'unconverged', flow };
  assert.ok(selectCoupledEulerStartup(precursor));
  flow.solverInput.upwind = { ...flow.solverInput.upwind, mcrit: .9 };
  const before = structuredClone(precursor);
  assert.equal(selectCoupledEulerStartup(precursor), null);
  assert.deepEqual(precursor, before);
});

test('a partial precursor cannot rank a first-order root against second-order residuals', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'momentum', upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const flow = solveStreamtubeIses(input, { maxIterations: 0, iterationGeometry: 'ises-sampled',
    stepAcceptance: 'armijo', firstOrderStartup: true, retainBestCheckpoint: true });
  flow.history.push({ iteration: 20, residual: 2 * flow.diagnostics.residual });
  flow.diagnostics.residual *= 2;
  flow.finalQuality.minCornerSine *= .5;
  const precursor = { status: 'unconverged', flow };
  assert.ok(selectCoupledEulerStartup(precursor));
  flow.solverInput.upwind = { ...flow.solverInput.upwind, mucon: 1 };
  const before = structuredClone(precursor);
  assert.equal(selectCoupledEulerStartup(precursor), null);
  assert.deepEqual(precursor, before);
  flow.bestCheckpoint.checkpoint.continuation.targetMucon = 2;
  assert.throws(() => selectCoupledEulerStartup(precursor), /different geometry, equations or update controls/);
});
