// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { solveStreamtubeIsesAutomatic } from '../src/euler/tests/streamtube-ises-automatic.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

let source;
function accepted() {
  if (!source) {
    const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2, mach: .2 }),
      flowModel: 'compressible', streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
      upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
    source = solveStreamtubeIses(input, { retainCheckpoint: true, maxIterations: 20,
      stepAcceptance: 'admissible', tolerance: 1e-10 });
    assert.equal(source.converged, true, source.reason);
  }
  return structuredClone(source.checkpoint);
}

test('same-Mach accepted checkpoint returns exactly with no target stage or LU', () => {
  const checkpoint = accepted(), before = structuredClone(checkpoint);
  const result = solveStreamtubeIsesAutomatic(.2, { initialCheckpoint: checkpoint, stageMaxIterations: 0,
    onStage: () => assert.fail('No target stage for identical Mach.') });
  assert.equal(result.converged, true);
  assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(result.continuation.reachedTarget, true);
  assert.deepEqual(result.continuation.attempts, []);
  assert.deepEqual(result.checkpoint, before); assert.deepEqual(checkpoint, before);
  assert.deepEqual(result.initialRedistribution.passages, []);
});

test('real hybrid .2 to .25 to .3 continuation preserves complete transfer and maintenance history', () => {
  let checkpoint = accepted();
  for (const target of [.25, .3]) {
    const before = structuredClone(checkpoint), meshes = [], checkpoints = [], iterations = [];
    const system = createStreamtubeBodySystem(checkpoint.input), ne = system.layout.densityCount;
    const result = solveStreamtubeIsesAutomatic(target, { initialCheckpoint: checkpoint,
      onIteration: h => iterations.push(h), onMesh: m => meshes.push(m),
      onCheckpoint: (c, details) => checkpoints.push({ checkpoint: c, details }) });
    assert.equal(result.converged, true, result.reason);
    assert.equal(result.stateConverged, true); assert.equal(result.status, 'research-euler-ises-converged');
    assert.equal(result.continuation.reachedTarget, true);
    assert.equal(result.continuation.currentMach, target);
    assert.equal(result.checkpoint.input.mach, target);
    assert.equal(result.finalQuality.valid, true);
    assert.equal(result.continuation.attempts.length, 1);
    assert.ok(result.continuation.attempts.every(a => !a.initialSMOVERepeated
      && a.transfer.physicalMassPreserved && a.transfer.physicalDensityPreserved && a.transfer.maximumGeometryChange === 0));
    assert.ok(iterations.length > 1 && meshes.length === iterations.length);
    assert.equal(checkpoints.filter(c => c.details.kind === 'iterate').length, iterations.length);
    assert.equal(checkpoints.filter(c => c.details.kind === 'accepted').length, 1);
    assert.equal(checkpoints.at(-1).details.stateConverged, true);
    assert.equal(checkpoints.at(-1).details.quality.valid, true);
    assert.deepEqual(checkpoints.at(-1).checkpoint, result.checkpoint);
    const transferred = checkpoints[0].checkpoint;
    assert.deepEqual(transferred.initialEuler.x.slice(0, ne), before.initialEuler.x.slice(0, ne));
    assert.deepEqual(transferred.initialEuler.nodes, before.initialEuler.nodes);
    assert.deepEqual(transferred.continuation, before.continuation);
    assert.notDeepEqual(transferred.residual, before.residual);
    assert.deepEqual(checkpoint, before);
    assert.ok(result.history.slice(1).every(h => h.maintenance.inlet));
    assert.match(result.formulation, /hybrid.*shared MSES speed upwinding/);
    assert.doesNotMatch(result.formulation, /no evolving BL or shock dissipation/);
    assert.doesNotThrow(() => structuredClone(meshes));
    checkpoint = result.checkpoint;
  }
});

test('bounded failed stages halve the increment and retain the last accepted source, not target success', () => {
  const checkpoint = accepted(), before = structuredClone(checkpoint);
  const result = solveStreamtubeIsesAutomatic(.3, { initialCheckpoint: checkpoint,
    stageMaxIterations: 0, maxSubdivisions: 2, maxStages: 8 });
  assert.equal(result.converged, false); assert.equal(result.continuation.reachedTarget, false);
  assert.equal(result.stateConverged, true); assert.equal(result.status, 'research-euler-ises-target-not-reached');
  assert.equal(result.continuation.currentMach, .2); assert.equal(result.continuation.targetMach, .3);
  const attempted = result.continuation.attempts.map(a => a.mach);
  assert.equal(attempted.length, 3);
  for (let i = 0; i < attempted.length; i++) assert.ok(Math.abs(attempted[i] - (.2 + .1 / 2 ** i)) < 1e-15);
  assert.ok(result.continuation.attempts.every(a => !a.converged && a.fromMach === .2 && a.linearSolves === 0));
  assert.deepEqual(result.checkpoint, before); assert.deepEqual(checkpoint, before);
  const capped = solveStreamtubeIsesAutomatic(.3, { initialCheckpoint: checkpoint, stageMaxIterations: 0, maxStages: 1 });
  assert.equal(capped.continuation.attempts.length, 1); assert.match(capped.reason, /stage limit/);
});

test('all observers propagate cancellation without silently starting another stage', () => {
  for (const name of ['onStage', 'onIteration', 'onMesh', 'onCheckpoint']) {
    let stages = 0, calls = 0;
    const error = new Error(`cancel ${name}`);
    const callbacks = { onStage: () => { stages++; } };
    callbacks[name] = (...args) => {
      if (name === 'onStage') stages++;
      calls++;
      const iteration = name === 'onCheckpoint' ? args[1].history.at(-1).iteration
        : name === 'onIteration' ? args[0].iteration : name === 'onMesh' ? args[0].iteration.iteration : 1;
      if (iteration >= 1) throw error;
    };
    assert.throws(() => solveStreamtubeIsesAutomatic(.25, { initialCheckpoint: accepted(), ...callbacks }), e => e === error);
    assert.equal(stages, 1); assert.ok(calls >= 1);
  }
  assert.throws(() => solveStreamtubeIsesAutomatic(.25, { initialCheckpoint: accepted(), onStage: () => { throw null; } }), e => e === null);
  const acceptedError = new Error('cancel after converged stage');
  assert.throws(() => solveStreamtubeIsesAutomatic(.25, { initialCheckpoint: accepted(),
    onCheckpoint: (_, details) => { if (details.kind === 'accepted') throw acceptedError; } }), e => e === acceptedError);
});

test('observer mutations cannot alter input, accepted data or subsequent stages', () => {
  const initialCheckpoint = accepted(), expected = solveStreamtubeIsesAutomatic(.25, { initialCheckpoint });
  const actual = solveStreamtubeIsesAutomatic(.25, { initialCheckpoint,
    onStage: stage => { initialCheckpoint.input.mach = 99; stage.mach = 99; },
    onIteration: h => { h.residual = 99; if (h.maintenance) h.maintenance.inlet.maxDisplacement = 99; },
    onMesh: m => { m.mesh.quality.valid = false; m.iteration.residual = 99; },
    onCheckpoint: (checkpoint, details) => {
      checkpoint.input.mach = 99; checkpoint.initialEuler.x.fill(99); checkpoint.initialEuler.nodes[0][0][0].x = 99;
      checkpoint.continuation.fractions[0][0] = 99; if (details.history) details.history.length = 0;
    } });
  assert.deepEqual(actual.checkpoint, expected.checkpoint);
  assert.deepEqual(actual.continuation, expected.continuation);
  assert.deepEqual(actual.history, expected.history);
});

test('invalid controls and corrupt or unconverged sources fail before recovery', () => {
  const checkpoint = accepted(); let stages = 0;
  const run = (initialCheckpoint = checkpoint, extra = {}, target = .25) => solveStreamtubeIsesAutomatic(target,
    { initialCheckpoint, onStage: () => { stages++; }, ...extra });
  for (const target of [0, 1, NaN]) assert.throws(() => run(checkpoint, {}, target), /controls/);
  for (const controls of [{ stageMaxIterations: -1 }, { maxStages: 0 }, { maxStages: 257 }, { maxSubdivisions: 21 }, { tolerance: 0 }])
    assert.throws(() => run(checkpoint, controls), /controls/);
  assert.throws(() => run({}), /complete accepted checkpoint/);
  const broken = mutate => { const c = structuredClone(checkpoint); mutate(c); return c; };
  assert.throws(() => run(broken(c => c.residual[0] += 1)), /residual does not replay exactly/);
  assert.throws(() => run(broken(c => c.initialEuler.nodes[0][1][1].y += 1e-4)), /does not replay exactly/);
  assert.throws(() => run(broken(c => c.continuation.fractions[0][1] = 0)), /inlet fractions/);
  assert.throws(() => run(broken(c => c.continuation.stepAcceptance = 'invalid')), /controls/);
  const unfinished = solveStreamtubeIses(checkpoint.input, { initialEuler: checkpoint.initialEuler,
    maxIterations: 0, retainCheckpoint: true, stepAcceptance: 'admissible' });
  assert.equal(unfinished.converged, false);
  assert.throws(() => run(unfinished.checkpoint), /source must be converged/);
  assert.equal(stages, 0);
});
