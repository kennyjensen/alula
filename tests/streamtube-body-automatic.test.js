// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { solveStreamtubeBodyAutomatic } from '../src/euler/tests/streamtube-body-automatic.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const controls = { tolerance: 1e-10, directMaxIterations: 15, stageMaxIterations: 15 };
const fixture = () => intrinsicBodyFixture({ bodySegments: 4, tubes: 2, elements: 1, alpha: 0, mach: .2 });
const maximum = values => Math.max(...Array.from(values, Math.abs));
function targetRoot(result, mach) {
  assert.ok(result.system.layout.n < 150, 'Keep automatic-wrapper Newton tests on a tiny controlled grid.');
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.continuation.reachedTarget, true);
  assert.equal(result.system.conditions.mach, mach);
  assert.equal(result.continuation.targetMach, mach);
  assert.equal(result.continuation.currentMach, mach);
  assert.equal(result.continuation.analyticShockInformationUsed, false);
  assert.equal(result.continuation.densityReinitializedOnWarmRestart, false);
  assert.equal(result.physicalAcceptance, false);
  assert.equal(result.fullSolver, false);
  assert.ok(maximum(result.residual) < controls.tolerance);
  assert.ok(result.x.every(Number.isFinite));
  assert.equal(streamtubeMeshSnapshot({ system: result.system, nodes: result.nodes }).quality.valid, true);
}

test('automatic body startup converges the requested low-Mach conservative body directly', () => {
  const input = fixture(), before = structuredClone(input);
  const result = solveStreamtubeBodyAutomatic(input, controls);
  targetRoot(result, input.mach);
  assert.equal(result.continuation.attempts.length, 1);
  assert.equal(result.continuation.attempts[0].label, 'target-direct');
  assert.equal(result.continuation.attempts[0].startup.initialGuessOnly, true);
  assert.equal(result.system.conditions.streamwiseMode, 'momentum');
  assert.equal(result.system.conditions.upwind.boundary.kind, 'unfiltered-first-two');
  assert.deepEqual(input, before);
});

test('automatic body budgets cannot label a failed or intermediate-Mach root as the requested solution', () => {
  const input = fixture();
  const zero = solveStreamtubeBodyAutomatic(input,
    { ...controls, directMaxIterations: 0, stageMaxIterations: 0, maxStages: 2 });
  assert.equal(zero.converged, false);
  assert.equal(zero.continuation.reachedTarget, false);
  assert.equal(zero.continuation.attempts.length, 2);
  assert.ok(zero.continuation.attempts.every(a => !a.converged && a.iterations === 0));
  assert.match(zero.reason, /bounded attempts/);

  const capped = solveStreamtubeBodyAutomatic(input, { ...controls, directMaxIterations: 0, maxStages: 2 });
  assert.equal(capped.converged, false);
  assert.equal(capped.continuation.reachedTarget, false);
  assert.equal(capped.continuation.attempts.length, 2);
  assert.equal(capped.continuation.attempts[1].label, 'subcritical-start');
  assert.equal(capped.continuation.attempts[1].converged, true);
  assert.ok(maximum(capped.residual) < controls.tolerance,
    'A converged intermediate residual must still report failure to reach the requested Mach.');
  assert.equal(capped.system.conditions.mach, input.mach / 2);
  assert.equal(capped.continuation.currentMach, input.mach / 2);
  assert.equal(capped.continuation.targetMach, input.mach);
  assert.match(capped.reason, /stage limit/);
});

test('automatic Mach fallback preserves the accepted physical flow and isolates scalar callbacks and caller edits', () => {
  const input = fixture(), before = structuredClone(input), stages = [], meshIterations = [];
  const result = solveStreamtubeBodyAutomatic(input, { ...controls, directMaxIterations: 0,
    initialFractionStep: .5, maxStages: 6,
    onStage: event => {
      stages.push({ ...event });
      event.label = 'observer mutation'; event.fraction = -1; event.mach = NaN; event.targetMach = NaN;
      // These arrays and scalars belong to the caller, not the private case
      // snapshot that must define every later continuation stage.
      input.mach = NaN; input.alpha = 45; input.bodies[0].points[0].y = 100;
      input.outerLower[0].x = NaN; input.weights[0][0] = -1;
    },
    onIteration: event => {
      event.residual = NaN; event.step = NaN; event.mach = NaN;
      event.targetMach = NaN; event.continuationFraction = -1;
    },
    onMesh: event => {
      meshIterations.push({ stage: event.stage, mach: event.mach, targetMach: event.targetMach,
        fraction: event.continuationFraction, iteration: { ...event.iteration } });
      for (const key of ['system', 'flow', 'nodes']) assert.equal(Object.hasOwn(event, key), false);
      assert.equal(event.mesh.quality.valid, true);
      event.mach = NaN; event.targetMach = NaN; event.continuationFraction = -1;
      event.iteration.residual = NaN; event.iteration.iteration = -1;
      // Only detached display data is exposed; corrupting every coordinate
      // and connectivity entry must not alter this or the next flow stage.
      event.mesh.vertices.forEach(p => { p.x = NaN; p.y = NaN; });
      event.mesh.cells.forEach(cell => cell.fill(-1));
      event.mesh.quality.valid = false;
    },
  });
  targetRoot(result, before.mach);
  const attempts = result.continuation.attempts;
  assert.equal(attempts.length, 4);
  assert.equal(attempts[0].converged, false);
  assert.ok(attempts.slice(1).every(a => a.converged));
  assert.deepEqual(attempts.map(a => a.mach), [.2, .1, .15000000000000002, .2]);
  assert.deepEqual(attempts.map(({ label, mach, fraction }) => ({ label, mach, fraction, targetMach: before.mach })), stages);
  assert.ok(attempts.every(a => a.history.every(h => Number.isFinite(h.residual))));
  assert.ok(meshIterations.length > 0);
  assert.ok(meshIterations.every(e => Number.isFinite(e.mach) && e.targetMach === before.mach
    && e.iteration.iteration > 0 && Number.isFinite(e.iteration.residual)));
  assert.equal(result.system.conditions.alpha, before.alpha);
  assert.deepEqual(result.system.layout.bodies[0].points, before.bodies[0].points);
  for (const attempt of attempts.filter(a => a.label === 'mach-continuation')) {
    assert.equal(attempt.startup.physicalMassPreserved, true);
    assert.equal(attempt.startup.physicalDensityPreserved, true);
    assert.match(attempt.startup.densityInitialization, /no isentropic inversion/);
    assert.equal(attempt.startup.targetMach, attempt.mach);
    assert.ok(attempt.startup.maximumGeometryChange < 1e-12);
    assert.ok(attempt.startup.maximumSpeedChange < 1e-12);
  }
});

test('body continuation propagates observer exceptions without starting a recovery stage', () => {
  for (const callback of ['onStage', 'onIteration', 'onMesh']) {
    const error = new Error(`cancel from ${callback}`);
    let stages = 0, calls = 0;
    const options = { ...controls, onStage: () => { stages++; } };
    options[callback] = () => { calls++; throw error; };
    assert.throws(() => solveStreamtubeBodyAutomatic(fixture(), options), caught => caught === error);
    assert.equal(calls, 1);
    assert.equal(stages, callback === 'onStage' ? 0 : 1);
  }
});

test('automatic body startup rejects invalid budgets and unsupported physical modes before solving', () => {
  for (const invalid of [
    { directMaxIterations: -1 }, { directMaxIterations: .5 }, { stageMaxIterations: -1 },
    { maxSubdivisions: -1 }, { maxSubdivisions: 21 }, { maxStages: 1 }, { maxStages: 2.5 },
    { initialFractionStep: 0 }, { initialFractionStep: 1.1 }, { tolerance: 0 }, { tolerance: NaN },
  ]) assert.throws(() => solveStreamtubeBodyAutomatic(fixture(), { ...controls, ...invalid }), /continuation controls/);
  for (const invalid of [{ flowModel: 'incompressible' }, { streamwiseMode: 'isentropic' }])
    assert.throws(() => solveStreamtubeBodyAutomatic({ ...fixture(), ...invalid }, controls), /compressible conservative momentum/);
  assert.throws(() => solveStreamtubeBodyAutomatic(fixture(), { ...controls, initialEuler: {} }), /geometry-only seed/);
});
