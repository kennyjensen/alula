// SPDX-License-Identifier: GPL-2.0-or-later
// Exercise the actual continuation controller with explicit fake numerical
// boundaries. These tests prove orchestration, not Euler/BL convergence.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { pathToFileURL } from 'node:url';
import path from 'node:path';

const source = fs.readFileSync('src/euler/tests/streamtube-coupled-automatic.js', 'utf8');
const clone = structuredClone;
const families = { euler: 1e-12, boundaryLayer: 2e-12, edgeMatching: 3e-12 };
const checkpoint = (mach = .2, marker = 'original') => ({
  version: 1, marker, families: { ...families },
  restart: { input: { mach, wakeGeometry: 'independent-banks', outerLower: Array(6), bodies: [{ trailingIndex: 3 }] },
    initialBL: [0, 1, 2, 3] },
  continuation: { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing' },
});
const result = (cp, valid = true) => ({
  conditions: { mach: cp.restart.input.mach }, checkpoint: clone(cp), converged: true,
  residual: [1e-12, 2e-12, 3e-12], families: { ...families }, status: 'converged',
  mesh: { quality: { valid, invalidCells: valid ? [] : [16] } },
  flow: { nodes: Array.from({ length: 2 }, () => Array.from({ length: 6 }, () => Array(3))),
    diagnostics: { maxMach: .3 } },
  linearDiagnostics: { solves: 0 }, initialRedistribution: { resumed: true }, history: [{}],
});
let id = 0;
async function controller(overrides = {}) {
  const calls = { transfers: [], stages: [], recoveries: [] };
  const fixture = {
    prepare(mach, cp) {
      calls.transfers.push(clone(cp));
      const target = clone(cp); target.restart.input.mach = mach;
      return { checkpoint: target, diagnostics: { sourceMach: cp.restart.input.mach } };
    },
    solve(_, options) {
      if (options.maxIterations === 0) return result(options.resume);
      calls.stages.push(clone(options.resume));
      const trial = overrides.solve?.(options, calls) ?? result(options.resume, false);
      if (!trial.mesh.quality.valid) {
        // Reproduce the newly eligible trigger: geometric failure before
        // numerical closure. This target must never become a recovery source.
        trial.residual[1] = 1.5; trial.families.boundaryLayer = 1.5; trial.converged = false;
      }
      return trial;
    },
    recover(cp, options) {
      calls.recoveries.push({ checkpoint: clone(cp), maxIterations: options.maxIterations });
      return overrides.recover?.(cp, options, calls)
        ?? { result: result(checkpoint(cp.restart.input.mach, 'recovered')), diagnostics: {} };
    },
    snapshot(value) { return clone(value.mesh); },
  };
  const key = `__wakeRecoveryControllerTest${++id}`;
  globalThis[key] = fixture;
  const proxy = 'data:text/javascript;base64,' + Buffer.from(`const h=globalThis[${JSON.stringify(key)}];
    export const initializeCoupledStreamtubeFromFlow=h.prepare,solveCoupledStreamtubeIses=h.solve,
      recoverCoupledWakeCorrespondence=h.recover,streamtubeMeshSnapshot=h.snapshot;
    export {coupledWakeRecoveryPlan} from ${JSON.stringify(pathToFileURL(path.resolve('src/euler/tests/streamtube-coupled-wake-recovery.js')).href)};`).toString('base64');
  try {
    const original = pathToFileURL(path.resolve('src/euler/tests/streamtube-coupled-automatic.js'));
    const controlled = new Set([
      './streamtube-coupled-flow-restart.js', '../streamtube-coupled-ises.js',
      '../streamtube-mesh-preview.js', './streamtube-coupled-wake-recovery.js',
    ]);
    const code = source.replace(/from '([^']+)'/g, (_, specifier) =>
      `from '${controlled.has(specifier) ? proxy : new URL(specifier, original).href}'`);
    const module = await import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
    return { run: module.solveCoupledStreamtubeAutomatic, calls };
  } finally { delete globalThis[key]; }
}

test('recovery has its own 12-update cap and actual-source progress, then retries the failed Mach from the recovered checkpoint', async () => {
  const initialCheckpoint = checkpoint(), before = clone(initialCheckpoint), events = [];
  const { run, calls } = await controller({
    solve: options => result(options.resume, options.resume.marker === 'recovered'),
    recover(cp, options) {
      const recovered = checkpoint(cp.restart.input.mach, 'recovered');
      options.onIteration({ iteration: 1 });
      options.onMesh({ mesh: { vertices: [{ x: 1 }] }, flow: { nodes: [1], bodyPressure() {} }, iteration: 1, coupledFamilies: families });
      options.onCheckpoint(recovered, {});
      return { result: result(recovered), diagnostics: {} };
    },
  });
  const terminal = run(.3, { initialCheckpoint, stageMaxIterations: 23, maxWakeRecoveries: 1, maxStages: 2,
    includeFlowState: true,
    onStage(info) { events.push(clone(info)); info.mach = 99; },
    onIteration(info) { assert.equal(info.mach, .2); assert.equal(info.fraction, 0); info.mach = 99; },
    onMesh(info, state) { assert.equal(info.mach, .2); assert.equal(state.flow.bodyPressure, undefined); state.flow.nodes[0] = 99; },
    onCheckpoint(cp, info) {
      if (info.kind === 'wake-initializer-iterate') { assert.equal(info.mach, .2); assert.equal(info.fraction, 0); }
      if (info.kind === 'accepted' && info.stage === 'coupled-wake-grid') assert.equal(info.reachedTarget, false);
      cp.marker = 'observer mutation';
    },
  });
  assert.equal(calls.recoveries[0].maxIterations, 12);
  assert.equal(calls.recoveries[0].checkpoint.restart.input.mach, .2);
  assert.deepEqual(calls.recoveries[0].checkpoint.families, families);
  assert.equal(calls.stages.length, 2);
  assert.deepEqual(calls.stages.map(cp => [cp.restart.input.mach, cp.marker]), [[.3, 'original'], [.3, 'recovered']]);
  assert.equal(terminal.converged, true); assert.equal(terminal.conditions.mach, .3);
  assert.equal(terminal.checkpoint.marker, 'recovered'); assert.deepEqual(initialCheckpoint, before);
  assert.equal(terminal.wakeRecovery.attempts[0].plan.residualConverged, false);
  const recovery = events.find(info => info.stage === 'coupled-wake-grid');
  assert.equal(recovery.mach, .2); assert.equal(recovery.targetMach, .3); assert.equal(recovery.fraction, 0);
});

test('recovery preserves the accepted source unless the replacement is converged, convex, same-Mach and within residual tolerance', async () => {
  for (const mutate of [r => { r.converged = false; }, r => { r.mesh.quality.valid = false; },
    r => { r.conditions.mach = .25; }, r => { r.residual[0] = 1e-3; }]) {
    const initialCheckpoint = checkpoint();
    const { run, calls } = await controller({ recover(cp) {
      const candidate = result(checkpoint(cp.restart.input.mach, 'unacceptable')); mutate(candidate);
      return { result: candidate, diagnostics: {} };
    } });
    const terminal = run(.3, { initialCheckpoint, maxWakeRecoveries: 1, maxStages: 3 });
    assert.equal(calls.recoveries.length, 1, 'A failed recovery consumes its bounded attempt.');
    assert.equal(terminal.converged, false); assert.equal(terminal.stateConverged, true);
    assert.equal(terminal.conditions.mach, .2); assert.equal(terminal.continuation.targetMach, .3);
    assert.deepEqual(terminal.checkpoint, initialCheckpoint);
    assert.ok(calls.stages.every(cp => cp.marker === 'original'));
  }
});

test('recovery callback cancellation propagates even when the numerical boundary swallows a thrown undefined', async () => {
  for (const hook of ['onStage', 'onIteration', 'onMesh', 'onCheckpoint']) {
    const initialCheckpoint = checkpoint(), before = clone(initialCheckpoint);
    const { run, calls } = await controller({ recover(cp, options) {
      try {
        if (hook === 'onIteration') options.onIteration({ iteration: 0 });
        if (hook === 'onMesh') options.onMesh({ mesh: {}, flow: {}, coupledFamilies: families });
        if (hook === 'onCheckpoint') options.onCheckpoint(cp, {});
      } catch { /* Simulate a driver converting an observer error to a rejected step. */ }
      return { result: result(checkpoint(.2, 'must not commit')), diagnostics: {} };
    } });
    let caught = false;
    try {
      run(.3, { initialCheckpoint, maxWakeRecoveries: 1, maxStages: 3,
        [hook]: (value, details) => {
          if (hook === 'onStage' && value.stage !== 'coupled-wake-grid') return;
          if (hook === 'onCheckpoint') assert.equal(details.kind, 'wake-initializer-iterate');
          throw undefined;
        } });
    } catch (error) { caught = true; assert.equal(error, undefined); }
    assert.equal(caught, true, hook); assert.equal(calls.stages.length, 1);
    assert.deepEqual(initialCheckpoint, before);
  }
});

test('recovery is opt-in and its bounded control rejects invalid values before source reconstruction', async () => {
  const initialCheckpoint = checkpoint(), { run, calls } = await controller();
  for (const maxWakeRecoveries of [-1, 1.5, 5, NaN]) {
    assert.throws(() => run(.3, { initialCheckpoint, maxWakeRecoveries }), /controls/);
  }
  assert.equal(calls.transfers.length, 0);
  const terminal = run(.3, { initialCheckpoint, maxStages: 1 });
  assert.equal(calls.recoveries.length, 0); assert.equal(terminal.wakeRecovery, undefined);
});
