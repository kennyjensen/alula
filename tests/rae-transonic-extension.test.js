import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { coupledStartupExtensionPlan } from '../src/euler/streamtube-coupled-startup.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

test('reported transonic grid finishes the exact requested equations after its productive iteration-40 handoff', () => {
  const { original } = JSON.parse(gunzipSync(fs.readFileSync(
    new URL('./fixtures/rae64x11-transonic-extension.json.gz', import.meta.url))));
  const cp = original.checkpoint, before = structuredClone(cp), options = cp.restart.options;
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = cp.continuation;
  const controls = { iterationGeometry, stepAcceptance, stagnationLimiter };
  const replay = solveCoupledStreamtubeIses(undefined, { ...controls, resume: cp, maxIterations: 0, tolerance: 1e-10 });
  assert.deepEqual(replay.families, original.families);
  replay.history = original.history;
  const plan = coupledStartupExtensionPlan(replay, { startupAttempt: 1, maxIterations: 40,
    requestedConditions: { mach: .74, reynolds: options.reynolds, ncrit: 4,
      transitionMode: 'automatic', tripFractions: [[1, 1]], ismom: 4 } });
  assert.ok(plan);
  assert.match(plan.selection, /two full Armijo steps/);
  assert(replay.history[40].boundaryLayer > replay.history[39].boundaryLayer);
  const result = solveCoupledStreamtubeIses(undefined, { ...controls, resume: plan.resume,
    maxIterations: plan.additionalIterations, tolerance: 1e-10 });
  assert.equal(result.converged, true);
  assert.equal(result.reason, 'residual');
  assert(result.history.at(-1).iteration <= plan.additionalIterations);
  assert(Math.max(...Object.values(result.families)) <= 1e-10);
  assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.mesh.cells.length, 5292);
  assert(result.mesh.quality.minCornerSine > .18);
  assert.equal(result.conditions.mach, .74);
  assert.equal(result.conditions.ncrit, 4);
  assert.equal(result.conditions.reynolds, options.reynolds);
  assert.equal(result.solverInput.hybrid.ismom, 4);
  assert.deepEqual(result.solverInput.upwind, { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } });
  assert.deepEqual(cp, before);
});

test('reported 64x7 transonic case retains productive progress beyond the old cap and finishes the requested equations', () => {
  const saved = JSON.parse(gunzipSync(fs.readFileSync(
    new URL('./fixtures/rae64x7-progress-limit.json.gz', import.meta.url))));
  const before = structuredClone(saved);
  const plain = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
  const planFor = packet => {
    const cp = packet.checkpoint, options = cp.restart.options;
    const replay = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp,
      maxIterations: 0, tolerance: 1e-10 });
    assert.deepEqual(replay.families, cp.families);
    Object.assign(replay, { history: packet.history, initialRedistribution: packet.initialRedistribution,
      startupExtension: packet.startupExtension });
    const plan = coupledStartupExtensionPlan(replay, { startupAttempt: 1, maxIterations: 40,
      extensionCount: packet.startupExtension.extensionCount,
      requestedConditions: { mach: .74, reynolds: options.reynolds, ncrit: 4,
        transitionMode: 'automatic', tripFractions: [[1, 1]], ismom: 4 } });
    assert.ok(plan, `Productive iteration ${packet.history.length - 1} must receive another bounded chunk.`);
    assert.equal(plan.additionalIterations, 20);
    assert.equal(plan.equationsChanged, false);
    assert.deepEqual(plain(plan.resume), cp);
    return plan;
  };
  assert.equal(saved.source.history.length, 121);
  planFor(saved.source);
  // The full cold-worker receipt checks every intermediate chunk. Keep the
  // numerical regression short by replaying the final productive checkpoint.
  const plan = planFor(saved.late), cp = plan.resume;
  const r = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp,
    maxIterations: plan.additionalIterations, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason);
  assert.equal(r.reason, 'residual');
  assert(Math.max(...Object.values(r.families)) <= 1e-10);
  assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.mesh.cells.length, 3360);
  assert(r.mesh.quality.minCornerSine > .15);
  assert.deepEqual(r.checkpoint.restart.input, cp.restart.input);
  assert.equal(r.conditions.mach, .74);
  assert.equal(r.conditions.reynolds, cp.restart.options.reynolds);
  assert.equal(r.conditions.ncrit, 4);
  assert.deepEqual(r.coupledOptions.tripFractions, [[1, 1]]);
  assert.deepEqual(saved, before);
});
