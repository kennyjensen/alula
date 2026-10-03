import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { coupledStartupExtensionPlan } from '../src/euler/streamtube-coupled-startup.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

const saved = JSON.parse(gunzipSync(fs.readFileSync(
  new URL('./fixtures/rae64x24-early-trip-extension.json.gz', import.meta.url))));

test('reported 64x24 early-trip endpoint earns continuation of the exact retained equations', () => {
  const { source } = saved, cp = source.checkpoint, before = structuredClone(cp);
  const replay = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 0 });
  assert.deepEqual(replay.families, source.families);
  assert.equal(replay.converged, false);
  assert.ok(replay.families.boundaryLayer > .06);
  Object.assign(replay, { history: source.history, initialRedistribution: source.initialRedistribution });
  const plan = coupledStartupExtensionPlan(replay, { startupAttempt: 1, maxIterations: 40,
    requestedConditions: { mach: .74, reynolds: cp.restart.options.reynolds, ncrit: 4,
      transitionMode: 'automatic', tripFractions: [[.03, .07]], ismom: 4 } });
  assert.ok(plan);
  assert.equal(plan.additionalIterations, 20);
  assert.equal(plan.equationsChanged, false);
  assert.match(plan.selection, /measured descent/);
  const plain = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
  assert.deepEqual(plain(plan.resume), cp);
  assert.deepEqual(cp, before);
  assert.equal(replay.mesh.cells.length, 9072);
  assert.equal(replay.mesh.quality.valid, true);
});

test('the cold worker endpoint after 53 updates also reaches the strict equation root', () => {
  const cp = saved.converged.checkpoint;
  assert.equal(saved.converged.iteration, 53);
  assert.equal(cp.convergence.converged, true);
  const result = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp,
    convergence: 'residual', maxIterations: 3, tolerance: 1e-10 });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.residualConverged, true);
  assert.ok(Math.max(...Object.values(result.families)) <= 1e-10);
  assert.equal(result.mesh.cells.length, 9072);
  assert.equal(result.mesh.quality.valid, true);
  assert.deepEqual(result.checkpoint.restart.input, cp.restart.input);
  assert.deepEqual(result.coupledOptions.tripFractions, [[.03, .07]]);
  assert.equal(result.conditions.mach, .74);
  assert.equal(result.conditions.reynolds, cp.restart.options.reynolds);
  assert.equal(result.conditions.ncrit, 4);
});
