import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

test('transonic section matching damps finite grid maintenance independently of Newton', () => {
  const cp = JSON.parse(gunzipSync(fs.readFileSync(new URL('./fixtures/rae-mach076-maintenance-stall.json.gz', import.meta.url))));
  const before = structuredClone(cp);
  const r = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp, maxIterations: 12, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason);
  assert.equal(r.conditions.edgeMatching, 'section-velocity');
  assert.equal(r.conditions.mach, .76);
  assert.equal(r.mesh.quality.valid, true);
  const first = r.history[1];
  assert.equal(first.step, 1);
  assert.ok(first.maintenance.correctionScale > 0 && first.maintenance.correctionScale < 1);
  assert.ok(first.residualDecrease.afterSquaredNorm < first.residualDecrease.beforeSquaredNorm);
  assert.ok(r.history.slice(1).every(h => h.step === 1));
  for (const residual of Object.values(r.families)) assert.ok(residual <= 1e-10);
  assert.deepEqual(r.checkpoint.restart.input, cp.restart.input);
  assert.deepEqual(cp, before);
});
