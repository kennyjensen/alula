import test from 'node:test';
import assert from 'node:assert/strict';
import { distributeFacingDensity } from '../src/geometry/facing-density-stations.js';
import { createSpacingEnvelope } from '../src/geometry/spacing-envelope.js';

const uniform = n => Array.from({ length: n + 1 }, (_, i) => i / n);
test('the finer uniform request sets the common count in physical-arc units', () => {
  const profiles = [uniform(8), uniform(16).map(u => 3 + .3 * u)], r = distributeFacingDensity({ profiles });
  assert.equal(r.intervals, 16); assert.equal(r.resolutionSatisfied, true);
  r.progress.forEach((u, i) => assert.ok(Math.abs(u - i / 16) < 2e-14));
  r.positions[1].forEach((s, i) => assert.ok(Math.abs(s - (3 + .3 * i / 16)) < 1e-14));
});

test('complementary end refinement requires more intervals than either input alone', () => {
  const profiles = [[0, .1, .2, .6, 1], [0, .4, .8, .9, 1]];
  // Analytic lower envelope: two constant end segments, two ramps and
  // a constant central segment. This is independent of the implementation.
  const expected = 2 * (.15 / .1 + Math.log(4) / 1.2) + .2 / .4;
  const r = distributeFacingDensity({ profiles });
  assert.ok(Math.abs(r.requestedIntervals - expected) < 1e-13);
  assert.equal(r.requiredIntervals, 6); assert.equal(r.intervals, 6);
  const fixed = distributeFacingDensity({ profiles, intervals: 4 });
  assert.equal(fixed.resolutionSatisfied, false); assert.ok(fixed.densityScale < 1);
  assert.equal(fixed.endpointJoinMatched, false);
  r.progress.forEach((u, i) => assert.ok(Math.abs(u + r.progress[6 - i] - 1) < 2e-15));
});

test('face ordering, physical units, origins and duplicate requests leave the result unchanged', () => {
  const profiles = [[0, .02, .1, .3, .6, 1], [0, .3, .7, .9, .98, 1]];
  const a = distributeFacingDensity({ profiles, intervals: 20 });
  const b = distributeFacingDensity({ profiles: profiles.map((p, k) => p.map(s => s * (k ? 2 : .3) + k)).reverse(), intervals: 20 });
  a.progress.forEach((u, i) => assert.ok(Math.abs(u - b.progress[i]) < 2e-14));
  assert.deepEqual(distributeFacingDensity({ profiles: [...profiles, profiles[0]], intervals: 20 }).progress, a.progress);
});

test('envelope mass inversion and increasing stations hold across an actual request crossover', () => {
  const profiles = [[0, .1, .4, .7, 1], [0, .3, .55, .8, 1]], e = createSpacingEnvelope({ profiles, start: 0, end: 1 });
  const r = distributeFacingDensity({ profiles, intervals: 51 });
  r.progress.forEach((u, i) => assert.ok(Math.abs(e.metric(u) - i * e.total / 51) < 2e-14));
  for (let i = 0; i <= 100; i++) assert.ok(Math.abs(e.inverse(e.metric(i / 100)) - i / 100) < 1e-14);
  assert.throws(() => e.metric(-.1), /bounds/); assert.throws(() => e.inverse(e.total + 1), /bounds/);
});

test('bad ordering and insufficient budgets cannot silently change the request', () => {
  assert.throws(() => distributeFacingDensity({ profiles: [[0, .5, .5, 1]] }), /ordered/);
  assert.throws(() => distributeFacingDensity({ profiles: [uniform(16)], maximumIntervals: 8 }), /budget/);
  assert.throws(() => distributeFacingDensity({ profiles: [uniform(8)], intervals: 2.5 }), /counts/);
});

test('a physical hit clips the original monitor without inventing refined endpoint intervals', () => {
  const r = distributeFacingDensity({ profiles: [uniform(8)], ranges: [[.237, .863]] });
  assert.ok(Math.abs(r.requestedIntervals - 8 * (.863 - .237)) < 1e-14);
  assert.equal(r.requiredIntervals, 6);
  assert.equal(r.positions[0][0], .237); assert.equal(r.positions[0].at(-1), .863);
  const full = createSpacingEnvelope({ profiles: [uniform(8)], start: 0, end: 1 });
  const part = createSpacingEnvelope({ profiles: [uniform(8)], start: .237, end: .863 });
  assert.ok(Math.abs(part.total - (full.metric(.863) - full.metric(.237))) < 1e-14);
  assert.throws(() => distributeFacingDensity({ profiles: [uniform(8)], ranges: [[-.1, .8]] }), /ranges/);
});

test('endpoint fitting preserves common physical progress and reports local loss of requested resolution', () => {
  const profiles = [uniform(16), uniform(8).map(s => 2 + .3 * s)];
  const r = distributeFacingDensity({ profiles, intervals: 16, endpointSpacings: { firstSpacing: .12, lastSpacing: .12 } });
  assert.ok(Math.abs(r.progress[1] - .12) < 1e-14);
  assert.ok(Math.abs(1 - r.progress.at(-2) - .12) < 1e-14);
  r.positions[1].forEach((s, i) => assert.ok(Math.abs((s - 2) / .3 - r.progress[i]) < 3e-15));
  assert.equal(r.countCapacitySatisfied, true);
  assert.equal(r.resolutionSatisfied, false);
  assert.ok(r.maximumDensityInterval >= 1.92 - 1e-13);
  assert.throws(() => distributeFacingDensity({ profiles, intervals: 16, endpointSpacings: { firstSpacing: .6, lastSpacing: .6 } }), /feasible/);
});
