import test from 'node:test';
import assert from 'node:assert/strict';
import { createMonotoneCubicMap } from '../src/numerics/monotone-cubic.js';
import { reconcileFacingArcSlopes } from '../src/geometry/facing-arc-spacing.js';

const prescribed = (map, slopes) => createMonotoneCubicMap(map.knots, map.values, { derivatives: 'prescribed', slopes });
const normalized = (map, rank, start, end) => (map.value(rank) - map.value(start)) / (map.value(end) - map.value(start));

test('one facing block gets a shared nonlinear arc correspondence with all physical hits and C1 joins fixed', () => {
  const a = createMonotoneCubicMap([0, 1, 2], [0, .8, 1]), b = createMonotoneCubicMap([1, 2, 3], [0, .5, 1]);
  const old = structuredClone({ a: a.slopes, b: b.slopes });
  assert.ok(Math.abs(normalized(a, 1.3, 1, 2) - normalized(b, 1.3, 1, 2)) > .05);
  const fit = reconcileFacingArcSlopes(a, b), x = prescribed(a, fit.first), y = prescribed(b, fit.second);
  for (const t of [0, .01, .1, .3, .7, .95, 1])
    assert.ok(Math.abs(normalized(x, 1 + t, 1, 2) - normalized(y, 1 + t, 1, 2)) < 1e-14);
  for (const [before, after] of [[a, x], [b, y]]) {
    before.knots.forEach((v, i) => assert.equal(after.value(v), before.values[i]));
    const knot = before.knots[1], h = 1e-9;
    assert.ok(Math.abs(after.evaluate(knot - h).derivative - after.evaluate(knot + h).derivative) < 1e-8);
  }
  assert.deepEqual({ a: a.slopes, b: b.slopes }, old);
  const repeated = reconcileFacingArcSlopes(x, y);
  repeated.first.forEach((s, i) => assert.ok(Math.abs(s - fit.first[i]) < 1e-15));
  repeated.second.forEach((s, i) => assert.ok(Math.abs(s - fit.second[i]) < 1e-15));
});

test('incompatible neighboring arc ratios receive the stated log-slope compromise without zero slopes', () => {
  const a = createMonotoneCubicMap([0, 1, 2], [0, 1, 5]), b = createMonotoneCubicMap([0, 1, 2], [0, 4, 5]);
  const fit = reconcileFacingArcSlopes(a, b), middle = fit.knots[1];
  assert.equal(middle.targetLogRatio, 0);
  assert.equal(fit.first[1], fit.second[1]);
  assert.ok(fit.first.every(s => s > 0) && fit.second.every(s => s > 0));
  const objective = logRatio => middle.incidentLogRatios.reduce((sum, r) => sum + (logRatio - r) ** 2, 0);
  assert.ok(objective(0) < objective(-.1) && objective(0) < objective(.1));
  prescribed(a, fit.first); prescribed(b, fit.second);
});

test('the common constrained adjustment preserves the fitted ratio and the positive Hermite cone', () => {
  const a = createMonotoneCubicMap([0, 1, 2], [0, .01, 1.01]), b = createMonotoneCubicMap([1, 2, 3], [0, 1, 1.01]);
  const fit = reconcileFacingArcSlopes(a, b);
  assert.ok(fit.knots.every(k => k.commonLogAdjustment < 0));
  for (const k of fit.knots) assert.ok(Math.abs(Math.log(k.after[0] / k.after[1]) - k.targetLogRatio) < 1e-14);
  assert.equal(fit.first[1], .02);
  prescribed(a, fit.first); prescribed(b, fit.second);
  const scaledA = createMonotoneCubicMap(a.knots.map(v => 3 + 2 * v), a.values.map(v => -2 + 7 * v));
  const scaledB = createMonotoneCubicMap(b.knots.map(v => 3 + 2 * v), b.values.map(v => 9 + 4 * v));
  const scaled = reconcileFacingArcSlopes(scaledA, scaledB);
  fit.first.forEach((s, i) => assert.ok(Math.abs(scaled.first[i] - 3.5 * s) < 1e-13));
  fit.second.forEach((s, i) => assert.ok(Math.abs(scaled.second[i] - 2 * s) < 1e-13));
});

test('only actual facing blocks are reconciled and missing physical hits fail explicitly', () => {
  const a = createMonotoneCubicMap([0, 1], [0, 2]), b = createMonotoneCubicMap([2, 3], [0, 1]);
  const fit = reconcileFacingArcSlopes(a, b);
  assert.deepEqual(fit.first, a.slopes); assert.deepEqual(fit.second, b.slopes); assert.deepEqual(fit.knots, []);
  const overlap = createMonotoneCubicMap([.5, 1.5], [0, 1]);
  assert.throws(() => reconcileFacingArcSlopes(a, overlap), /shared physical hit/);
  assert.throws(() => reconcileFacingArcSlopes({ ...a, slopes: [0, 1] }, b), /Invalid/);
});
