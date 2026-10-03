import test from 'node:test';
import assert from 'node:assert/strict';
import { createMonotoneCubicMap } from '../src/numerics/monotone-cubic.js';

test('monotone block interpolation preserves anchors and affine data and joins derivatives continuously', () => {
  const x = [0, .2, 1, 4], y = x.map(v => 3 + 2 * v), map = createMonotoneCubicMap(x, y);
  for (const v of x) assert.equal(map.value(v), 3 + 2 * v);
  for (const v of [.03, .15, .7, 2, 3.8]) {
    assert.ok(Math.abs(map.value(v) - 3 - 2 * v) < 2e-15);
    assert.ok(Math.abs(map.evaluate(v).derivative - 2) < 2e-14);
  }
  const curved = createMonotoneCubicMap([0, 2, 3, 8], [0, .1, 1, 2]);
  for (let i = 1; i < 3; i++) {
    const v = curved.knots[i];
    assert.ok(Math.abs(curved.evaluate(v - 1e-8).derivative - curved.evaluate(v + 1e-8).derivative) < 1e-7);
    assert.equal(curved.evaluate(v).derivative, curved.slopes[i]);
  }
  assert.deepEqual(x, [0, .2, 1, 4]); assert.deepEqual(y, [3, 3.4, 5, 11]);
});

test('strongly unequal block slopes have positive derivatives throughout every interval without overshoot', () => {
  const x = [0, .001, 1, 1.01, 7], y = [0, 2, 2.001, 4, 4.1], map = createMonotoneCubicMap(x, y);
  for (let i = 0; i < x.length - 1; i++) {
    const secant = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    const a = map.slopes[i] / secant, b = map.slopes[i + 1] / secant;
    // Independent minimum of the exact quadratic derivative on [0,1].
    const c1 = 6 - 4 * a - 2 * b, c2 = 3 * (a + b - 2);
    const probes = [0, 1]; if (c2 > 0 && -c1 / (2 * c2) > 0 && -c1 / (2 * c2) < 1) probes.push(-c1 / (2 * c2));
    for (const t of probes) assert.ok(a + c1 * t + c2 * t * t > 0);
    let previous = y[i];
    for (let k = 1; k <= 100; k++) {
      const p = map.evaluate(k === 100 ? x[i + 1] : x[i] + (x[i + 1] - x[i]) * k / 100);
      assert.ok(p.value > previous && p.value <= y[i + 1]); assert.ok(p.derivative > 0); previous = p.value;
    }
  }
  const scaled = createMonotoneCubicMap(x.map(v => 3 + 2 * v), y.map(v => -5 + .7 * v));
  for (const v of [.0004, .3, 1.005, 2, 5]) assert.ok(Math.abs(scaled.value(3 + 2 * v) - (-5 + .7 * map.value(v))) < 2e-12);
  for (const args of [[[0, 0], [0, 1]], [[0, 1], [1, 0]], [[0, 1], [0, Infinity]]])
    assert.throws(() => createMonotoneCubicMap(...args), /increasing/);
  assert.throws(() => map.value(-1), /outside/);
});
