import test from 'node:test';
import assert from 'node:assert/strict';
import { createMonotoneCubicMap } from '../src/numerics/monotone-cubic.js';

const close = (a, b, tolerance = 2e-13) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const pchip = (x, y) => createMonotoneCubicMap(x, y, { derivatives: 'pchip' });

test('PCHIP recovers three-point polynomial endpoint derivatives and independent Hermite values and curvature', () => {
  // y=x^2 at x=0,1,3. Three-point endpoint derivatives are exactly 0,6;
  // the weighted harmonic interior derivative is 3/2 rather than y'(1).
  const map = pchip([0, 1, 3], [0, 1, 9]);
  assert.equal(map.slopes[0], 0); close(map.slopes[1], 1.5); assert.equal(map.slopes[2], 6);
  for (const x of [0, .07, .4, .93]) {
    const actual = map.evaluate(x);
    close(actual.value, 1.5 * x * x - .5 * x ** 3);
    close(actual.derivative, 3 * x - 1.5 * x * x);
    close(actual.secondDerivative, 3 - 3 * x);
  }
  for (const x of [1, 1.13, 1.7, 2.91, 3]) {
    const v = x - 1, actual = map.evaluate(x);
    close(actual.value, 1 + 1.5 * v + 1.5 * v * v - .125 * v ** 3);
    close(actual.derivative, 1.5 + 3 * v - .375 * v * v);
    close(actual.secondDerivative, 3 - .75 * v);
  }
  close(map.evaluate(1).secondDerivative, 3); // Interior-knot curvature uses the right interval.
  for (const x of [0, 1, 3]) assert.equal(map.value(x), x * x);
});

test('nonuniform PCHIP derivatives match an independent rational reference while the default harmonic rule is unchanged', () => {
  // hL=1,hR=3 and slopes 2,1/3 give d=12/(7/2+15)=24/37.
  const x = [0, 1, 4], y = [0, 2, 3], map = pchip(x, y);
  close(map.slopes[0], 29 / 12); close(map.slopes[1], 24 / 37); assert.equal(map.slopes[2], 0);
  const legacy = createMonotoneCubicMap(x, y), explicit = createMonotoneCubicMap(x, y, { derivatives: 'harmonic' });
  assert.deepEqual(legacy.slopes, [2, (1 / 3) * (2 / (1 + (1 / 3) / 2)), 1 / 3]);
  assert.deepEqual(explicit.slopes, legacy.slopes);
  for (const q of [0, .1, .5, 1, 1.9, 3.7, 4]) assert.deepEqual(explicit.evaluate(q), legacy.evaluate(q));
  assert.deepEqual(x, [0, 1, 4]); assert.deepEqual(y, [0, 2, 3]);
});

test('PCHIP preserves exact affine and two-point data and transforms correctly under units and reflection', () => {
  const line = pchip([0, .25, 1, 4], [3, 3.5, 5, 11]);
  for (const x of [0, .125, .25, .5, 1, 2, 4])
    assert.deepEqual(line.evaluate(x), { value: 3 + 2 * x, derivative: 2, secondDerivative: 0 });
  const two = pchip([2, 5], [-1, 8]);
  for (const x of [2, 2.125, 3, 4.5, 5]) assert.deepEqual(two.evaluate(x), { value: 3 * x - 7, derivative: 3, secondDerivative: 0 });
  const x = [0, .1, .7, 2, 5], y = [1, 1.2, 3, 4, 9], original = pchip(x, y);
  const scaled = pchip(x.map(v => 7 + 4 * v), y.map(v => -3 + 2 * v));
  const reflected = pchip(x.toReversed().map(v => 10 - v), y.toReversed().map(v => 20 - v));
  for (const q of [.03, .3, 1.2, 3.7]) {
    const value = original.evaluate(q), units = scaled.evaluate(7 + 4 * q), reverse = reflected.evaluate(10 - q);
    close(units.value, -3 + 2 * value.value); close(units.derivative, .5 * value.derivative);
    close(units.secondDerivative, .125 * value.secondDerivative);
    close(reverse.value, 20 - value.value); close(reverse.derivative, value.derivative);
    close(reverse.secondDerivative, -value.secondDerivative);
  }
});

test('PCHIP derivative quadratics are nonnegative throughout strongly nonuniform increasing intervals', () => {
  for (const [x, y] of [[[0, .001, 1, 1.01, 7], [0, 2, 2.001, 4, 4.1]],
    [[0, 1, 4], [0, 2, 3]], [[0, 1, 2, 5], [0, .001, 100, 101]]]) {
    const map = pchip(x, y);
    for (let i = 0; i < x.length - 1; i++) {
      const h = x[i + 1] - x[i], secant = (y[i + 1] - y[i]) / h;
      const a = map.slopes[i] / secant, b = map.slopes[i + 1] / secant;
      const linear = 6 - 4 * a - 2 * b, quadratic = 3 * (a + b - 2);
      // Endpoints and the quadratic's sole interior minimum suffice; this
      // checks the whole polynomial rather than relying on a sampling grid.
      assert.ok(a >= 0 && b >= 0);
      if (quadratic > 0) {
        const t = -linear / (2 * quadratic);
        if (t > 0 && t < 1) assert.ok(a + linear * t + quadratic * t * t > 0);
      }
      for (const t of [.1, .37, .8]) {
        const actual = map.evaluate(x[i] + t * h);
        assert.ok(actual.value > y[i] && actual.value < y[i + 1]);
        close(actual.derivative, secant * (a + linear * t + quadratic * t * t));
        close(actual.secondDerivative, secant / h * (linear + 2 * quadratic * t));
      }
    }
  }
});

test('PCHIP options and strictly increasing-data requirements reject unsupported inputs', () => {
  for (const option of [null, 'pchip', [], { derivatives: null }, { derivatives: 'unknown' }, { derivative: 'pchip' }])
    assert.throws(() => createMonotoneCubicMap([0, 1], [0, 1], option), /options/);
  for (const [x, y] of [[[0, 0], [0, 1]], [[0, 1], [0, 0]], [[0, 1], [1, 0]], [[0, 1], [0, Infinity]]])
    assert.throws(() => pchip(x, y), /increasing/);
  const map = pchip([0, 1, 2], [0, 1, 4]);
  assert.throws(() => map.evaluate(-.1), /outside/);
  assert.throws(() => map.evaluate(NaN), /outside/);
});
