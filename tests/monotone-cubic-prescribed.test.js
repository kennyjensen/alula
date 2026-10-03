import test from 'node:test';
import assert from 'node:assert/strict';
import { createMonotoneCubicMap } from '../src/numerics/monotone-cubic.js';

const close = (a, b, tolerance = 2e-13) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const prescribed = (x, y, slopes) => createMonotoneCubicMap(x, y, { derivatives: 'prescribed', slopes });
const polynomial = x => x + .2 * x * x + .1 * x ** 3;
const first = x => 1 + .4 * x + .3 * x * x;
const second = x => .4 + .6 * x;

test('prescribed Hermite slopes recover an independent cubic, fixed hits, and continuous first derivatives', () => {
  const x = [0, .5, 1, 2], y = x.map(polynomial), slopes = x.map(first), map = prescribed(x, y, slopes);
  x.forEach((knot, i) => {
    assert.equal(map.value(knot), y[i]); close(map.evaluate(knot).derivative, slopes[i]);
  });
  for (const q of [.03, .2, .49, .51, .77, 1.01, 1.7, 2]) {
    const value = map.evaluate(q);
    close(value.value, polynomial(q)); close(value.derivative, first(q)); close(value.secondDerivative, second(q));
  }
  const before = map.evaluate(.73);
  x[1] = .4; y[1] = -1; slopes[1] = -100;
  assert.deepEqual(map.evaluate(.73), before, 'Caller array mutations must not change the interpolation.');
});

test('the entire derivative quadratic remains strictly positive throughout the prescribed two-secant cone', () => {
  for (const a of [.01, .5, 1, 2]) for (const b of [.01, .5, 1, 2]) {
    const map = prescribed([0, 1], [0, 1], [a, b]);
    // q(t)=a+(6-4a-2b)t+3(a+b-2)t^2. Its endpoints and sole
    // possible interior minimum exhaust the interval, without a sample grid.
    const linear = 6 - 4 * a - 2 * b, quadratic = 3 * (a + b - 2), locations = [0, 1];
    if (quadratic > 0) {
      const t = -linear / (2 * quadratic);
      if (t > 0 && t < 1) locations.push(t);
    }
    for (const t of locations) {
      assert.ok(a + linear * t + quadratic * t * t > 0);
      assert.ok(map.evaluate(t).derivative > 0);
    }
  }
  close(prescribed([0, 1], [0, 1], [.5, 2]).value(.5), .3125,
    1e-15); // Two-point prescribed interpolation need not be linear.
});

test('equal facing endpoint slope ratios produce identical nonuniform normalized arc cubics', () => {
  // Arc/rank secants are 2 and 5. Equal normalized endpoint derivatives
  // (.5,1.5) give the nonuniform exact progress law (t+t^2)/2.
  const a = prescribed([0, 3], [1, 7], [1, 3]), b = prescribed([0, 3], [-4, 11], [2.5, 7.5]);
  for (const rank of [0, .13, .7, 1.5, 2.6, 3]) {
    const t = rank / 3, pa = a.evaluate(rank), pb = b.evaluate(rank);
    close((pa.value - 1) / 6, (t + t * t) / 2);
    close((pa.value - 1) / 6, (pb.value + 4) / 15);
    close(pa.derivative / 2, pb.derivative / 5);
    close(pa.secondDerivative / 2, pb.secondDerivative / 5);
  }
});

test('prescribed slopes preserve geometric units and reflection', () => {
  const x = [0, .5, 1, 2], y = x.map(polynomial), slopes = x.map(first);
  const map = prescribed(x, y, slopes);
  const scaled = prescribed(x.map(v => 7 + 4 * v), y.map(v => -3 + 2 * v), slopes.map(v => .5 * v));
  const reflected = prescribed(x.toReversed().map(v => 5 - v), y.toReversed().map(v => 10 - v), slopes.toReversed());
  for (const q of [.13, .8, 1.7]) {
    const base = map.evaluate(q), units = scaled.evaluate(7 + 4 * q), reverse = reflected.evaluate(5 - q);
    close(units.value, -3 + 2 * base.value); close(units.derivative, .5 * base.derivative);
    close(units.secondDerivative, .125 * base.secondDerivative);
    close(reverse.value, 10 - base.value); close(reverse.derivative, base.derivative);
    close(reverse.secondDerivative, -base.secondDerivative);
  }
});

test('prescribed data reject missing, nonpositive, nonfinite, or either-adjacent-cone violations', () => {
  for (const slopes of [undefined, null, [1], new Array(2), [1, 0], [1, -1], [1, Infinity], [NaN, 1], [1, 2.01]])
    assert.throws(() => prescribed([0, 1], [0, 1], slopes), /monotonicity cone/);
  // The interior slope .3 satisfies the left secant's cap but violates
  // the smaller right secant's cap of approximately .2.
  assert.throws(() => prescribed([0, 1, 2], [0, 1, 1.1], [1, .3, .1]), /monotonicity cone/);
  assert.throws(() => prescribed([0, 1, 2], [0, .1, 1.1], [.1, .3, 1]), /monotonicity cone/);
  for (const derivatives of [undefined, 'harmonic', 'pchip'])
    assert.throws(() => createMonotoneCubicMap([0, 1], [0, 1], { derivatives, slopes: [1, 1] }), /require/);
});
