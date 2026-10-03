import test from 'node:test';
import assert from 'node:assert/strict';
import { monotoneRationalQuadratic as value, inverseMonotoneRationalQuadratic as inverse } from '../src/numerics/rational-quadratic.js';

test('rational Hermite interpolation reproduces affine and quadratic mass fields including zero stagnation slopes', () => {
  for (const a of [0, .3, 1, 1.7, 2]) for (let k = 0; k <= 100; k++) {
    const t = k / 100, y = a * t + (1 - a) * t * t;
    assert.ok(Math.abs(value(t, a, 2 - a) - y) < 4e-16);
    assert.ok(Math.abs(inverse(y, a, 2 - a) - t) < 2e-14);
  }
  assert.equal(inverse(.25, 0, 2), .5);
});

test('nonnegative slopes give bounded monotone segments and a stable inverse over large slope ratios', () => {
  for (const a of [0, 1e-4, .5, 1, 10, 1e4]) for (const b of [0, 1e-4, .5, 1, 10, 1e4]) {
    let before = -1;
    for (let k = 0; k <= 40; k++) {
      const t = k / 40, y = value(t, a, b);
      assert.ok(y >= 0 && y <= 1 && y > before); before = y;
      assert.ok(Math.abs(inverse(y, a, b) - t) < 3e-10, `${a},${b},${t}`);
    }
  }
  const t = inverse(.01, 1e308, 0); assert.ok(t > 0 && t < 1e-308);
  assert.ok(Math.abs(value(t, 1e308, 0) - .01) < 1e-12);
});

test('invalid or unrepresentable interpolation requests are rejected', () => {
  for (const fn of [value, inverse]) for (const args of [[-.1, 1, 1], [1.1, 1, 1], [.5, -1, 1], [.5, 1, Infinity], [NaN, 1, 1]])
    assert.throws(() => fn(...args), /Invalid/);
  assert.throws(() => inverse(.5, 0, 1e308), /Unresolved/);
});
