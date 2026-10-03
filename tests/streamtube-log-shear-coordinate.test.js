import test from 'node:test';
import assert from 'node:assert/strict';
import { logarithmicShearIncrement as advance } from '../src/euler/streamtube-coupled-xfoil-update.js';

test('logarithmic shear has the physical Newton tangent across positive shear scales', () => {
  for (const current of [1e-12, 1e-8, 1e-4, .03, .24])
    for (const relativeDirection of [-114.47366889630831, -.01, .1, 5]) {
      const direction = current * relativeDirection;
      const error = step => Math.abs((advance(current, direction, step) - current) / step - direction)
        / Math.max(current, Math.abs(direction));
      assert(error(1e-6) < 6e-5);
      assert(error(5e-7) < 3e-5);
    }
});

test('logarithmic coordinate scaling preserves a nonzero shear-equation Newton correction', () => {
  // This transport-shaped row is deliberately evaluated away from its root.
  const current = .002, upstream = .03, A = .2, B = .7, C = .13, upstreamDirection = -.004;
  const equation = (shear, previous) => A - B * shear - C * Math.log(shear / previous);
  const residual = equation(current, upstream), derivative = -B - C / current, upstreamDerivative = C / upstream;
  const direction = (-residual - upstreamDerivative * upstreamDirection) / derivative;
  assert(Math.abs(residual) > .1);
  assert(Math.abs(derivative * current * (direction / current) + upstreamDerivative * upstreamDirection + residual) < 1e-15);
  for (const step of [1e-5, 1e-6]) {
    const observed = (equation(advance(current, direction, step), upstream + step * upstreamDirection) - residual) / step;
    assert(Math.abs(observed + residual) < 1e-5);
  }
});

test('zero logarithmic increments are exact and unusable positive-domain values reject', () => {
  const current = .00014339274139684038;
  assert.equal(advance(current, -1, 0), current);
  assert.equal(advance(current, 0, 1), current);
  const direction = -.016414693200795867, step = .02;
  assert(current + step * direction < 0);
  assert(advance(current, direction, step) > 0);
  for (const values of [[0, 1, .1], [-1, 1, .1], [1, Infinity, .1], [1, 1, NaN],
    [1, 1, -1], [1, 1, 2], [1, 1000, 1], [1, -1000, 1]])
    assert.throws(() => advance(...values));
});

test('subnormal directions and finite results with extreme exponent intermediates remain representable', () => {
  assert.equal(advance(Number.MIN_VALUE, 2 * Number.MIN_VALUE, .4), 2 * Number.MIN_VALUE);
  const decay = advance(1e300, -1e303, .75);
  assert(decay > 0 && Number.isFinite(decay));
  assert(Math.abs(Math.log(decay) - (Math.log(1e300) - 750)) < 2e-13);
  const growth = advance(1e-300, 750e-300, 1);
  assert(growth > 0 && Number.isFinite(growth));
  assert(Math.abs(Math.log(growth) - (Math.log(1e-300) + 750)) < 2e-13);
});
