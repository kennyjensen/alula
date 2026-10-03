import test from 'node:test';
import assert from 'node:assert/strict';
import { boundedResidualPartial } from '../src/numerics/bounded-partial.js';

test('bounded residual partials reproduce independent polynomial derivatives without leaving their domain', () => {
  for (const x of [0, 1e-12, .4, 1 - 1e-12, 1]) {
    const probes = [], f = y => { assert.ok(y >= 0 && y <= 1); probes.push(y); return [y * y, 3 * y - 7, 1]; };
    const d = boundedResidualPartial(f, x, { step: 1e-4, lower: 0, upper: 1 });
    for (const [i, expected] of [2 * x, 3, 0].entries()) assert.ok(Math.abs(d[i] - expected) < 1e-10);
    assert.ok(probes.length <= 3);
  }
  assert.throws(() => boundedResidualPartial(x => [x], 2, { step: .1, lower: 0, upper: 1 }), /domain/);
  assert.throws(() => boundedResidualPartial(() => [NaN], 0, { step: .1 }), /Nonfinite/);
});

test('one-sided boundary derivatives retain second-order accuracy under step refinement', () => {
  for (const [x, sign] of [[0, 1], [1, -1]]) {
    const errors = [.02, .01, .005].map(step => Math.abs(boundedResidualPartial(y => [Math.exp(y)], x,
      { step, lower: 0, upper: 1 })[0] - Math.exp(x)));
    assert.ok(errors[0] / errors[1] > 3.8); assert.ok(errors[1] / errors[2] > 3.8);
  }
});
