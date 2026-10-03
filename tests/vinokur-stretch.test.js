// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createVinokurStretch } from '../src/numerics/vinokur-stretch.js';

function close(actual, expected, relative = 8e-13, absolute = 0) {
  assert.ok(Math.abs(actual - expected) <= absolute + relative * Math.abs(expected), `${actual} != ${expected}`);
}
function prescribed(a, branch, length = 7, logicalWidth = 3) {
  const ratio = branch === 'hyperbolic' ? Math.sinh(2 * a) / (2 * a) : Math.sin(2 * a) / (2 * a);
  return { length, logicalWidth, initialDerivative: length / logicalWidth / ratio };
}

test('Vinokur known hyperbolic and trigonometric parameters match the published forward maps', () => {
  for (const [branch, parameters] of [['hyperbolic', [0.03, 0.5, 2, 5]], ['trigonometric', [0.03, 0.5, 1, 1.5]]]) {
    for (const a of parameters) {
      const input = prescribed(a, branch), map = createVinokurStretch(input);
      assert.equal(map.branch, branch);
      close(map.parameter, a, 2e-12);
      assert.equal(map.value(0), 0);
      assert.equal(map.value(1), input.length);
      assert.equal(map.evaluate(0).derivative, input.initialDerivative);
      assert.equal(map.evaluate(1).secondDerivative, 0);
      const tangent = branch === 'hyperbolic' ? Math.tanh : Math.tan;
      for (const t of [0.1, 0.3, 0.8]) {
        const q = map.evaluate(t), u = a * (1 - t);
        close(q.value, input.length * (1 - tangent(u) / tangent(a)), 2e-12);
        const cosine = branch === 'hyperbolic' ? Math.cosh : Math.cos;
        close(q.derivative, input.length / input.logicalWidth * a / tangent(a) / cosine(u) ** 2);
        close(q.secondDerivative, (branch === 'hyperbolic' ? 1 : -1) * 2 * a / input.logicalWidth * q.derivative * tangent(u));
        assert.ok(q.derivative > 0);
      }
      const samples = Array.from({ length: 21 }, (_, i) => map.value(i / 20));
      assert.ok(samples.slice(1).every((v, i) => v > samples[i]));
    }
  }
});

test('Vinokur derivatives use the logical coordinate and agree with independent finite differences', () => {
  for (const branch of ['hyperbolic', 'trigonometric']) {
    const input = prescribed(0.9, branch, 5, 2), map = createVinokurStretch(input), h = 2e-5;
    const physical = x => map.value(x / input.logicalWidth);
    for (const t of [0.15, 0.5, 0.85]) {
      const x = t * input.logicalWidth, q = map.evaluate(t);
      close((physical(x + h) - physical(x - h)) / (2 * h), q.derivative, 2e-9);
      close((physical(x + h) - 2 * physical(x) + physical(x - h)) / h ** 2, q.secondDerivative, 2e-5);
    }
    // Independent backward derivative at the far end has no linear term:
    // f'(1-h)-f'(1)=O(h^2), because the specified far curvature is zero.
    const end = map.evaluate(1).derivative;
    const coarse = Math.abs(map.evaluate(0.99).derivative - end);
    const fine = Math.abs(map.evaluate(0.995).derivative - end);
    close(coarse / fine, 4, 1e-4);
  }
});

test('Vinokur affine and one-ULP near-uniform data preserve the supplied slope', () => {
  const affine = createVinokurStretch({ length: 12, logicalWidth: 3, initialDerivative: 4 });
  assert.equal(affine.branch, 'affine'); assert.equal(affine.parameter, 0);
  assert.deepEqual(affine.evaluate(0.25), { value: 3, derivative: 4, secondDerivative: 0 });
  assert.equal(affine.value(0.3), 12 * 0.3);
  for (const change of [-Number.EPSILON, Number.EPSILON, -1e-8, 1e-8]) {
    const m = 1 + change, map = createVinokurStretch({ length: 1, logicalWidth: 1, initialDerivative: m });
    assert.equal(map.branch, change < 0 ? 'hyperbolic' : 'trigonometric');
    assert.equal(map.evaluate(0).derivative, m);
    assert.ok(map.parameter > 0);
    // log(sinh(2a)/(2a)) or -log(sin(2a)/(2a)) = 2a^2/3 + O(a^4).
    close(map.parameter ** 2, 1.5 * Math.abs(Math.log1p(m - 1)), 1e-8);
    close(map.value(0.4), 0.4, 2e-8);
    assert.equal(map.evaluate(1).secondDerivative, 0);
  }
});

test('Vinokur unit scaling and reflection preserve the endpoint constraints', () => {
  for (const branch of ['hyperbolic', 'trigonometric']) {
    const input = prescribed(0.8, branch), base = createVinokurStretch(input);
    const scaled = createVinokurStretch({ length: 5 * input.length, logicalWidth: 7 * input.logicalWidth,
      initialDerivative: 5 / 7 * input.initialDerivative });
    close(scaled.parameter, base.parameter);
    const a = base.evaluate(0.3), b = scaled.evaluate(0.3);
    close(b.value, 5 * a.value); close(b.derivative, 5 / 7 * a.derivative);
    close(b.secondDerivative, 5 / 49 * a.secondDerivative);
    const reflected = t => {
      const q = base.evaluate(1 - t);
      return { value: input.length - q.value, derivative: q.derivative, secondDerivative: -q.secondDerivative };
    };
    assert.equal(reflected(0).value, 0); assert.equal(reflected(1).value, input.length);
    assert.equal(reflected(1).derivative, input.initialDerivative);
    assert.ok(reflected(0).secondDerivative === 0);
    close(reflected(0.7).value, input.length - a.value);
  }
});

test('Vinokur hyperbolic evaluation avoids large sinh and endpoint cancellation', () => {
  const input = prescribed(30, 'hyperbolic', 2, 3), map = createVinokurStretch(input);
  close(map.parameter, 30);
  close(map.evaluate(1).derivative, 20);
  // At small t, f(t)=f'(0)*t to relative O(a*t), without 1-tanh()/tanh().
  close(map.value(1e-12), input.initialDerivative * input.logicalWidth * 1e-12, 1e-10);
  assert.ok(map.value(1e-12) > 0);
});

test('Vinokur invalid and nonrepresentable requests fail explicitly without changing inputs', () => {
  const input = Object.freeze({ length: 3, logicalWidth: 2, initialDerivative: 0.8 });
  const map = createVinokurStretch(input);
  assert.deepEqual(input, { length: 3, logicalWidth: 2, initialDerivative: 0.8 });
  for (const key of Object.keys(input)) for (const bad of [0, -1, Infinity, NaN, undefined, '2'])
    assert.throws(() => createVinokurStretch({ ...input, [key]: bad }), /finite and positive/);
  assert.throws(() => createVinokurStretch(), /finite and positive/);
  assert.throws(() => createVinokurStretch({ length: Number.MAX_VALUE, logicalWidth: Number.MIN_VALUE, initialDerivative: 1 }), /not representable/);
  assert.throws(() => createVinokurStretch({ length: 1, logicalWidth: 1, initialDerivative: 1e17 }), /representable parameter/);
  assert.throws(() => createVinokurStretch({ length: 1, logicalWidth: 1e-200, initialDerivative: 5e199 }), /second derivative is not representable/);
  for (const t of [-0.1, 1.1, NaN, Infinity, undefined]) {
    assert.throws(() => map.value(t), /\[0,1\]/);
    assert.throws(() => map.evaluate(t), /\[0,1\]/);
  }
});
