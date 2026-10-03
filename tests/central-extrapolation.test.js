import test from 'node:test';
import assert from 'node:assert/strict';
import { extrapolateCentral } from '../scripts/validation/central-extrapolation.js';

test('eighth-order centered extrapolation differentiates degree-seven polynomials and has the degree-nine remainder', () => {
  const derivative = (h, ninth) => {
    const f = t => 2 * t + 3 * t ** 3 + 5 * t ** 5 + 7 * t ** 7 + ninth * t ** 9;
    const samples = Array.from({ length: 5 }, (_, i) => {
      const step = h / 2 ** i; return Float64Array.of((f(step) - f(-step)) / (2 * step));
    });
    return extrapolateCentral(samples).at(-1).estimates.map(d => d[0]);
  };
  for (const value of derivative(.8, 0)) assert.ok(Math.abs(value - 2) < 2e-14);
  const values = derivative(.8, 11);
  values.forEach((v, i) => assert.ok(Math.abs(v - (2 - 11 * (.8 / 2 ** i) ** 8 / 4096)) < 2e-14));
  assert.ok(Math.abs((2 - values[0]) / (2 - values[1]) - 256) < 2e-6);
});

test('resolved smooth derivatives expose an incorrect Jacobian and preserve unavailable domain samples', () => {
  const f = t => Float64Array.of(Math.exp(t), Math.sin(2 * t), 1 / (1 - t));
  const samples = Array.from({ length: 6 }, (_, i) => {
    const h = .2 / 2 ** i, p = f(h), m = f(-h); return p.map((v, k) => (v - m[k]) / (2 * h));
  });
  const exact = [1, 2, 1], resolved = extrapolateCentral(samples).at(-1).estimates;
  for (const v of resolved) v.forEach((q, k) => assert.ok(Math.abs(q - exact[k]) < 1e-8));
  assert.ok(Math.abs(resolved[0][1] - 2.001) > 5e-6);
  const missing = extrapolateCentral([undefined, ...samples]).at(-1).estimates;
  assert.equal(missing[0], undefined); assert.deepEqual(missing.slice(1), resolved);
});
