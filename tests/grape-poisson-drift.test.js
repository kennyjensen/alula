import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { grapePoissonDrift } from '../src/geometry/grape-poisson-drift.js';

const close = (a, b) => assert.ok(Math.abs(a - b) < 3e-14 * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
test('GRAPE drift reproduces the executed original Fortran coefficient loop', () => {
  const reference = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/grape-drift.json', import.meta.url)));
  assert.equal(createHash('sha256').update(fs.readFileSync(reference.source)).digest('hex'), reference.sourceSha256);
  for (const { input: [p, r, j2, rdx, alpha, gamma, rdesq], expected } of reference.cases) {
    const { coefficients } = grapePoissonDrift({ drift: j2 * (p + r), leftDistance: 1 / rdx, rightDistance: 1 / rdx });
    const base = [alpha * rdx ** 2, -2 * alpha * rdx ** 2 - 2 * gamma * rdesq, alpha * rdx ** 2];
    expected.forEach((v, k) => close(base[k] + coefficients[k], v));
  }
});

test('nonuniform GRAPE source retains affine consistency and favorable signs at arbitrarily strong drift', () => {
  for (const a of [.001, .1, .4]) for (const b of [.002, .2, .6]) for (const drift of [-1e6, -3, 0, 3, 1e6]) {
    const { coefficients: [l, d, u], first } = grapePoissonDrift({ drift, leftDistance: a, rightDistance: b });
    assert.ok(l >= 0 && u >= 0 && d <= 0); close(l + d + u, 0);
    close(l * -a + u * b, drift); close(first[0] * -a + first[2] * b, 1);
    const reflected = grapePoissonDrift({ drift: -drift, leftDistance: b, rightDistance: a }).coefficients;
    [u, d, l].forEach((v, k) => close(v, reflected[k]));
  }
  for (const input of [{ drift: NaN, leftDistance: 1, rightDistance: 1 }, { drift: 1, leftDistance: 0, rightDistance: 1 }])
    assert.throws(() => grapePoissonDrift(input), /GRAPE drift/);
});

test('the source derivative is explicitly first order and selects a one-sided derivative at zero drift', () => {
  const errors = [16, 32, 64].map(n => {
    const h = 1 / n, { coefficients: [l, , u] } = grapePoissonDrift({ drift: 2, leftDistance: h, rightDistance: h });
    return Math.abs(l * Math.expm1(-h) + u * Math.expm1(h) - 2);
  });
  assert.ok(errors[1] < .51 * errors[0] && errors[2] < .51 * errors[1]);
  assert.ok(errors[2] > .49 * errors[1]);
  assert.deepEqual(grapePoissonDrift({ drift: 0, leftDistance: .1, rightDistance: .2 }).first, [0, -5, 5]);
});
