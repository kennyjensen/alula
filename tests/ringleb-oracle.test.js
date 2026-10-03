// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { ringlebPoint, ringlebLocalTube, ringlebAtPoint } from './oracles/ringleb.js';

const close = (a, b, tol) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const sonic = Math.sqrt(5 / 6);
const points = [[.82, .63], [sonic, .7], [1.02, .77]];
test('Ringleb map derivatives match independent fourth-order differences across sonic speed', () => {
  for (const [q, psi] of points) for (const axis of ['q', 'psi']) for (const h of [2e-4, 1e-4]) {
    const at = offset => ringlebPoint(q + (axis === 'q' ? offset * h : 0), psi + (axis === 'psi' ? offset * h : 0));
    const [mm, m, p, pp] = [-2, -1, 1, 2].map(at), s = ringlebPoint(q, psi);
    for (const key of ['x', 'y', 'rho', 'p', 'u', 'v'])
      close(s.derivatives[key][axis], (mm[key] - 8 * m[key] + 8 * p[key] - pp[key]) / (12 * h), 3e-10);
  }
});
test('Exact Ringleb field satisfies continuity, both Euler momenta and irrotationality', () => {
  for (const [q, psi] of points) {
    const s = ringlebPoint(q, psi), d = s.derivatives;
    const [dr, du, dv, dp] = ['rho', 'u', 'v', 'p'].map(key => s.spatialGradient(d[key]));
    close(s.u * dr.x + s.v * dr.y + s.rho * (du.x + dv.y), 0, 3e-14);
    close(s.rho * (s.u * du.x + s.v * du.y) + dp.x, 0, 3e-14);
    close(s.rho * (s.u * dv.x + s.v * dv.y) + dp.y, 0, 3e-14);
    close(dv.x - du.y, 0, 3e-14);
    close(s.enthalpy + .5 * (s.u * s.u + s.v * s.v), 2.5, 2e-15);
    close(s.entropyOverR, 0, 2e-15);
  }
});
test('Streamline orientation and exact tube mass use psi, including at Mach one', () => {
  for (const [q, psi] of points) {
    const s = ringlebPoint(q, psi), d = s.derivatives;
    close(s.u * d.y.q - s.v * d.x.q, 0, 2e-15);
    close(s.rho * (s.u * d.y.psi - s.v * d.x.psi), 1, 2e-15);
    assert.ok(s.u * d.x.q + s.v * d.y.q < 0, 'q decreases downstream');
    assert.ok(s.determinant < 0 && s.envelopeMargin > .1);
  }
  close(ringlebPoint(sonic, .7).mach, 1, 2e-15);
});
test('Local nested samples retain the exact physical center and halve both widths', () => {
  const qCenter = 1.04 / Math.sqrt(1 + .2 * 1.04 ** 2), center = ringlebPoint(qCenter, .7);
  for (const level of [0, 1, 2]) {
    const h = 2 ** -level, tube = ringlebLocalTube({ qCenter, qHalfStep: .012 * h, psiWidth: .04 * h });
    assert.equal(tube.qCenter, center.q); assert.equal(tube.psiCenter, center.psi);
    assert.equal(tube.massFlow, .04 * h); assert.equal(tube.lower[3].q, qCenter);
    assert.equal(tube.upper[3].q, qCenter); assert.equal(tube.exactSections.length, 4);
  }
  assert.throws(() => ringlebPoint(1.1, .5), /envelope/);
  assert.throws(() => ringlebPoint(1, 1), /upper branch/);
});
test('Bounded inverse hodograph mapping recovers exact physical samples across sonic speed', () => {
  for (const [q, psi] of points) {
    const target = ringlebPoint(q, psi), s = ringlebAtPoint(target.x, target.y, { guess: { q: q + .003, psi: psi - .002 } });
    close(s.q, q, 8e-14); close(s.psi, psi, 8e-14);
    assert.ok(s.mapping.iterations <= 4 && s.mapping.residual <= s.mapping.tolerance);
    assert.ok(s.mapping.maxCondition < 100);
  }
  const s = ringlebPoint(.9, .7);
  assert.throws(() => ringlebAtPoint(s.x, s.y, { guess: { q: .91, psi: .7 }, maxIterations: 0 }), /exceeded/);
});
