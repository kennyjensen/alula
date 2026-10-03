// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createXfoilDeadAirGap } from '../src/viscous/xfoil-dead-air-gap.js';

const close = (a, b, tol = 2e-13) => assert.ok(Math.abs(a - b) <= tol, `${a} != ${b}`);
const geometry = (slope, normalGap = .002) => {
  // With upper=(-1,0), choose lower so its normalized cross is
  // slope/sqrt(1+slope^2), hence the uncapped XICALC slope is known.
  const c = slope / Math.sqrt(1 + slope * slope);
  return { normalGap, upperDerivative: { x: -1, y: 0 }, lowerDerivative: { x: Math.sqrt(1 - c * c), y: -c } };
};

test('XICALC cubic matches its TE width/slope and closes C1 at 2.5 widths', () => {
  for (const slope of [-2, -.3, 0, .7, 2]) {
    const m = createXfoilDeadAirGap(geometry(slope)), expectedSlope = Math.max(-1.2, Math.min(1.2, slope));
    close(m.at(0).gap, .002, 2e-18); close(m.at(0).dDistance, expectedSlope, 2e-15);
    assert.equal(m.closureDistance, .005); assert.equal(m.lengthToGapRatio, 2.5);
    assert.deepEqual(m.at(.005), { gap: 0, dDistance: 0 }); assert.deepEqual(m.at(.05), { gap: 0, dDistance: 0 });
    const near = m.at(.005 * (1 - 1e-7)); assert.ok(near.gap >= 0 && near.gap < 2e-15); assert.ok(Math.abs(near.dDistance) < 5e-7);
  }
  const expanding = createXfoilDeadAirGap(geometry(2));
  // Positive limited slope deliberately allows the documented polynomial
  // to expand before closing. No artificial cap at the TE width is allowed.
  close(expanding.at(.001).gap, 1.28 * .002, 2e-18);
});

test('distance derivative agrees with independent fourth-order differences inside the cubic', () => {
  for (const slope of [-3, -.4, 0, .6, 3]) {
    const m = createXfoilDeadAirGap(geometry(slope));
    for (const fraction of [.09, .36, .73, .94]) {
      const x = fraction * m.closureDistance, h = 1e-4 * m.closureDistance;
      const f = d => m.at(d).gap;
      const fd = (f(x - 2 * h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h)) / (12 * h);
      close(m.at(x).dDistance, fd, 3e-12);
    }
  }
});

test('gap scaling and rotation preserve the physical value and dimensionless slope', () => {
  const g = geometry(-.4), a = createXfoilDeadAirGap(g);
  const angle = .72, rot = p => ({ x: p.x * Math.cos(angle) - p.y * Math.sin(angle), y: p.x * Math.sin(angle) + p.y * Math.cos(angle) });
  for (const factor of [.01, 3, 1e3]) {
    const b = createXfoilDeadAirGap({ ...g, normalGap: factor * g.normalGap,
      upperDerivative: rot(g.upperDerivative), lowerDerivative: rot(g.lowerDerivative) });
    for (const x of [0, .0003, .003, .006]) {
      close(b.at(factor * x).gap, factor * a.at(x).gap, 2e-15 * Math.max(1, factor));
      close(b.at(factor * x).dDistance, a.at(x).dDistance, 2e-15);
    }
  }
});

test('zero-gap, sharp and perpendicular derivative limits are explicit and finite', () => {
  for (const g of [{ ...geometry(.3), normalGap: 0 }, { ...geometry(.3), sharp: true }]) {
    const m = createXfoilDeadAirGap(g);
    assert.equal(m.closureDistance, 0);
    for (const x of [0, .001, 1]) assert.deepEqual(m.at(x), { gap: 0, dDistance: 0 });
  }
  for (const sign of [-1, 1]) {
    const m = createXfoilDeadAirGap({ normalGap: .002, upperDerivative: { x: 1, y: 0 }, lowerDerivative: { x: 0, y: sign } });
    assert.equal(m.slope, sign * 1.2); assert.equal(m.rawSlope, null); assert.equal(m.rawSlopeUnbounded, true);
    assert.ok(Object.values(m.at(.001)).every(Number.isFinite));
  }
  assert.throws(() => createXfoilDeadAirGap({ ...geometry(0), normalGap: -1 }), /Invalid/);
  assert.throws(() => createXfoilDeadAirGap({ ...geometry(0), upperDerivative: { x: 0, y: 0 } }), /nonzero/);
  assert.throws(() => createXfoilDeadAirGap(geometry(0)).at(-1), /nonnegative/);
  assert.throws(() => createXfoilDeadAirGap(geometry(0)).at(NaN), /finite/);
});

test('standalone gap values replay the original Fortran XICALC fixture', () => {
  const fixture = JSON.parse(fs.readFileSync('tests/fixtures/fortran/dead-air-gap.json'));
  assert.equal(fixture.passed, true); assert.equal(fixture.cases.length, 3);
  for (const c of fixture.cases) {
    const m = createXfoilDeadAirGap(c.input);
    for (const p of c.samples) close(m.at(p.distance).gap, p.expectedGap, 2e-14 * c.input.normalGap);
  }
});
