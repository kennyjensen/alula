import test from 'node:test';
import assert from 'node:assert/strict';
import { distributeCurvatureSurface } from '../src/geometry/curvature-surface-spacing.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const close = (a, b, tol) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);
const parabola = { length: 2.8, knots: [0, .4, .9, 1.4, 2, 2.8], evaluate: s => ({ point: { x: s, y: s * s },
  derivative: { x: 1, y: 2 * s }, secondDerivative: { x: 0, y: 2 } }) };
const exactArc = s => .5 * s * Math.sqrt(1 + 4 * s * s) + .25 * Math.asinh(2 * s);
const integrate = (curve, a, b, fn) => {
  const knots = [a, ...curve.knots.filter(s => s > a && s < b), b]; let sum = 0;
  for (let k = 1; k < knots.length; k++) {
    const left = knots[k - 1], right = knots[k], h = (right - left) / 64; let subtotal = 0;
    for (let i = 0; i <= 64; i++) subtotal += (i === 0 || i === 64 ? 1 : i % 2 ? 4 : 2) * fn(curve.evaluate(i === 64 ? right : left + i * h));
    sum += subtotal * h / 3;
  }
  return sum;
};
const speed = v => Math.hypot(v.derivative.x, v.derivative.y);

test('curvature monitor equidistribution agrees with analytic parabola arc and tangent-turn integrals on both branches', () => {
  for (const side of ['upper', 'lower']) for (const curvatureWeight of [0, .4]) {
    const r = distributeCurvatureSurface({ curve: parabola, side, stagnation: 1.4, count: 19, exponent: 1, curvatureWeight });
    const reference = exactArc(2.8) / 2, q = s => exactArc(s) + curvatureWeight * reference * Math.atan(2 * s);
    const end = side === 'upper' ? 0 : 2.8, total = Math.abs(q(end) - q(1.4));
    r.parameters.forEach((s, i) => close(Math.abs(q(s) - q(1.4)), total * i / 18, 3e-10));
    r.spacings.forEach((v, i) => close(v, Math.abs(exactArc(r.parameters[i + 1]) - exactArc(r.parameters[i])), 3e-10));
    assert.equal(r.fractions[0], 0); assert.equal(r.fractions.at(-1), 1);
    assert.equal(r.parameters[0], 1.4); assert.equal(r.parameters.at(-1), end);
  }
});

test('constant curvature yields uniform physical arc; impossible LE control is rejected instead of using negative density', () => {
  const circle = { length: 2 * Math.PI, knots: [0, Math.PI / 2, Math.PI, 3 * Math.PI / 2, 2 * Math.PI], evaluate: s => ({
    point: { x: Math.cos(s), y: Math.sin(s) }, derivative: { x: -Math.sin(s), y: Math.cos(s) }, secondDerivative: { x: -Math.cos(s), y: -Math.sin(s) } }) };
  const input = { curve: circle, side: 'lower', stagnation: Math.PI, count: 17, exponent: .5 };
  const r = distributeCurvatureSurface({ ...input, curvatureWeight: 3 });
  r.fractions.forEach((v, i) => close(v, i / 16, 1e-13));
  assert.throws(() => distributeCurvatureSurface({ ...input, leadingSpacingRatio: .15, artificialLeadingCurvature: false }), /negative|artificial/);
  for (const extra of [{ count: 3 }, { exponent: -1 }, { curvatureWeight: NaN }, { curvatureWeight: 1, leadingSpacingRatio: .2 }])
    assert.throws(() => distributeCurvatureSurface({ ...input, curvatureWeight: 1, ...extra }));
});

test('LE and TE spacing controls match independent physical-arc integration and equidistribute the complete monitor', () => {
  const curve = createContourCurve(naca4('0012', 160)), stagnation = curve.length / 2;
  for (const side of ['upper', 'lower']) {
    const r = distributeCurvatureSurface({ curve, side, stagnation, count: 33, exponent: .5, leadingSpacingRatio: .2, trailingSpacingRatio: .4 });
    const totalArc = integrate(curve, 0, curve.length, speed), reference = totalArc / 2, average = reference / 33;
    const values = r.parameters.slice(1).map((s, i) => integrate(curve, Math.min(s, r.parameters[i]), Math.max(s, r.parameters[i]), speed));
    close(values[0] / average, .2, 2e-9); close(values.at(-1) / average, .4, 2e-9);
    const { curvatureWeight: a, trailingBump: te, totalMetric } = r.diagnostics;
    assert.ok(te.weight > 0 && a > 0); let arc = 0;
    r.parameters.slice(1).forEach((s, i) => {
      const left = arc, right = arc + values[i], h = (right - left) / 256; let bump = 0;
      // Independent Simpson integral of the density bump in physical arc.
      for (let k = 0; k <= 256; k++) {
        const t = (totalArc / 2 - (left + k * h)) / te.width;
        const b = t >= 0 && t <= 1 ? (1 - t) ** 4 * (1 + 4 * t) : 0;
        bump += (k === 0 || k === 256 ? 1 : k % 2 ? 4 : 2) * b;
      }
      bump *= h / 3;
      const curvature = integrate(curve, Math.min(s, r.parameters[i]), Math.max(s, r.parameters[i]), v => {
        const d = v.derivative, dd = v.secondDerivative, q = speed(v);
        return q * Math.sqrt(reference * Math.abs(d.x * dd.y - d.y * dd.x) / q ** 3);
      });
      close(values[i] + a * curvature + te.weight * bump, totalMetric / 32, 3e-7);
      arc = right;
    });
  }
});

test('surface distributions retain rigid-transform and scale invariance on a cambered foil with distinct branches', () => {
  const points = naca4('2412', 160), curve = createContourCurve(points), scale = 2.3;
  const moved = createContourCurve(transform(points, { chord: scale, angle: 27, x: -2, y: .6 }));
  const controls = { count: 33, exponent: .7, leadingSpacingRatio: .2, trailingSpacingRatio: .4 };
  const before = structuredClone(points), results = {};
  for (const side of ['upper', 'lower']) {
    const a = distributeCurvatureSurface({ curve, side, stagnation: .51 * curve.length, ...controls });
    const b = distributeCurvatureSurface({ curve: moved, side, stagnation: .51 * moved.length, ...controls });
    a.fractions.forEach((f, i) => close(f, b.fractions[i], 2e-10));
    a.spacings.forEach((h, i) => close(scale * h, b.spacings[i], 2e-10));
    close(a.diagnostics.curvatureWeight, b.diagnostics.curvatureWeight, 2e-9);
    results[side] = a.fractions;
  }
  assert.notDeepEqual(results.upper, results.lower); assert.deepEqual(points, before);
});

test('artificial LE curvature meets endpoint targets on constant curvature and equidistributes an independently integrated positive monitor', () => {
  const circle = { length: 2 * Math.PI, knots: [0, Math.PI, 2 * Math.PI], evaluate: s => ({
    point: { x: Math.cos(s), y: Math.sin(s) }, derivative: { x: -Math.sin(s), y: Math.cos(s) }, secondDerivative: { x: -Math.cos(s), y: -Math.sin(s) } }) };
  for (const side of ['upper', 'lower']) for (const trailingSpacingRatio of [undefined, .4]) {
    const r = distributeCurvatureSurface({ curve: circle, side, stagnation: Math.PI, count: 17, exponent: .5,
      leadingSpacingRatio: .2, trailingSpacingRatio });
    const d = r.diagnostics;
    assert.equal(d.artificialLeadingCurvature, true);
    close(Math.abs(r.parameters[1] - Math.PI) / (Math.PI / 17), .2, 1e-12);
    if (trailingSpacingRatio) close(Math.abs(r.parameters.at(-1) - r.parameters.at(-2)) / (Math.PI / 17), .4, 1e-12);
    const bump = t => t >= 0 && t <= 1 ? (1 - t) ** 4 * (1 + 4 * t) : 0;
    const metric = s => 1 + d.curvatureWeight * Math.sqrt(Math.PI)
      + d.leadingBump.weight * bump(s / d.leadingBump.width)
      + (d.trailingBump.width ? d.trailingBump.weight * bump((Math.PI - s) / d.trailingBump.width) : 0);
    let last = 0;
    r.parameters.slice(1).forEach(p => {
      const next = Math.abs(p - Math.PI), h = (next - last) / 1024; let sum = 0;
      for (let k = 0; k <= 1024; k++) sum += (k === 0 || k === 1024 ? 1 : k % 2 ? 4 : 2) * metric(last + k * h);
      close(sum * h / 3, d.totalMetric / 16, 2e-9); last = next;
    });
    assert.ok(r.spacings.every(h => h > 0));
  }
});
