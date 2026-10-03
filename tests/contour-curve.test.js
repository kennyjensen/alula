import test from 'node:test';
import assert from 'node:assert/strict';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const close = (a, b, tolerance = 2e-9) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);

test('body contour spline interpolates the supplied geometry and has independently checked continuous derivatives', () => {
  const points = naca4('0012', 80), curve = createContourCurve(points);
  curve.knots.forEach((s, i) => { const p = curve.evaluate(s).point; close(p.x, points[i].x, 1e-14); close(p.y, points[i].y, 1e-14); });
  for (const fraction of [.11, .36, .51, .71, .92]) {
    const s = fraction * curve.length, h = 1e-6, value = curve.evaluate(s), plus = curve.evaluate(s + h), minus = curve.evaluate(s - h);
    for (const key of ['x', 'y']) {
      close(value.derivative[key], (plus.point[key] - minus.point[key]) / (2 * h));
      close(value.secondDerivative[key], (plus.derivative[key] - minus.derivative[key]) / (2 * h), 1e-7);
    }
  }
  for (const s of curve.knots.slice(1, -1)) {
    const a = curve.evaluate(s - 1e-10), b = curve.evaluate(s + 1e-10);
    for (const key of ['x', 'y']) { close(a.derivative[key], b.derivative[key], 1e-7); close(a.secondDerivative[key], b.secondDerivative[key], 1e-5); }
  }
  assert.throws(() => curve.evaluate(-1), /outside/);
});

test('stagnation motion redistributes both surface branches and preserves the sharp trailing edge', () => {
  const curve = createContourCurve(naca4('0012', 100)), stagnation = curve.length * .48;
  for (const side of ['upper', 'lower']) for (const fraction of [0, .1, .6, 1]) {
    const h = 1e-6, r = curve.branch(side, fraction, stagnation);
    const plus = curve.branch(side, fraction, stagnation + h), minus = curve.branch(side, fraction, stagnation - h);
    for (const key of ['x', 'y']) close(r.stagnationDerivative[key], (plus.point[key] - minus.point[key]) / (2 * h), 1e-8);
    if (fraction === 1) { close(r.point.x, 1, 1e-14); close(r.point.y, 0, 1e-14); close(r.stagnationDerivative.x, 0, 1e-14); }
  }
  assert.deepEqual(curve.branch('upper', 0, stagnation).point, curve.branch('lower', 0, stagnation).point);
  const moved = curve.branch('upper', .1, stagnation + .01).point, original = curve.branch('upper', .1, stagnation).point;
  assert.ok(Math.hypot(moved.x - original.x, moved.y - original.y) > .001);
  assert.throws(() => curve.branch('upper', .5, 0), /Invalid/);
});

test('parametric body interpolation is invariant under rigid transforms and approaches an independent circle under refinement', () => {
  const points = naca4('2412', 80), map = { chord: 2.1, angle: 31, x: 3, y: -2 };
  const original = createContourCurve(points), moved = createContourCurve(transform(points, map));
  for (const fraction of [.03, .31, .52, .78, .95]) {
    const expected = transform([original.evaluate(fraction * original.length).point], map)[0], actual = moved.evaluate(fraction * moved.length).point;
    close(expected.x, actual.x, 1e-12); close(expected.y, actual.y, 1e-12);
  }
  const errors = [];
  for (const n of [32, 64, 128]) {
    const points = Array.from({ length: n + 1 }, (_, i) => ({ x: Math.cos(2 * Math.PI * i / n), y: Math.sin(2 * Math.PI * i / n) }));
    const curve = createContourCurve(points); let error = 0;
    // Exclude the deliberately independent TE ends/natural endpoint curvature.
    for (let i = n / 4; i < 3 * n / 4; i++) { const p = curve.evaluate((i + .5) * curve.length / n).point; error = Math.max(error, Math.abs(Math.hypot(p.x, p.y) - 1)); }
    errors.push(error);
  }
  assert.ok(errors[1] < .08 * errors[0], String(errors)); assert.ok(errors[2] < .08 * errors[1], String(errors));
});
