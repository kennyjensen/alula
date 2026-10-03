import test from 'node:test';
import assert from 'node:assert/strict';
import { createContourArc } from '../src/geometry/contour-arc.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const close = (a, b, tolerance) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
test('physical arc quadrature recovers an analytic parabola and its speed', () => {
  const curve = { knots: [0, .15, .55, 1.4], length: 1.4,
    evaluate: s => ({ point: { x: s, y: s * s }, derivative: { x: 1, y: 2 * s } }) };
  const arc = createContourArc(curve), exact = s => .5 * s * Math.sqrt(1 + 4 * s * s) + .25 * Math.asinh(2 * s);
  for (const s of [0, .08, .15, .4, .55, .81, 1.4]) {
    close(arc.at(s), exact(s), 5e-10); close(arc.speed(s), Math.sqrt(1 + 4 * s * s), 1e-14);
  }
  close(arc.length, exact(1.4), 5e-10);
  assert.throws(() => arc.at(-.1), /outside/); assert.throws(() => arc.at(1.5), /outside/);
});

test('foil physical arc agrees with independent composite Simpson integration and transforms with geometry', () => {
  const points = naca4('2412', 80), curve = createContourCurve(points), arc = createContourArc(curve);
  const mapped = createContourCurve(transform(points, { chord: 2.3, angle: 29, x: -2, y: .7 })), moved = createContourArc(mapped);
  const simpson = end => {
    let total = 0;
    for (let i = 0; i < curve.knots.length - 1 && curve.knots[i] < end; i++) {
      const a = curve.knots[i], b = Math.min(end, curve.knots[i + 1]), h = (b - a) / 64;
      let sum = 0;
      for (let k = 0; k <= 64; k++) {
        const d = curve.evaluate(k === 64 ? b : a + k * h).derivative;
        sum += (k === 0 || k === 64 ? 1 : k % 2 ? 4 : 2) * Math.hypot(d.x, d.y);
      }
      total += sum * h / 3;
    }
    return total;
  };
  for (const f of [.13, .49, .51, .73, 1]) {
    const s = f * curve.length;
    close(arc.at(s), simpson(s), 5e-11);
    close(moved.at(f * mapped.length), 2.3 * arc.at(s), 2e-12);
    if (f < 1) {
      const h = 1e-6;
      close((arc.at(s + h) - arc.at(s - h)) / (2 * h), arc.speed(s), 5e-8);
    }
  }
  assert.ok(arc.length > curve.length, 'physical spline length exceeds the inscribed polygon length');
});
