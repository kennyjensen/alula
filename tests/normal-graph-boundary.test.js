// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createNormalGraphBoundary } from '../src/geometry/normal-graph-boundary.js';

const close = (a, b, tolerance = 3e-13) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const graph = (xs, f, df) => createNormalGraphBoundary({ points: xs.map(x => ({ x, y: f(x) })), slopes: xs.map(df) });

test('normal graph exactly preserves Hermite knots, endpoint slopes, and the C1 interpolant', () => {
  const points = [{ x: -1, y: .1 }, { x: .2, y: -.03 }, { x: 2, y: .2 }], slopes = [-.05, .02, .08];
  const curve = createNormalGraphBoundary({ points, slopes });
  points.forEach((p, i) => {
    assert.deepEqual(curve.evaluate(p.x).point, p);
    assert.equal(curve.evaluate(p.x).derivative.y, slopes[i]);
  });
  const h = 1e-6;
  for (const sign of [-1, 1]) {
    const coarse = Math.abs(curve.evaluate(.2 + sign * h).derivative.y - slopes[1]);
    const fine = Math.abs(curve.evaluate(.2 + sign * h / 2).derivative.y - slopes[1]);
    assert.ok(coarse < 1e-6);
    close(coarse / fine, 2, 5e-6);
  }
  for (const x of [-.5, .6, 1.5]) {
    const q = curve.evaluate(x), a = x < .2 ? 0 : 1, width = points[a + 1].x - points[a].x, t = (x - points[a].x) / width;
    const expected = (2 * t ** 3 - 3 * t ** 2 + 1) * points[a].y + (t ** 3 - 2 * t ** 2 + t) * width * slopes[a]
      + (-2 * t ** 3 + 3 * t ** 2) * points[a + 1].y + (t ** 3 - t ** 2) * width * slopes[a + 1];
    close(q.point.y, expected);
    close((curve.evaluate(x + h).point.y - curve.evaluate(x - h).point.y) / (2 * h), q.derivative.y, 3e-11);
    close((curve.evaluate(x + h).derivative.y - curve.evaluate(x - h).derivative.y) / (2 * h), q.secondDerivative.y, 3e-11);
  }
});

test('normal projection on horizontal and sloped lines matches the independent orthogonal foot', () => {
  const horizontal = graph([-2, 0, 3], () => .5, () => 0);
  for (const x of [-2, -.3, 1.7, 3]) assert.deepEqual(horizontal.projectNormal({ x, y: -4 }), { x, y: .5 });
  for (const slope of [-2, .5, 3]) {
    const curve = graph([-3, 0, 4], x => 2 + slope * x, () => slope);
    for (const p of [{ x: .2, y: 3 }, { x: -1, y: 1 }]) {
      const x = (p.x + slope * (p.y - 2)) / (1 + slope * slope), q = curve.projectNormal(p);
      close(q.x, x); close(q.y, 2 + slope * x);
      close((q.x - p.x) + (q.y - p.y) * slope, 0, 2e-12);
    }
  }
  const angle = .4, rotate = p => ({ x: p.x * Math.cos(angle) - p.y * Math.sin(angle), y: p.x * Math.sin(angle) + p.y * Math.cos(angle) });
  const line = createNormalGraphBoundary({ points: [-2, 0, 3].map(x => rotate({ x, y: .2 })), slopes: Array(3).fill(Math.tan(angle)) });
  const foot = line.projectNormal(rotate({ x: 1.1, y: .7 })), exact = rotate({ x: 1.1, y: .2 });
  close(foot.x, exact.x); close(foot.y, exact.y);
});

test('curved normal feet solve the independent tangent equation on both sides and at C1 knots', () => {
  const f = x => .02 * x ** 3 - .05 * x, df = x => .06 * x * x - .05;
  const curve = graph([-2, -1, 0, .7, 2], f, df);
  for (const x of [-1.3, -1, 0, .7, 1.4]) for (const distance of [-.4, .3]) {
    const p = { x: x - distance * df(x), y: f(x) + distance }, q = curve.projectNormal(p);
    close(q.x, x); close(q.y, f(x));
    close((q.x - p.x) + (q.y - p.y) * df(q.x), 0, 2e-13);
    assert.ok(q.x > -2 && q.x < 2);
  }
});

test('normal graph copies immutable descriptors and respects geometric units and translation', () => {
  const input = { points: [{ x: -1, y: .1 }, { x: 0, y: 0 }, { x: 1, y: .1 }], slopes: [-.2, 0, .2] };
  const saved = structuredClone(input), curve = createNormalGraphBoundary(input);
  input.points[0].y = 100; input.slopes[1] = 100;
  assert.deepEqual(curve.descriptor, { kind: 'piecewise-cubic-hermite-graph', ...saved });
  assert.deepEqual(createNormalGraphBoundary(curve.descriptor).evaluate(.2), curve.evaluate(.2));
  assert.ok(Object.isFrozen(curve) && Object.isFrozen(curve.descriptor) && Object.isFrozen(curve.points)
    && curve.points.every(Object.isFrozen) && Object.isFrozen(curve.slopes));
  assert.throws(() => { curve.points[0].x = 2; }, TypeError);
  const p = { x: .28, y: .37 }, original = curve.projectNormal(p);
  for (const scale of [.1, 7]) {
    const transform = p => ({ x: 3 + scale * p.x, y: -2 + scale * p.y });
    const scaled = createNormalGraphBoundary({ points: saved.points.map(transform), slopes: saved.slopes });
    const q = scaled.projectNormal(transform(p)), expected = transform(original);
    close(q.x, expected.x); close(q.y, expected.y);
    close(scaled.evaluate(3 + scale * .2).derivative.y, curve.evaluate(.2).derivative.y);
    close(scaled.evaluate(3 + scale * .2).secondDerivative.y, curve.evaluate(.2).secondDerivative.y / scale, 5e-12);
  }
});

test('normal graph rejects invalid data, unproved uniqueness, and missing root brackets without snapping', () => {
  const input = { points: [{ x: -1, y: 0 }, { x: 1, y: 0 }], slopes: [0, 0] };
  for (const invalid of [{}, { ...input, points: [input.points[0]] }, { ...input, slopes: [0] },
    { ...input, slopes: [0, NaN] }, { ...input, slopes: [, 0] }, { ...input, points: [input.points[1], input.points[0]] },
    { ...input, points: [input.points[0], input.points[0]] }, { ...input, points: [input.points[0], { x: 1, y: Infinity }] }])
    assert.throws(() => createNormalGraphBoundary(invalid), /finite|matching|increasing/);
  const curve = createNormalGraphBoundary(input);
  for (const x of [-1.01, 1.01, NaN, Infinity]) assert.throws(() => curve.evaluate(x), /outside/);
  assert.throws(() => curve.projectNormal({ x: 1 + Number.EPSILON, y: .1 }), /not bracketed/);
  assert.throws(() => curve.projectNormal({ x: -2, y: .1 }), /not bracketed/);
  assert.throws(() => curve.projectNormal({ x: 0, y: NaN }), /finite/);
  const parabola = graph([-1, 0, 1], x => x * x, x => 2 * x);
  // F=2*x^3-x has three roots. The sufficient derivative test must refuse.
  assert.throws(() => parabola.projectNormal({ x: 0, y: 1 }), /Cannot establish a unique/);
  // This root actually is unique (F=2*x^3+3*x), but the deliberately
  // conservative absolute-curvature bound cannot establish that fact.
  assert.throws(() => parabola.projectNormal({ x: 0, y: -1 }), /Cannot establish a unique/);
});
