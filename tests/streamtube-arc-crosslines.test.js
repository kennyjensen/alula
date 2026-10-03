import test from 'node:test';
import assert from 'node:assert/strict';
import { createArcCrosslineMap } from '../src/geometry/tests/streamtube-arc-crosslines.js';
import { matchPanelPotentialGuides } from '../src/euler/streamtube-potential-guides.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

test('C1 arc interpolation retains stagnation spacing at linear order and inverts its actual cubic', () => {
  // Manufactured stagnation potential phi=s^2, including an interior block.
  const map = createArcCrosslineMap({ anchors: [0, .25, 1], arcAt: Math.sqrt });
  for (const s of [.001, .1, .25, .5, .75, 1]) {
    const chi = map.inverse(s);
    assert.ok(Math.abs(map.value(chi) - s) < 1e-14);
  }
  const h = 1e-7;
  assert.ok(Math.abs(map.inverse(2 * h) / map.inverse(h) - 2) < 1e-6);
  assert.ok(Math.abs(map.inverse(h) / h - .5) < 1e-6);
  assert.deepEqual(map.coordinate, [0, .25, 1]); assert.deepEqual(map.arc, [0, .5, 1]);
  for (let i = 0; i < map.coordinate.length; i++) {
    assert.equal(map.value(map.coordinate[i]), map.arc[i]);
    assert.equal(map.inverse(map.arc[i]), map.coordinate[i]);
  }
  const shifted = createArcCrosslineMap({ anchors: [3, 3.5, 5], arcAt: chi => 7 + 4 * Math.sqrt((chi - 3) / 2) });
  for (const chi of [0, .1, .25, .8, 1]) assert.ok(Math.abs(shifted.value(3 + 2 * chi) - 7 - 4 * map.value(chi)) < 1e-13);
  assert.throws(() => map.value(-.001), /outside/);
  assert.throws(() => map.inverse(-.001), /outside/);
  assert.throws(() => map.inverse(NaN), /outside/);
  assert.throws(() => createArcCrosslineMap({ anchors: [0, 1], arcAt: () => 0 }), /increase/);
});

test('foreign block anchors preserve positive continuous arc derivatives even with strongly unequal secants', () => {
  const anchors = [0, .0001, 1], map = createArcCrosslineMap({ anchors, arcAt: Math.sqrt });
  // The former piecewise-affine map had slopes 100 and 1/1.01: a 101x
  // derivative jump at this foreign anchor despite continuous input arc.
  const anchor = anchors[1], derivative = map.evaluate(anchor).derivative;
  assert.ok(derivative > 0);
  for (const side of [-1, 1]) {
    const nearby = map.evaluate(anchor + side * 1e-12).derivative;
    assert.ok(Math.abs(nearby / derivative - 1) < 1e-5);
  }
  for (let block = 1; block < anchors.length; block++) {
    let previous = -Infinity;
    for (let i = 0; i <= 32; i++) {
      const chi = i === 32 ? anchors[block] : anchors[block - 1] + (anchors[block] - anchors[block - 1]) * i / 32;
      const { value, derivative } = map.evaluate(chi);
      assert.ok(derivative > 0); assert.ok(value > previous); previous = value;
      assert.ok(Math.abs(map.inverse(value) - chi) < 3e-14);
    }
  }
  assert.deepEqual(anchors, [0, .0001, 1]);
  assert.equal(map.value(anchor), .01); assert.equal(map.inverse(.01), anchor);
});

test('farfield LE/TE anchors stay fixed when a different surface is refined', () => {
  const curves = [createContourCurve(naca4('0012', 40)),
    createContourCurve(transform(naca4('0012', 40), { chord: .3, x: .9, y: -.1 }))];
  const profiles = curves.map((curve, b) => {
    const stag = curve.length / 2, phiStag = b ? .7 : 0;
    return { curve, stag, phiStag, phase: s => phiStag + (s - stag) ** 2,
      upstream: [{ x: -4, y: b, potential: -4 }], wake: [{ x: 6, y: b, potential: 5 }] };
  });
  const sample = (path, potential) => ({ x: potential, y: path[0].y, potential });
  const build = count => matchPanelPotentialGuides({ input: { bodies: [{}, {}], gridSpacing: { inlet: { intervals: 12 }, outlet: { intervals: 12 } } }, profiles,
    surfaceFractions: [16, count].map(n => Object.fromEntries(['upper', 'lower'].map(side => [side,
      Array.from({ length: n + 1 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / n)))]))),
    outerPaths: [-2, 2].map(y => [{ x: -4, y, potential: -4 }, { x: 6, y, potential: 6 }]),
    sample, sampleCut: sample, growth: .25, maxIntervals: 1000, outerSpread: .85, resolveTurns: true, surfaceMetric: 'arc' });
  const a = build(16), b = build(32);
  assert.notEqual(a.diagnostics.x.length, b.diagnostics.x.length);
  assert.deepEqual(a.diagnostics.anchors, b.diagnostics.anchors);
  for (const anchor of a.diagnostics.anchors) {
    const ia = a.diagnostics.x.indexOf(anchor), ib = b.diagnostics.x.indexOf(anchor);
    for (let side = 0; side < 2; side++) assert.deepEqual(a.outer[side][ia], b.outer[side][ib]);
  }
  assert.equal(a.diagnostics.outerSpreading, 'C1 monotone through LE/TE anchors');
});
