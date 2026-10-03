// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createContourTopology, createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const fixture = JSON.parse(fs.readFileSync('third_party/airfoils/nlr7301/geometry.json'));
const close = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const pointClose = (a, b, tolerance) => ['x', 'y'].forEach(k => close(a[k], b[k], tolerance));
const rotate = (p, angle) => transform([p], { angle })[0];
const metadata = lowerIndex => ({ trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex } });

test('NLR explicit TE indices partition every original surface and base panel exactly once', () => {
  const elements = getBenchmarkAirfoil('nlr7301').elements;
  for (const [k, element] of elements.entries()) {
    const source = fixture.elements[k], p = element.points;
    const topology = createContourTopology(p, element), lowerIndex = k === 0 ? 424 : 216;
    assert.deepEqual(element.trailingEdge, { kind: 'finite-base', upperIndex: 0, lowerIndex });
    assert.deepEqual(topology.points, source.points.map(([x, y]) => ({ x, y })));
    assert.deepEqual(topology.surface.points, p.slice(0, lowerIndex + 1));
    assert.deepEqual(topology.base.points, p.slice(lowerIndex));
    assert.deepEqual(topology.surface.points.concat(topology.base.points.slice(1)), p);
    assert.deepEqual(topology.surface.indices, Array.from({ length: lowerIndex + 1 }, (_, i) => i));
    assert.deepEqual(topology.base.indices, Array.from({ length: 33 }, (_, i) => lowerIndex + i));
    assert.equal(topology.base.panels.length, 32);
    assert.equal(topology.surface.panels.length, lowerIndex);
    const panels = [...topology.surface.panels, ...topology.base.panels];
    assert.deepEqual(panels.map(q => q.sourcePanelIndex), Array.from({ length: p.length - 1 }, (_, i) => i));
    for (const q of panels) {
      assert.deepEqual(q.start, p[q.sourcePanelIndex]); assert.deepEqual(q.end, p[q.sourcePanelIndex + 1]);
      // hypot, division and two squared components each round; the observed
      // 3-ulp norm error is arithmetic, not a geometric direction change.
      close(q.tangent.x ** 2 + q.tangent.y ** 2, 1, 8 * Number.EPSILON);
      close(q.outwardNormal.x * q.tangent.x + q.outwardNormal.y * q.tangent.y, 0, 5e-16);
      assert.deepEqual(q.outwardNormal, { x: q.tangent.y, y: -q.tangent.x });
    }
    assert.deepEqual(topology.trailingEdge.upper.point, { x: source.upperTrailingEdge[0], y: source.upperTrailingEdge[1] });
    assert.deepEqual(topology.trailingEdge.lower.point, { x: source.lowerTrailingEdge[0], y: source.lowerTrailingEdge[1] });
    close(topology.trailingEdge.gapLength, source.trailingEdgeBaseLength, 1e-16);
    assert.ok(topology.trailingEdge.normalGap > 0);
    assert.equal(topology.units.lengths, 'input length units');
  }
});

test('finite-base C2 fit excludes the base and keeps distinct stationary branch endpoints', () => {
  for (const element of getBenchmarkAirfoil('nlr7301').elements) {
    const topology = createContourTopology(element.points, element);
    const curve = createSurfaceContourCurve(element.points, element), stagnation = .48 * curve.length;
    assert.equal(curve.knots.length, topology.surface.points.length);
    close(curve.length, topology.surface.length, 3e-15);
    assert.ok(curve.length < topology.surface.length + topology.base.length);
    curve.knots.forEach((s, i) => pointClose(curve.evaluate(s).point, topology.surface.points[i], 1e-15));
    for (const side of ['upper', 'lower']) {
      const end = curve.branch(side, 1, stagnation), movedEnd = curve.branch(side, 1, stagnation + .01 * curve.length);
      assert.deepEqual(end.point, topology.trailingEdge[side].point);
      assert.deepEqual(movedEnd.point, end.point);
      pointClose(end.stagnationDerivative, { x: 0, y: 0 }, 0);
      const f = .73, h = 1e-6 * curve.length, value = curve.branch(side, f, stagnation);
      const plus = curve.branch(side, f, stagnation + h), minus = curve.branch(side, f, stagnation - h);
      pointClose(value.stagnationDerivative,
        { x: (plus.point.x - minus.point.x) / (2 * h), y: (plus.point.y - minus.point.y) / (2 * h) }, 2e-9);
    }
    assert.notDeepEqual(curve.branch('upper', 1, stagnation).point, curve.branch('lower', 1, stagnation).point);
    assert.deepEqual(curve.branch('upper', 0, stagnation).point, curve.branch('lower', 0, stagnation).point);
  }
});

test('explicit finite-base topology and surface interpolation respect rotation, translation and scale', () => {
  for (const element of getBenchmarkAirfoil('nlr7301').elements) {
    const a = createContourTopology(element.points, element), ca = createSurfaceContourCurve(element.points, element);
    for (const map of [{ chord: 2.7, angle: 37, x: -1.2, y: .9 }, { chord: .25, angle: -63, x: .4, y: -.6 }]) {
      const p = transform(element.points, map), b = createContourTopology(p, element), cb = createSurfaceContourCurve(p, element);
      assert.deepEqual(b.points, p); assert.deepEqual(b.surface.indices, a.surface.indices); assert.deepEqual(b.base.indices, a.base.indices);
      close(b.base.length, map.chord * a.base.length, 2e-14);
      close(b.trailingEdge.gapLength, map.chord * a.trailingEdge.gapLength, 2e-14);
      close(b.trailingEdge.normalGap, map.chord * a.trailingEdge.normalGap, 2e-14);
      close(b.trailingEdge.tangentialOffset, map.chord * a.trailingEdge.tangentialOffset, 2e-14);
      pointClose(b.trailingEdge.center, transform([a.trailingEdge.center], map)[0], 2e-15);
      pointClose(b.trailingEdge.wakeTangent, rotate(a.trailingEdge.wakeTangent, map.angle), 2e-11);
      pointClose(b.trailingEdge.wakeNormal, rotate(a.trailingEdge.wakeNormal, map.angle), 2e-11);
      for (const side of ['upper', 'lower']) {
        pointClose(b.trailingEdge[side].outwardNormal, rotate(a.trailingEdge[side].outwardNormal, map.angle), 2e-11);
        pointClose(b.trailingEdge[side].downstreamTangent, rotate(a.trailingEdge[side].downstreamTangent, map.angle), 2e-11);
      }
      for (const f of [0, .13, .42, .81, 1]) pointClose(cb.evaluate(f * cb.length).point,
        transform([ca.evaluate(f * ca.length).point], map)[0], 3e-13);
    }
  }
});

test('a curved finite base and a cyclic source start retain their full polylines', () => {
  const points = [[1, .1], [.8, .12], [.3, .15], [0, 0], [.3, -.15], [.8, -.12],
    [1, -.1], [1.03, -.05], [1.05, 0], [1.03, .05], [1, .1]].map(([x, y]) => ({ x, y }));
  const a = createContourTopology(points, metadata(6));
  assert.deepEqual(a.base.points, points.slice(6)); assert.equal(a.base.panels.length, 4);
  assert.ok(a.base.length > a.trailingEdge.gapLength);
  assert.ok(a.base.points.some(p => p.x > a.trailingEdge.center.x));
  const unique = points.slice(0, -1), shifted = unique.slice(3).concat(unique.slice(0, 3)); shifted.push({ ...shifted[0] });
  const b = createContourTopology(shifted, { trailingEdge: { kind: 'finite-base', upperIndex: 7, lowerIndex: 3 } });
  assert.deepEqual(b.points, shifted); assert.deepEqual(b.surface.points, a.surface.points); assert.deepEqual(b.base.points, a.base.points);
  assert.equal(new Set([...b.surface.panels, ...b.base.panels].map(q => q.sourcePanelIndex)).size, points.length - 1);
});

test('finite-base indices must be explicit and unambiguous; source coordinates are never snapped or reversed', () => {
  const e = getBenchmarkAirfoil('nlr7301').elements[0];
  assert.throws(() => createContourTopology(e.points), /sharp trailing edge/);
  assert.throws(() => createContourTopology(e.points, { trailingEdge: { kind: 'blunt' } }), /explicit/);
  for (const lowerIndex of [undefined, 0, 1, 456, 1.5]) assert.throws(() => createContourTopology(e.points, metadata(lowerIndex)), /indices|four points/);
  const reversed = e.points.slice().reverse();
  assert.throws(() => createContourTopology(reversed, e), /counterclockwise/);
  const open = structuredClone(e.points); open.at(-1).x += 1e-13;
  assert.throws(() => createContourTopology(open, e), /repeated closing point/);
  const openSurface = e.points.slice(0, e.trailingEdge.lowerIndex + 1);
  assert.throws(() => createContourCurve(openSurface), /sharp trailing edge/);
  assert.throws(() => createContourCurve(openSurface, { allowOpenEndpoints: 1 }), /Invalid/);
});

test('sharp-TE default curve results are exactly equal to the archived implementation', async () => {
  const archive = 'docs/nlr-finite-base/before/geometry/contour-curve.js.txt';
  const manifest = JSON.parse(fs.readFileSync('docs/nlr-finite-base/geometry-source-before.json'));
  const expected = manifest.files.find(f => f.path === 'src/geometry/contour-curve.js');
  assert.equal(createHash('sha256').update(fs.readFileSync(archive)).digest('hex'), expected.sha256);
  const old = await import('data:text/javascript;base64,' + fs.readFileSync(archive).toString('base64'));
  for (const points of [naca4('2412', 40), getBenchmarkAirfoil('rae2822').elements[0].points,
    getBenchmarkAirfoil('rae2822-mses').elements[0].points]) {
    const a = old.createContourCurve(points), b = createContourCurve(points), c = createSurfaceContourCurve(points);
    const t = createContourTopology(points);
    assert.equal(t.kind, 'sharp'); assert.deepEqual(t.points, points); assert.deepEqual(t.surface.points, points);
    assert.equal(t.base.panels.length, 0); assert.equal(t.base.length, 0); assert.equal(t.trailingEdge.gapLength, 0);
    assert.equal(a.length, b.length); assert.equal(a.length, c.length); assert.deepEqual(a.knots, b.knots); assert.deepEqual(a.knots, c.knots);
    for (const f of [0, .03, .37, .62, .95, 1]) {
      assert.deepEqual(a.evaluate(f * a.length), b.evaluate(f * b.length));
      assert.deepEqual(a.evaluate(f * a.length), c.evaluate(f * c.length));
      for (const side of ['upper', 'lower']) assert.deepEqual(a.branch(side, f, .48 * a.length), c.branch(side, f, .48 * c.length));
    }
  }
});

test('topology data and preset metadata are detached from their callers', () => {
  const e = getBenchmarkAirfoil('nlr7301').elements[0], original = structuredClone(e);
  const a = createContourTopology(e.points, e);
  a.points[0].x = 999; a.surface.points[0].y = 777; a.base.panels[0].start.x = -111;
  a.trailingEdge.upper.point.x = 888;
  assert.deepEqual(e, original);
  e.trailingEdge.lowerIndex = 1;
  assert.equal(getBenchmarkAirfoil('nlr7301').elements[0].trailingEdge.lowerIndex, 424);
});
