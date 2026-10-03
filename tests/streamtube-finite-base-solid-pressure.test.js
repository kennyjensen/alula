// SPDX-License-Identifier: GPL-2.0-or-later
// Manufactured pressure fields only: no Euler/BL residual or solve imports.
import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeSolidPressureForces } from '../src/euler/streamtube-forces.js';

const keys = ['cx', 'cy', 'cl', 'cm', 'pressureIntegralDrag'];
const close = (a, b) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) < 3e-12 * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
const contour = [
  { x: 2, y: .2 }, { x: 1, y: .5 }, { x: 0, y: 0 }, { x: 1, y: -.5 }, { x: 2, y: -.2 },
  { x: 2.1, y: -.1 }, { x: 2.16, y: 0 }, { x: 2.1, y: .1 }, { x: 2, y: .2 },
];
const pressureAt = (body, side, i) => i === 1 ? 3.8 : side === 'upper' ? (i === 2 ? 4.4 : 5.2) : (i === 2 ? 4.9 : 3.1);

function fixture(pressure = () => 3.5, count = 1) {
  const nx = 4, tubes = Array(count + 1).fill(1);
  const bodies = Array.from({ length: count }, (_, b) => ({ leadingIndex: 1, trailingIndex: 3, element: 10 - b,
    points: contour.map(p => ({ x: p.x, y: p.y + 3 * b })), trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 4 } }));
  const surface = (b, side, i) => i === 0 || i === nx ? { x: i - 1, y: 3 * b } : bodies[b].points[
    i === 1 ? 2 : side === 'upper' ? (i === 2 ? 1 : 0) : (i === 2 ? 3 : 4)];
  const nodes = tubes.map((_, g) => Array.from({ length: nx + 1 }, (_, i) => [
    g === 0 ? { x: i - 1, y: -2 } : { ...surface(g - 1, 'upper', i) },
    g === count ? { x: i - 1, y: 3 * count } : { ...surface(g, 'lower', i) },
  ]));
  const cells = Array.from({ length: nx - 1 }, () => tubes.map(() => [{ interfacePressure: { lower: 3.5, upper: 3.5 } }]));
  bodies.forEach((body, b) => {
    for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
      cells[i - 1][b][0].interfacePressure.upper = pressure(b, 'lower', i);
      cells[i - 1][b + 1][0].interfacePressure.lower = pressure(b, 'upper', i);
    }
  });
  return { flow: { nodes, cells }, layout: { nx, tubes, bodies }, conditions: { flowModel: 'compressible', pInf: 2, alpha: 0 } };
}

// Independent two-point Gauss integration of traction and moment. The
// integrand is quadratic at most, so this is exact on each straight edge.
// Do not import the production polygon quadrature or topology helper.
function integrate(points, pressures, { referenceChord = 1, momentReference = { x: referenceChord / 4, y: 0 }, alpha = 0, pInf = 2 } = {}) {
  let cx = 0, cy = 0, cm = 0;
  const ts = [.5 - .5 / Math.sqrt(3), .5 + .5 / Math.sqrt(3)];
  for (let k = 1; k < points.length; k++) {
    const a = points[k - 1], b = points[k], dx = (b.x - a.x) / referenceChord, dy = (b.y - a.y) / referenceChord;
    for (const t of ts) {
      const cp = 2 * ((1 - t) * pressures[k - 1] + t * pressures[k] - pInf);
      const x = ((1 - t) * a.x + t * b.x - momentReference.x) / referenceChord;
      const y = ((1 - t) * a.y + t * b.y - momentReference.y) / referenceChord;
      cx -= .5 * cp * dy; cy += .5 * cp * dx; cm -= .5 * cp * (x * dx + y * dy);
    }
  }
  const angle = alpha * Math.PI / 180;
  return { cx, cy, cm, cl: cy * Math.cos(angle) - cx * Math.sin(angle), pressureIntegralDrag: cx * Math.cos(angle) + cy * Math.sin(angle) };
}

function explicitPressures(input, b = 0) {
  const body = input.layout.bodies[b], cells = input.flow.cells;
  const p = side => [1, 2, 3].map(i => side === 'upper' ? cells[i - 1][b + 1][0].interfacePressure.lower : cells[i - 1][b][0].interfacePressure.upper);
  const upper = p('upper'), lower = p('lower'), result = [upper[2], upper[1], .5 * (upper[0] + lower[0]), lower[1], lower[2]];
  const lengths = body.points.slice(5).map((p, i) => distance(body.points[i + 4], p));
  const total = lengths.reduce((a, b) => a + b, 0); let arc = 0;
  for (const length of lengths) { arc += length; result.push(lower[2] + (upper[2] - lower[2]) * arc / total); }
  return result;
}

test('finite-base uniform pressure cancels on both entire closed solid contours', () => {
  const input = fixture(() => 8.5, 2), before = structuredClone(input);
  const result = streamtubeSolidPressureForces({ ...input, referenceChord: 1.7, momentReference: { x: -.3, y: .8 } });
  for (const key of keys) { close(result[key], 0); result.perBody.forEach(body => close(body[key], 0)); }
  result.perBody.forEach(body => { close(body.enclosedSolidArea, 1.236); assert.equal(body.basePressureModel.panelCount, 4); });
  assert.deepEqual(result.perBody.map(b => b.element), [10, 9]);
  assert.deepEqual(input, before);
  assert.equal(result.physicalAcceptance, false); assert.equal(result.includesSkinFriction, false);
  assert.equal(result.includesExitDefects, false); assert.equal(Object.hasOwn(result, 'cd'), false);
});

test('unequal TE pressures use retained bowed-base arc interpolation and independent exact quadrature', () => {
  const input = fixture(pressureAt), before = structuredClone(input), options = { referenceChord: .7, momentReference: { x: -.4, y: .6 } };
  input.conditions.alpha = 17;
  const actual = streamtubeSolidPressureForces({ ...input, ...options });
  const pressures = explicitPressures(input), points = input.layout.bodies[0].points;
  const expected = integrate(points, pressures, { ...options, ...input.conditions });
  for (const key of keys) close(actual[key], expected[key]);
  assert.deepEqual(actual.perBody[0].basePressureModel, { kind: 'TE-pressure interpolation on retained base', panelCount: 4, upper: 5.2, lower: 3.1 });
  const collapsed = integrate([...points.slice(0, 5), points[0]], [...pressures.slice(0, 5), pressures[0]], { ...options, ...input.conditions });
  assert.ok(Math.abs(actual.cy - collapsed.cy) > 1e-3, 'A single closing edge must not substitute for the bowed base.');
  const meanBase = pressures.map((p, i) => i >= 4 ? .5 * (pressures[0] + pressures[4]) : p);
  const meanBaseOnly = integrate(points.slice(4), meanBase.slice(4), { ...options, ...input.conditions });
  const arcBaseOnly = integrate(points.slice(4), pressures.slice(4), { ...options, ...input.conditions });
  assert.ok(Math.abs(meanBaseOnly.cy - arcBaseOnly.cy) > 1e-3, 'The unequal-TE model is arc interpolation, not constant mean pressure.');
  assert.deepEqual(input.flow, before.flow); assert.deepEqual(input.layout, before.layout);
});

test('symmetric finite-base pressures give zero lift and moment; reflection reverses asymmetric loads', () => {
  const symmetric = fixture((b, side, i) => [0, 3.8, 4.4, 3.1][i]);
  const zero = streamtubeSolidPressureForces(symmetric); close(zero.cl, 0); close(zero.cm, 0);
  const original = streamtubeSolidPressureForces(fixture(pressureAt));
  const reflected = streamtubeSolidPressureForces(fixture((b, side, i) => pressureAt(b, side === 'upper' ? 'lower' : 'upper', i)));
  close(reflected.cx, original.cx); close(reflected.pressureIntegralDrag, original.pressureIntegralDrag);
  close(reflected.cy, -original.cy); close(reflected.cl, -original.cl); close(reflected.cm, -original.cm);
});

test('fresh inviscid base wake and total viscous displacement both integrate the preserved solid', () => {
  const input = fixture(pressureAt), expected = streamtubeSolidPressureForces(input);
  for (const wallThickness of [0, .08]) {
    const displaced = structuredClone(input);
    displaced.layout.displacedBoundaries = true;
    displaced.flow.displacement = { surfaces: [{ upper: Array(3).fill(wallThickness), lower: Array(3).fill(2 * wallThickness) }], wakes: [[.4 + 3 * wallThickness]] };
    displaced.flow.undisplacedNodes = structuredClone(input.flow.nodes);
    for (let i = 1; i <= 3; i++) {
      displaced.flow.nodes[0][i][1].y -= 2 * wallThickness;
      displaced.flow.nodes[1][i][0].y += wallThickness;
    }
    displaced.flow.nodes[0][4][1].y -= .2 + 2 * wallThickness;
    displaced.flow.nodes[1][4][0].y += .2 + wallThickness;
    displaced.flow.boundaryLayer = { ue: NaN, deltaStar: wallThickness };
    const before = structuredClone(displaced), actual = streamtubeSolidPressureForces(displaced);
    for (const key of keys) assert.equal(actual[key], expected[key]);
    assert.deepEqual(actual.perBody, expected.perBody);
    assert.equal(actual.geometry, 'undisplaced-solid-grid'); assert.deepEqual(displaced, before);
    delete displaced.flow.undisplacedNodes;
    assert.throws(() => streamtubeSolidPressureForces(displaced), /undisplaced solid grid/);
  }
});

test('finite-base loads transform with geometry and references, including every retained base vertex', () => {
  const input = fixture(pressureAt), baseReference = { x: .3, y: -.2 };
  const baseline = streamtubeSolidPressureForces({ ...input, momentReference: baseReference });
  const scale = 2.3, angle = 31 * Math.PI / 180, c = Math.cos(angle), s = Math.sin(angle);
  const move = p => ({ x: 1.2 + scale * (c * p.x - s * p.y), y: -.8 + scale * (s * p.x + c * p.y) });
  const moved = structuredClone(input);
  moved.flow.nodes = moved.flow.nodes.map(group => group.map(row => row.map(move)));
  moved.layout.bodies.forEach(body => { body.points = body.points.map(move); }); moved.conditions.alpha = 31;
  const actual = streamtubeSolidPressureForces({ ...moved, referenceChord: scale, momentReference: move(baseReference) });
  close(actual.cx, c * baseline.cx - s * baseline.cy); close(actual.cy, s * baseline.cx + c * baseline.cy);
  close(actual.cl, baseline.cl); close(actual.cm, baseline.cm); close(actual.pressureIntegralDrag, baseline.pressureIntegralDrag);
  close(actual.perBody[0].enclosedSolidArea, 1.236 * scale ** 2);
});

test('finite-base force guards reject shifted solid corners or omitted finite-base topology', () => {
  const input = fixture(pressureAt), shifted = structuredClone(input);
  shifted.flow.undisplacedNodes = structuredClone(input.flow.nodes);
  shifted.flow.undisplacedNodes[1][3][0].x += 1e-8;
  assert.throws(() => streamtubeSolidPressureForces(shifted), /both original solid TE corners/);
  const missing = structuredClone(input); delete missing.layout.bodies[0].trailingEdge;
  assert.throws(() => streamtubeSolidPressureForces(missing), /joined leading and trailing edges/);
  for (const p of [NaN, Infinity, 0, -1]) {
    const bad = structuredClone(input); bad.flow.cells[2][1][0].interfacePressure.lower = p;
    assert.throws(() => streamtubeSolidPressureForces(bad), /interface pressure/);
  }
});
