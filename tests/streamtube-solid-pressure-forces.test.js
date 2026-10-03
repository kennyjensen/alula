// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeSolidPressureForces } from '../src/euler/streamtube-forces.js';

const close = (a, b, tolerance = 2e-12) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) < tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const keys = ['cx', 'cy', 'cl', 'cm', 'pressureIntegralDrag'];

// Diamond area=1 and centroid=(1,0). The fixture contains pressures only;
// it neither calls the Euler equations nor obtains forces from a BL speed.
function fixture(pressure = () => 3, { count = 1, alpha = 0 } = {}) {
  const nx = 4, tubes = Array(count + 1).fill(1);
  const bodies = Array.from({ length: count }, (_, body) => ({ leadingIndex: 1, trailingIndex: 3, element: count - body - 1 }));
  const lower = b => [{ x: -1, y: 3 * b }, { x: 0, y: 3 * b }, { x: 1, y: 3 * b - .5 }, { x: 2, y: 3 * b }, { x: 3, y: 3 * b }];
  const upper = b => lower(b).map((p, i) => ({ x: p.x, y: i === 2 ? 3 * b + .5 : p.y }));
  const nodes = tubes.map((_, group) => Array.from({ length: nx + 1 }, (_, i) => [
    group === 0 ? { x: i - 1, y: -2 } : upper(group - 1)[i],
    group === count ? { x: i - 1, y: 3 * count } : lower(group)[i]
  ]));
  const cells = Array.from({ length: nx - 1 }, () => tubes.map(() => [{ interfacePressure: { lower: 3, upper: 3 } }]));
  bodies.forEach((b, body) => {
    for (let i = b.leadingIndex; i <= b.trailingIndex; i++) {
      cells[i - 1][body][0].interfacePressure.upper = pressure(lower(body)[i], body, 'lower', i);
      cells[i - 1][body + 1][0].interfacePressure.lower = pressure(upper(body)[i], body, 'upper', i);
    }
  });
  return { flow: { nodes, cells }, layout: { nx, tubes, bodies }, conditions: { pInf: 2, alpha, flowModel: 'compressible' } };
}

test('constant physical pressure produces zero load on every closed body', () => {
  const input = fixture(() => 8.5, { count: 2, alpha: 13 }), before = structuredClone(input);
  const result = streamtubeSolidPressureForces({ ...input, referenceChord: 1.7, momentReference: { x: -.3, y: .8 } });
  for (const key of keys) { close(result[key], 0); result.perBody.forEach(body => close(body[key], 0)); }
  assert.deepEqual(result.perBody.map(body => body.element), [1, 0]);
  assert.equal(Object.hasOwn(result, 'cd'), false);
  assert.equal(result.dragKind, 'pressure-contribution'); assert.equal(result.includesSkinFriction, false);
  assert.equal(result.includesExitDefects, false); assert.equal(result.physicalAcceptance, false);
  assert.deepEqual(input, before);
});

test('linear physical pressure agrees with area-integral force and nose-up moment at arbitrary references', () => {
  const a = .3, b = -.2, alpha = 17;
  const input = fixture(p => 4 + a * p.x + b * p.y, { count: 2, alpha });
  for (const referenceChord of [.7, 2.3]) {
    const momentReference = { x: -.4, y: .8 };
    const result = streamtubeSolidPressureForces({ ...input, referenceChord, momentReference });
    const angle = alpha * Math.PI / 180;
    let totalMoment = 0;
    result.perBody.forEach((body, k) => {
      // Divergence theorem: F=-2*A*grad(p)/c, and the nose-up
      // moment is 2*A*((xc-x0)*p_y-(yc-y0)*p_x)/c^2.
      const cx = -2 * a / referenceChord, cy = -2 * b / referenceChord;
      const cm = 2 * ((1 - momentReference.x) * b - (3 * k - momentReference.y) * a) / referenceChord ** 2;
      close(body.cx, cx); close(body.cy, cy); close(body.cm, cm);
      close(body.cl, cy * Math.cos(angle) - cx * Math.sin(angle));
      close(body.pressureIntegralDrag, cx * Math.cos(angle) + cy * Math.sin(angle));
      close(body.enclosedSolidArea, 1); totalMoment += cm;
    });
    close(result.cm, totalMoment);
    for (const key of keys) close(result[key], result.perBody.reduce((sum, body) => sum + body[key], 0));
  }
});

test('physical pressure changes affect loads while BL speed and local Mach do not enter the integration', () => {
  const baseline = fixture(p => 4 + .2 * p.x), initial = streamtubeSolidPressureForces(baseline);
  const tagged = structuredClone(baseline);
  tagged.flow.boundaryLayer = { stations: [{ ue: NaN }, { ue: 900 }] };
  for (const mach of [.8, 1, 1.4]) {
    const q = 1.2, p = 2, enthalpy = q ** 2 / (.4 * mach ** 2), rho = 1.4 * p / (.4 * enthalpy);
    // Each tagged sample obeys p=(gamma-1)/gamma*rho*h and
    // M^2=q^2/((gamma-1)h), on either side of sonic speed.
    tagged.flow.sections = [[[{ rho, q, p, enthalpy, machSquared: mach ** 2 }]]];
    assert.deepEqual(streamtubeSolidPressureForces(tagged), initial);
  }
  // Change actual wall pressure while keeping the deliberately irrelevant
  // BL data identical: a new transverse gradient creates a lift change.
  const changed = fixture(p => 4 + .2 * p.x - .15 * p.y);
  changed.flow.boundaryLayer = tagged.flow.boundaryLayer;
  const after = streamtubeSolidPressureForces(changed);
  close(after.cy - initial.cy, .3); close(after.cl - initial.cl, .3);
});

test('displaced flow integrates the preserved solid contour and rejects a missing solid grid', () => {
  const input = fixture(p => 4 - .2 * p.y), baseline = streamtubeSolidPressureForces(input);
  const shifted = structuredClone(input);
  shifted.layout.displacedBoundaries = true;
  shifted.flow.undisplacedNodes = structuredClone(input.flow.nodes);
  shifted.flow.nodes[0][2][1].y -= .4; shifted.flow.nodes[1][2][0].y += .4;
  const result = streamtubeSolidPressureForces(shifted);
  for (const key of keys) assert.equal(result[key], baseline[key]);
  assert.equal(result.geometry, 'undisplaced-solid-grid');
  delete shifted.flow.undisplacedNodes;
  assert.throws(() => streamtubeSolidPressureForces(shifted), /undisplaced solid grid/);
  const flagged = structuredClone(input); flagged.flow.displacement = {};
  assert.throws(() => streamtubeSolidPressureForces(flagged), /undisplaced solid grid/);
});

test('unconverged leading-edge pressures use their explicit physical mean', () => {
  const input = fixture(() => 4);
  input.flow.cells[0][0][0].interfacePressure.upper = 3;
  input.flow.cells[0][1][0].interfacePressure.lower = 5;
  const result = streamtubeSolidPressureForces(input);
  for (const key of keys) close(result[key], 0);
  assert.equal(result.perBody[0].stagnationPressure, 4);
  assert.deepEqual(result.perBody[0].leadingPressures, { upper: 5, lower: 3 });
  assert.equal(result.perBody[0].leadingPressureMismatch, 2);
});

test('pressure-load adapter rejects malformed, nonphysical and unclosed input', () => {
  const input = fixture();
  for (const change of [{ referenceChord: 0 }, { referenceChord: Infinity }, { momentReference: { x: NaN, y: 0 } }])
    assert.throws(() => streamtubeSolidPressureForces({ ...input, ...change }), /reference/);
  for (const pInf of [0, -1, NaN]) assert.throws(() => streamtubeSolidPressureForces({ ...input,
    conditions: { ...input.conditions, pInf } }), /reference/);
  for (const pressure of [0, -1, NaN, Infinity]) {
    const bad = fixture(); bad.flow.cells[1][1][0].interfacePressure.lower = pressure;
    assert.throws(() => streamtubeSolidPressureForces(bad), /interface pressure/);
  }
  for (const mutate of [x => x.flow.nodes.pop(), x => x.flow.cells.pop(), x => { x.flow.nodes[0][0][0].x = NaN; },
    x => { x.layout.tubes[0] = 0; }]) {
    const bad = fixture(); mutate(bad); assert.throws(() => streamtubeSolidPressureForces(bad), /dimensions|nodes and cells/);
  }
  const open = fixture(); open.flow.nodes[0][3][1].y += .01;
  assert.throws(() => streamtubeSolidPressureForces(open), /joined leading and trailing edges/);
  const clockwise = fixture(); clockwise.flow.nodes[0][2][1].y = .5; clockwise.flow.nodes[1][2][0].y = -.5;
  assert.throws(() => streamtubeSolidPressureForces(clockwise), /positive counterclockwise area/);
});
