import test from 'node:test';
import assert from 'node:assert/strict';
import { mapCoupledBLProfiles, prepareCoupledCoarseProfile } from '../src/euler/tests/streamtube-coupled-coarse-profile.js';
import { createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { streamtubeBaseGeometry } from '../src/euler/streamtube-geometry.js';
import { createXfoilDeadAirGap } from '../src/viscous/xfoil-dead-air-gap.js';

const conditions = { lengthScale: 1, massScale: 1, center: { x: .5, y: 0 }, mach: .2, alpha: 2,
  gamma: 1.4, flowModel: 'compressible', streamwiseMode: 'isentropic', reynolds: 1e6, ncrit: 9,
  edgeMatching: 'section-velocity', transitionMode: 'automatic', pressureCorrectionFactor: .1,
  normalStencil: 'body-stations', stagnationMotion: 'walls-only', geometryDomain: 'positive-simple',
  wakeGeometry: 'independent-banks', wakeOutlet: 'banks' };
const body = element => ({ element, points: [[1, 0], [.5, .1], [0, 0], [.5, -.1], [1, 0]]
  .map(([x, y]) => ({ x: x + element * 2, y })) });
function source(order = [0]) {
  const d = { conditions: structuredClone(conditions), scale: .001, bodies: order.map(body),
    trips: order.map(() => [1, 1]), states: [], surfaces: [], wakes: [] };
  const put = p => { const id = d.states.length; d.states.push(p); return id; };
  order.forEach((element, b) => {
    const factor = element + 1;
    for (const side of ['upper', 'lower']) d.surfaces.push({ body: b, side, transition: 2, transitionS: .55,
      ids: [.1, .4, .7, 1].map((s, j) => put({ s, theta: .001 * factor * (1 + s), deltaStar: .002 * factor * (1 + s),
        ue: .2 + .8 * s, aux: [0, 2, .04, .05][j] })) });
    d.wakes.push({ body: b, ids: [1, 1.5, 2].map((s, j) => put({ s, theta: .004 * factor + .001 * factor * (s - 1),
      deltaStar: .008 * factor - .001 * factor * (s - 1), ue: .9 + .05 * (s - 1), aux: .05 + .01 * (s - 1) })) });
  });
  return d;
}
function target(from, order = from.bodies.map(b => b.element), surfaceArc = [.05, .1, .25, .5, .55, .6, .8, 1], wakeArc = [1, 1.25, 1.75, 2]) {
  const d = { conditions: structuredClone(from.conditions), scale: from.scale, bodies: order.map(body),
    trips: order.map(() => [1, 1]), coordinates: [], surfaces: [], wakes: [] };
  const put = s => { const id = d.coordinates.length; d.coordinates.push({ s }); return id; };
  order.forEach((element, b) => {
    for (const side of ['upper', 'lower']) d.surfaces.push({ body: b, side, ids: surfaceArc.map(put) });
    d.wakes.push({ body: b, ids: wakeArc.map(put) });
  });
  return d;
}
const state = (r, t, id) => ({ aux: r.initialBL[4 * id], theta: t.scale * r.initialBL[4 * id + 1],
  deltaStar: t.scale * r.initialBL[4 * id + 2], ue: r.initialBL[4 * id + 3] });
const close = (a, b, tolerance = 2e-17) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);

// Two distinct physical bases and independently sampled XICALC tails. These
// are analytic BL profiles on manufactured solids, not converged flows.
function withFiniteBases(d, elements = d.bodies.map(b => b.element)) {
  d.bodies = d.bodies.map(b => {
    if (!elements.includes(b.element)) return b;
    const h = .004 * (b.element + 1);
    return { element: b.element, trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 8 },
      points: [[1, h], [.8, .032], [.5, .075], [.2, .06], [0, 0], [.2, -.06], [.5, -.075], [.8, -.032], [1, -h], [1, 0], [1, h]]
        .map(([x, y]) => ({ x: x + 2 * b.element, y })) };
  });
  for (const wake of d.wakes) {
    const b = d.bodies[wake.body]; if (b.trailingEdge?.kind !== 'finite-base') continue;
    const base = streamtubeBaseGeometry([b], [createSurfaceContourCurve(b.points, b)])[0];
    const model = createXfoilDeadAirGap({ normalGap: base.width, upperDerivative: base.upperDerivative, lowerDerivative: base.lowerDerivative });
    const data = d.states ?? d.coordinates, start = data[wake.ids[0]].s, length = data[wake.ids.at(-1)].s - start;
    for (const id of wake.ids) {
      const p = data[id], fraction = (p.s - start) / length;
      p.wakeGap = model.at(p.s - start).gap;
      if (d.states) {
        p.theta = .001 * (1 + fraction) * (1 + b.element);
        p.deltaStar = .003 * (1 + .5 * fraction) * (1 + b.element) + p.wakeGap;
      }
    }
  }
  return d;
}

test('same-grid mapping preserves every physical field, phase, endpoint and independent caller data', () => {
  const s = source([0, 1]), t = target(s, [0, 1], [.1, .4, .7, 1], [1, 1.5, 2]);
  const before = structuredClone({ s, t }), r = mapCoupledBLProfiles({ source: s, target: t });
  s.states.forEach((p, id) => {
    assert.equal(r.initialBL[4 * id], p.aux); assert.equal(r.initialBL[4 * id + 1], p.theta / t.scale);
    assert.equal(r.initialBL[4 * id + 2], p.deltaStar / t.scale); assert.equal(r.initialBL[4 * id + 3], p.ue);
  });
  assert.deepEqual(r.transitionState, [2, 2, 2, 2]); assert.deepEqual({ s, t }, before);
  r.initialBL.fill(NaN); r.mapping.surfaces[0].rows[0].s = NaN;
  assert.deepEqual({ s, t }, before);
});

test('wake displacement chart is locked across profile transfer with legacy omission meaning fixed', () => {
  const s = source(), t = target(s), baseline = mapCoupledBLProfiles({ source: s, target: t });
  t.conditions.wakeDisplacementMotion = 'fixed';
  assert.deepEqual(mapCoupledBLProfiles({ source: s, target: t }), baseline);
  t.conditions.wakeDisplacementMotion = 'te-center';
  assert.throws(() => mapCoupledBLProfiles({ source: s, target: t }), /wakeDisplacementMotion/);
  s.conditions.wakeDisplacementMotion = 'te-center';
  assert.deepEqual(mapCoupledBLProfiles({ source: s, target: t }), baseline);
  t.conditions.wakeDisplacementMotion = 'unknown';
  assert.throws(() => mapCoupledBLProfiles({ source: s, target: t }), /physical\/model/);
});

test('nonuniform physical-arc interpolation is linear while the leading extension retains finite thickness and linear speed', () => {
  const s = source(), t = target(s), r = mapCoupledBLProfiles({ source: s, target: t });
  const surface = t.surfaces[0], first = state(r, t, surface.ids[0]);
  close(first.theta, .0011); close(first.deltaStar, .0022); close(first.ue, .14); assert.equal(first.aux, 0);
  for (const id of surface.ids.slice(1)) {
    const x = t.coordinates[id].s, p = state(r, t, id);
    close(p.theta, .001 * (1 + x)); close(p.deltaStar, .002 * (1 + x)); close(p.ue, .2 + .8 * x, 3e-16);
  }
  const stretched = target(s, [0], [.05, .1, .25, .5, .55, .6, .8, 1].map(x => 1.2 * x));
  const z = mapCoupledBLProfiles({ source: s, target: stretched });
  close(state(z, stretched, 0).theta, first.theta);
  close(state(z, stretched, 0).ue / stretched.coordinates[0].s, first.ue / t.coordinates[0].s / 1.2, 1e-14);
});

test('shared auxiliary slot is never blended between N and Ctau across natural transition', () => {
  const s = source(), t = target(s), r = mapCoupledBLProfiles({ source: s, target: t }), ids = t.surfaces[0].ids;
  assert.equal(r.transitionState[0], 4);
  assert.equal(state(r, t, ids[3]).aux, 2); // s=.5: laminar side of the bracket.
  assert.equal(state(r, t, ids[4]).aux, .04); // s=.55: transition/shear side.
  assert.equal(state(r, t, ids[5]).aux, .04);
  assert.ok(r.mapping.surfaces[0].rows.slice(3, 6).every(p => p.crossedAuxiliaryPhase));
});

test('wake mapping uses TE-to-outlet fraction and preserves merged endpoint values', () => {
  const s = source(), t = target(s, [0], [.1, .4, .7, 1], [2, 2.5, 3.5, 4]);
  const r = mapCoupledBLProfiles({ source: s, target: t }), ids = t.wakes[0].ids;
  close(state(r, t, ids[0]).theta, .004); close(state(r, t, ids.at(-1)).theta, .005);
  close(state(r, t, ids[1]).theta, .00425); close(state(r, t, ids[2]).deltaStar, .00725);
  const sum = t.surfaces.reduce((a, p) => a + state(r, t, p.ids.at(-1)).theta, 0);
  close(state(r, t, ids[0]).theta, sum);
  assert.equal(r.mapping.wakes[0].sourceLength, 1); assert.equal(r.mapping.wakes[0].targetLength, 2);
});

test('initial geometry guard choice may differ while actual Euler stencil and physical conditions must match', () => {
  const s = source(), t = target(s); t.conditions.geometryDomain = 'convex';
  const reference = mapCoupledBLProfiles({ source: s, target: target(s) });
  assert.deepEqual(mapCoupledBLProfiles({ source: s, target: t }).initialBL, reference.initialBL);
  t.conditions.normalStencil = 'centered';
  assert.throws(() => mapCoupledBLProfiles({ source: s, target: t }), /condition normalStencil/);
});

test('explicit element identities support reordered bodies, surfaces and wake arrays', () => {
  const s = source([0, 1]), t = target(s, [1, 0]);
  s.surfaces.reverse(); s.wakes.reverse(); t.surfaces.reverse(); t.wakes.reverse();
  const before = structuredClone({ s, t }), r = mapCoupledBLProfiles({ source: s, target: t });
  assert.deepEqual(r.mapping.bodyMap, [1, 0]);
  for (const p of t.surfaces) close(state(r, t, p.ids.at(-1)).theta, .002 * (t.bodies[p.body].element + 1));
  for (const p of t.wakes) close(state(r, t, p.ids[0]).theta, .004 * (t.bodies[p.body].element + 1));
  assert.deepEqual({ s, t }, before);
});

test('new leading turbulent station uses the supplied native closure rather than interpolated amplification', () => {
  const s = source();
  for (const p of s.surfaces) {
    p.transition = 0; p.transitionS = .02; p.ids.forEach((id, j) => { s.states[id].aux = .04 + .01 * j; });
  }
  const t = target(s); let calls = 0;
  assert.throws(() => mapCoupledBLProfiles({ source: s, target: t }), /native shear closure/);
  const r = mapCoupledBLProfiles({ source: s, target: t, transitionShear: p => {
    calls++; close(p.theta, .0011); close(p.ue, .14); assert.equal(p.aux, .03); return .081;
  } });
  assert.equal(calls, 2); assert.equal(r.initialBL[0], .081); assert.deepEqual(r.transitionState, [0, 0]);
});

test('unsupported controls, malformed bases, identity errors and incomplete station coverage fail without mutation', () => {
  const changes = [
    [x => { x.t.conditions.transitionMode = 'fixed-trip'; }, /automatic/],
    [x => { x.s.trips[0][0] = .9; }, /terminal/],
    [x => { x.t.conditions.mach = .3; }, /condition mach/],
    [x => { x.t.conditions.lengthScale = 2; }, /condition lengthScale/],
    [x => { x.t.conditions.ncrit = 4; }, /condition ncrit/],
    [x => { x.t.conditions.blThermodynamics = 'historical-common-isentrope'; }, /condition blThermodynamics/],
    [x => { delete x.t.conditions.normalStencil; }, /complete finite physical/],
    [x => { x.t.conditions.mach = NaN; }, /complete finite physical/],
    [x => { x.t.scale *= 2; }, /normalization/],
    [x => { x.t.bodies[0].trailingEdge = { kind: 'finite-base' }; }, /Finite-base topology/],
    [x => { x.s.states.at(-1).wakeGap = .001; }, /zero base gap/],
    [x => { x.t.bodies[0].points[1].y += .1; }, /contour/],
    [x => { delete x.t.bodies[0].element; }, /identities/],
    [x => { x.t.wakes[0].ids[0] = x.t.surfaces[0].ids[0]; }, /overlap/],
    [x => { x.t.coordinates.push({ s: 3 }); }, /omitted/],
    [x => { x.t.coordinates[1].s = x.t.coordinates[0].s; }, /increase/],
    [x => { x.s.surfaces[0].transitionS = NaN; }, /transition/],
    [x => { x.s.states[0].ue = 0; }, /physical profile/],
  ];
  for (const [change, expression] of changes) {
    const x = { s: source() }; x.t = target(x.s); change(x); const before = structuredClone(x);
    assert.throws(() => mapCoupledBLProfiles({ source: x.s, target: x.t }), expression); assert.deepEqual(x, before);
  }
});

test('finite-base mapping interpolates fluid displacement and reconstructs the independently resolved target dead-air tail', () => {
  const s = source();
  s.wakes[0].ids.forEach((id, j) => { s.states[id].s = [1, 1.01, 1.02][j]; });
  withFiniteBases(s);
  const t = withFiniteBases(target(s, [0], [.1, .4, .7, 1], [2, 2.002, 2.007, 2.015, 2.04]));
  const before = structuredClone({ s, t }), r = mapCoupledBLProfiles({ source: s, target: t });
  for (const id of t.wakes[0].ids) {
    const fraction = (t.coordinates[id].s - 2) / .04, p = state(r, t, id);
    close(p.theta, .001 * (1 + fraction));
    close(p.deltaStar - t.coordinates[id].wakeGap, .003 * (1 + .5 * fraction));
    assert.ok(p.deltaStar - t.coordinates[id].wakeGap > p.theta);
  }
  // This station is inside the target's cubic tail. Interpolating the
  // source TOTAL displacement would preserve the wrong geometric width.
  const id = t.wakes[0].ids[2], row = r.mapping.wakes[0].rows[2];
  const a = s.states[s.wakes[0].ids[row.left]], b = s.states[s.wakes[0].ids[row.right]];
  assert.ok(Math.abs(state(r, t, id).deltaStar - (a.deltaStar + row.t * (b.deltaStar - a.deltaStar))) > .001);
  assert.equal(r.fluidWakeDisplacement.length, t.wakes[0].ids.length);
  assert.deepEqual({ s, t }, before);
});

test('same-grid finite-base mapping preserves total fields exactly and mixed bodies keep separate material wakes', () => {
  const s = source([0, 1]);
  for (const w of s.wakes) w.ids.forEach((id, j) => { s.states[id].s = [1, 1.004, 1.012][j]; });
  withFiniteBases(s);
  const t = withFiniteBases(target(s, [1, 0], [.1, .4, .7, 1], [1, 1.004, 1.012]));
  const r = mapCoupledBLProfiles({ source: s, target: t });
  for (const w of t.wakes) {
    const parent = s.wakes.find(p => s.bodies[p.body].element === t.bodies[w.body].element);
    w.ids.forEach((id, j) => {
      const p = s.states[parent.ids[j]];
      assert.equal(r.initialBL[4 * id + 1], p.theta / t.scale);
      assert.equal(r.initialBL[4 * id + 2], p.deltaStar / t.scale);
      assert.equal(r.initialBL[4 * id], p.aux); assert.equal(r.initialBL[4 * id + 3], p.ue);
    });
  }
  assert.deepEqual(r.transitionState, [2, 2, 2, 2]);
  const mixed = withFiniteBases(source([0, 1]), [1]);
  const targetMixed = withFiniteBases(target(mixed, [1, 0]), [1]);
  const transferred = mapCoupledBLProfiles({ source: mixed, target: targetMixed });
  assert.equal(transferred.fluidWakeDisplacement.length, targetMixed.wakes[0].ids.length);
  assert.ok(transferred.fluidWakeDisplacement.every(p => targetMixed.wakes[0].ids.includes(p.id)));
});

test('finite-base transfer rejects missing geometric gaps, fluid-domain failures, incompatible bank mode and changed TE corners', () => {
  const mutations = [
    [x => { delete x.s.states[x.s.wakes[0].ids[1]].wakeGap; }, /every decoded/],
    [x => { x.t.coordinates[x.t.wakes[0].ids[1]].wakeGap = -1; }, /nonnegative/],
    [x => { x.s.states[x.s.wakes[0].ids[1]].wakeGap = NaN; }, /nonnegative/],
    [x => { const p = x.s.states[x.s.wakes[0].ids[0]]; p.deltaStar = p.wakeGap + .5 * p.theta; }, /physical profile/],
    [x => { x.s.states[x.s.wakes[0].ids[0]].wakeGap = 0; }, /positive solid-base/],
    [x => { x.t.conditions.wakeGeometry = 'centerline'; }, /independent/],
    [x => { x.t.bodies[0].trailingEdge = { kind: 'finite-base', upperIndex: 1, lowerIndex: 7 }; }, /trailing-edge topology/],
    [x => { x.t.coordinates[x.t.wakes[0].ids[0]].wakeGap *= 2; }, /base width/],
  ];
  for (const [mutate, expression] of mutations) {
    const s = withFiniteBases(source()), t = withFiniteBases(target(s)), x = { s, t }; mutate(x);
    const before = structuredClone(x);
    assert.throws(() => mapCoupledBLProfiles({ source: s, target: t }), expression); assert.deepEqual(x, before);
  }
});

test('wrapper rejects unaccepted sources before constructing any target or evaluating the source', () => {
  let evaluations = 0; const sourceSystem = { bl: {}, evaluate: () => { evaluations++; throw new Error('unexpected'); } };
  for (const sourceResult of [undefined, { converged: false }, { converged: true, mesh: { quality: { valid: true } }, families: { euler: 1, boundaryLayer: 0, edgeMatching: 0 } }])
    assert.throws(() => prepareCoupledCoarseProfile({ sourceSystem, sourceResult }), /accepted source/);
  assert.equal(evaluations, 0);
});
