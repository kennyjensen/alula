import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createStreamtubeBoundaryLayers } from '../src/euler/streamtube-boundary-layers.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { scaleStreamtubeBLThicknesses } from '../src/euler/streamtube-coupled-initializer.js';
import { createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { createXfoilDeadAirGap } from '../src/viscous/xfoil-dead-air-gap.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

// Manufactured finite-base geometry, not a modification of measured data.
// Three near-wake stations resolve the prescribed cubic; no PDE is solved.
function fixture(finite = true) {
  // A low-Mach reservoir keeps the deliberately curved, very short first
  // wake volumes inside the Euler side-pressure domain without a flow solve.
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 2, contourPanels: 40, mach: .03 });
  if (finite) {
    const body = input.bodies[0], lowerIndex = body.points.length - 1;
    const points = body.points.map((p, i) => ({ x: p.x,
      y: p.y + (i < lowerIndex / 2 ? 1 : -1) * .002 * p.x }));
    body.points = [...points, { x: 1, y: 0 }, { ...points[0] }];
    body.trailingEdge = { kind: 'finite-base', upperIndex: 0, lowerIndex };
    body.stagnationParameter = createSurfaceContourCurve(body.points, body).length / 2;
    for (const row of [input.outerLower, input.outerUpper, ...input.cutPaths])
      [.0015, .004, .007].forEach((dx, k) => { row[body.trailingIndex + k + 1].x = 1 + dx; });
  }
  input.wakeGeometry = 'independent-banks';
  input.wakeOutlet = 'banks';
  input.displacement = { surfaces: [{ upper: Array(5).fill(0), lower: Array(5).fill(0) }],
    wakes: [Array(input.outerLower.length - 1 - input.bodies[0].trailingIndex).fill(finite ? .004 : 0)] };
  return input;
}

function prepared(finite = true) {
  const input = fixture(finite), euler = createStreamtubeBodySystem(input);
  const state = euler.initial.slice();
  // Curved wake centers give nonzero arc-length sensitivities in the normal
  // chart, independent of the kernel derivative that is being checked.
  if (finite) for (const position of euler.layout.positions) if (position.i > input.bodies[0].trailingIndex)
    state[position.column] += .00015 * Math.sin(position.i * 1.7);
  const bl = createStreamtubeBoundaryLayers(euler, state);
  const geo = bl.geometry(state), values = new Float64Array(4 * bl.stations.length);
  for (const s of bl.stations) {
    const theta = s.kind === 'wake' ? .0004 : .0002;
    const delta = 2.5 * theta + (geo.coordinates[s.id].wakeGap ?? 0);
    values.set([['similarity', 'laminar'].includes(s.regime) ? .2 : .03, theta / bl.scale, delta / bl.scale, 1 + .008 * Math.sin(s.id)], 4 * s.id);
  }
  return { input, euler, state, bl, geo, values };
}

const close = (a, b, tol = 1e-11) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const fd4 = (evaluate, x, d, h) => {
  const at = scale => evaluate(x.map((v, k) => v + scale * h * d[k]));
  const pp = at(2), p = at(1), m = at(-1), mm = at(-2);
  return p.map((_, k) => (-pp[k] + 8 * p[k] - 8 * m[k] + mm[k]) / (12 * h));
};

test('native finite-gap partial follows BLPRV and explicit HWA, independently checked by FD4', t => {
  const kernel = createIntegralKernel({ mach: .3, exactJacobian: true });
  const input = { regime: 'wake', upstream: { s: 1, ue: 1.04, theta: .0005, deltaStar: .0017, wakeGap: .0004, aux: .035 },
    downstream: { s: 1.02, ue: 1.01, theta: .00053, deltaStar: .00168, wakeGap: .0003, aux: .032 } };
  const v = kernel.interval(input); let maximum = 0;
  for (const side of ['upstream', 'downstream']) {
    const d = v[side].map(row => -row[2]), explicit = Math.log(input.downstream.ue / input.upstream.ue) / (2 * input[side].theta);
    d[1] += explicit; d[2] -= explicit;
    const h = 2e-8, at = q => kernel.interval({ ...input, [side]: { ...input[side], wakeGap: input[side].wakeGap + q * h } }).residual;
    const pp = at(2), p = at(1), m = at(-1), mm = at(-2);
    d.forEach((exact, k) => { const fd = (-pp[k] + 8 * p[k] - 8 * m[k] + mm[k]) / (12 * h);
      maximum = Math.max(maximum, Math.abs(fd - exact) / Math.max(1, Math.abs(fd), Math.abs(exact))); close(exact, fd, 2e-8); });
  }
  t.diagnostic(JSON.stringify({ maximum }));
});

test('finite-base BL geometry uses physical TE midpoint distance, total displacement and an invariant base contribution', () => {
  const { euler, state, bl, geo, values } = prepared(), base = euler.baseGeometry[0], w = bl.wakes[0];
  const model = createXfoilDeadAirGap({ normalGap: base.width, upperDerivative: base.upperDerivative, lowerDerivative: base.lowerDerivative });
  const nodes = euler.decode(state).nodes; let previous = base.center, distance = 0;
  for (const id of w.ids) {
    const i = bl.stations[id].i;
    if (id !== w.ids[0]) {
      const p = { x: .5 * (nodes[0][i].at(-1).x + nodes[1][i][0].x), y: .5 * (nodes[0][i].at(-1).y + nodes[1][i][0].y) };
      distance += Math.hypot(p.x - previous.x, p.y - previous.y); previous = p;
    }
    close(geo.coordinates[id].wakeDistance * euler.conditions.lengthScale, distance, 1e-15);
    close(geo.coordinates[id].wakeGap * euler.conditions.lengthScale, model.at(distance).gap, 1e-15);
  }
  const v = bl.evaluate(values, state), upper = bl.surfaces[0].ids.at(-1), lower = bl.surfaces[1].ids.at(-1), id = w.ids[0];
  close(v.residual[4 * id + 2], 0, 1e-13);
  close(v.states[id].deltaStar - v.states[upper].deltaStar - v.states[lower].deltaStar, base.width / euler.conditions.lengthScale, 1e-15);
  const half = scaleStreamtubeBLThicknesses(bl, values, .125);
  for (const s of bl.stations) {
    const k = 4 * s.id, gap = bl.initialWakeGaps[s.id] / bl.scale;
    close(half[k + 2] - gap, .125 * (values[k + 2] - gap), 1e-15);
    assert.equal(half[k + 1], .125 * values[k + 1]); assert.equal(half[k], values[k]); assert.equal(half[k + 3], values[k + 3]);
  }
  assert.deepEqual(scaleStreamtubeBLThicknesses(bl, values, 1), values);
});

test('finite-base complete coupled Jacobian and physical-domain gradients include moving dead-air arc', t => {
  const p = prepared(), input = { ...p.input }; delete input.displacement;
  const system = createCoupledStreamtubeBody(input, { initialEuler: { x: p.state, nodes: p.euler.decode(p.state).nodes }, initialBL: p.values,
    edgeMatching: 'section-velocity' });
  const x = system.initial, before = x.slice(), j = system.jacobian(x, { sparse: false });
  const geo = system.bl.geometry(x.subarray(0, system.ne), true);
  assert.ok(geo.coordinates.some(p => [...p.wakeGapDerivatives?.values() ?? []].some(v => Math.abs(v) > 1e-6)));
  let maximum = 0, constraintMaximum = 0;
  for (const phase of [.53, 1.07, 1.83]) {
    const direction = x.map((v, k) => Math.sin((k + 1) * phase) * (k < system.ne ? .03 : Math.max(.01, Math.abs(v))));
    const fd = fd4(a => system.residual(a), x, direction, 2e-7);
    for (let row = 0; row < system.n; row++) {
      let exact = 0; for (let col = 0; col < system.n; col++) exact += j[row * system.n + col] * direction[col];
      maximum = Math.max(maximum, Math.abs(exact - fd[row]) / Math.max(1, Math.abs(exact), Math.abs(fd[row])));
      close(exact, fd[row], 5e-6);
    }
    const constraints = system.stepConstraints(x), cf = fd4(a => system.constraintValues(a), x, direction, 2e-6);
    constraints.forEach((c, row) => {
      if (c.kind !== 'kinematic-shape') return;
      let exact = 0; for (const [col, d] of c.gradient) exact += d * direction[col];
      constraintMaximum = Math.max(constraintMaximum, Math.abs(exact - cf[row])); close(exact, cf[row], 5e-7);
    });
  }
  assert.deepEqual(x, before);
  const bad = x.slice(), id = system.bl.wakes[0].ids[0], k = system.ne + 4 * id;
  bad[k + 2] = .5 * geo.coordinates[id].wakeGap / system.bl.scale + bad[k + 1];
  assert.equal(system.bl.admissible(bad.subarray(system.ne)), true, 'legacy positive total-H precheck alone would admit this state');
  assert.equal(system.admissible(bad), false, 'physical fluid H must subtract the prescribed gap');
  t.diagnostic(JSON.stringify({ unknowns: system.n, maximum, constraintMaximum }));
});

test('one local finite-base wake initializer closes native TE matching without a PDE solve', () => {
  const p = prepared(), values = p.bl.initialize(() => 1), value = p.bl.evaluate(values, p.state);
  for (const w of p.bl.wakes) for (const id of w.ids) {
    const s = value.states[id]; assert.ok(s.deltaStar - s.wakeGap > s.theta);
    for (let k = 0; k < 3; k++) close(value.residual[4 * id + k], 0, 3e-9);
  }
});

for (const finite of [false, true]) for (const chart of ['fixed', 'te-center'])
test(`${finite ? 'finite' : 'sharp'} wake starts at the displaced TE in the ${chart} chart`, () => {
  const p = prepared(finite), input = { ...p.input, wakeDisplacementMotion: chart }; delete input.displacement;
  const values = p.values.slice();
  // Unequal surface thicknesses move the wake origin away from the solid TE.
  const upper = p.bl.surfaces[0].ids.at(-1), lower = p.bl.surfaces[1].ids.at(-1);
  values[4 * upper + 2] *= 2;
  const system = createCoupledStreamtubeBody(input, { initialEuler: { x: p.state, nodes: p.euler.decode(p.state).nodes },
    initialBL: values, edgeMatching: 'section-velocity' });
  const x = system.initial, ne = system.ne, w = system.bl.wakes[0];
  const sample = state => {
    system.euler.setDisplacement(system.bl.thicknesses(state.subarray(ne)));
    const flow = system.euler.decode(state.subarray(0, ne)), geo = system.bl.geometry(state.subarray(0, ne), true);
    const center = i => ({ x: .5 * (flow.nodes[0][i].at(-1).x + flow.nodes[1][i][0].x),
      y: .5 * (flow.nodes[0][i].at(-1).y + flow.nodes[1][i][0].y) });
    let previous = center(system.euler.layout.bodies[0].trailingIndex), distance = 0;
    for (const id of w.ids.slice(1)) {
      const point = center(system.bl.stations[id].i);
      distance += Math.hypot(point.x - previous.x, point.y - previous.y); previous = point;
      close(geo.coordinates[id].s - geo.coordinates[w.ids[0]].s, distance / system.euler.conditions.lengthScale, 2e-15);
    }
    return geo;
  };
  const geo = sample(x), id = w.ids[1];
  for (const station of [upper, lower]) {
    const col = ne + 4 * station + 2, h = 1e-4;
    const delta = new Float64Array(x.length); delta[col] = 1;
    const fd = fd4(state => [sample(state).coordinates[id].s], x, delta, h)[0];
    const exact = geo.coordinates[id].derivatives.get(col) ?? 0;
    close(exact, fd, 3e-11);
    if (chart === 'te-center') close(exact, 0, 1e-14, 'rigid wake translation preserves its arclength');
    else assert.ok(Math.abs(exact) > 1e-7, 'a fixed wake has a moving TE origin');
  }
});
