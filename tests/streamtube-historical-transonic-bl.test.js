// SPDX-License-Identifier: GPL-2.0-or-later
// Synthetic local states and frozen coupled evaluations, not solved shocks.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { createPhysicalBLDomain } from '../src/viscous/physical-domain.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody, coupledStreamtubeResult } from '../src/euler/streamtube-coupled.js';
import { createStreamtubeBoundaryLayers, initializeStreamtubeBLStation } from '../src/euler/streamtube-boundary-layers.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const model = 'historical-common-isentrope', gamma = 1.4, reynolds = 1e6;
const upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
const close = (a, b, limit = 2e-12) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) < limit * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const qAtMach = (edgeMach, mach) => Math.sqrt((gamma - 1) * edgeMach ** 2
  * (1 / ((gamma - 1) * mach ** 2) + .5) / (1 + .5 * (gamma - 1) * edgeMach ** 2));
const station = (edgeMach, mach, { s = .2, theta = .0002, hk = 2.2, aux = 0 } = {}) => ({
  s, theta, deltaStar: theta * (hk * (1 + .113 * edgeMach ** 2) + .29 * edgeMach ** 2),
  ue: qAtMach(edgeMach, mach), aux
});

test('historic BL gas and physical domain remain regular below, at and above local sonic speed', () => {
  const mach = .4, kernel = createIntegralKernel({ reynolds, mach, gamma, exactJacobian: true });
  const domain = createPhysicalBLDomain({ mach, gamma });
  for (const edgeMach of [.8, 1, 1.2, 1.6]) {
    const state = station(edgeMach, mach), p = kernel.station(state);
    const t = 1 + .5 * (gamma - 1) * mach ** 2 * (1 - state.ue ** 2);
    const t0 = 1 + .5 * (gamma - 1) * mach ** 2, sutherland = .35 * t0;
    close(p.machSquared, edgeMach ** 2); close(p.rawHk, 2.2);
    close(p.rho, t ** (1 / (gamma - 1)));
    close(p.viscosity, t ** 1.5 * (1 + sutherland) / (t + sutherland) / reynolds);
    close(p.reTheta, p.rho * state.ue * state.theta / p.viscosity);
    assert.ok(domain(state).enthalpy > 0 && domain(state).shape > 0);
  }
  const valid = station(1.2, mach), h0 = 1 / ((gamma - 1) * mach ** 2) + .5;
  const thermal = { ...valid, ue: Math.sqrt(2 * h0) * 1.001 };
  assert.ok(domain(thermal).enthalpy < 0);
  assert.throws(() => kernel.station(thermal), /compressible BL edge/);
  const shape = { ...valid, deltaStar: 1.01 * valid.theta };
  assert.ok(domain(shape).shape < 0); assert.throws(() => kernel.station(shape), /compressible BL edge/);
});

test('unchanged native laminar, turbulent and wake blocks have correct derivatives through sonic speed', () => {
  const mach = .4, kernel = createIntegralKernel({ reynolds, mach, gamma, exactJacobian: true });
  for (const regime of ['laminar', 'turbulent', 'wake']) {
    const hk = regime === 'laminar' ? 2.2 : regime === 'turbulent' ? 1.6 : 1.3;
    const input = { regime, upstream: station(.95, mach, { hk, aux: regime === 'laminar' ? .2 : .03 }),
      downstream: station(1.05, mach, { s: .23, theta: .00021, hk, aux: regime === 'laminar' ? .25 : .031 }) };
    const value = kernel.interval(input);
    for (const side of ['upstream', 'downstream']) for (const [k, key] of ['aux', 'theta', 'deltaStar', 'ue', 's'].entries()) {
      const base = input[side][key], h = 1e-5 * Math.max(Math.abs(base), key === 'aux' ? .01 : 1e-5);
      const samples = [-2, -1, 1, 2].map(f => kernel.interval({ ...input, [side]: { ...input[side], [key]: base + f * h } }).residual);
      for (let row = 0; row < 3; row++) {
        const fd = (samples[0][row] - 8 * samples[1][row] + 8 * samples[2][row] - samples[3][row]) / (12 * h);
        close(value[side][row][k], fd, 2e-6);
      }
    }
  }
});

// Prescribe a smooth positive grid, independent BL unknowns and a small
// band of physical Euler speeds. No constructor may run its missing-seed
// fallback solve or the BL initializer; all state arrays are supplied.
function frozenFixture(edgeMach = 1.2, mach = .5) {
  const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2, mach }),
    streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 }, upwind,
    normalStencil: 'body-stations', stagnationMotion: 'walls-only' };
  const euler = createStreamtubeBodySystem(input), x = euler.initial.slice(), base = euler.evaluate(x);
  const body = euler.layout.bodies[0], speed = qAtMach(edgeMach, mach);
  for (const i of [body.leadingIndex + 1, body.leadingIndex + 2]) {
    const g = 1, j = 0, section = base.sections[i][g][j];
    x[euler.layout.densityIndex(i, g, j)] = Math.log(section.rho * section.q / speed);
  }
  const bl = createStreamtubeBoundaryLayers(euler, x, { reynolds });
  const initialBL = new Float64Array(4 * bl.stations.length);
  for (const s of bl.stations) {
    const theta = s.kind === 'wake' ? 8e-7 : 4e-7;
    const state = station(edgeMach, mach, { theta, hk: 2.2 });
    const turbulent = ['transition', 'turbulent', 'trailing-edge', 'wake'].includes(s.regime);
    initialBL.set([turbulent ? .03 : 0, theta / bl.scale, state.deltaStar / bl.scale, state.ue], 4 * s.id);
  }
  return { input, options: { reynolds, edgeMatching: 'section-velocity', initialEuler: { x, nodes: base.nodes }, initialBL } };
}

test('explicit hybrid coupling admits supersonic physical edges while the default retains its sonic rejection', () => {
  const fixture = frozenFixture(), before = structuredClone(fixture);
  const legacy = createCoupledStreamtubeBody(fixture.input, fixture.options);
  assert.throws(() => legacy.evaluate(legacy.initial), /Corrected BL edge velocity is sonic/);
  const explicit = createCoupledStreamtubeBody(fixture.input, { ...fixture.options, blThermodynamics: model });
  const value = explicit.evaluate(explicit.initial);
  assert.equal(explicit.conditions.blThermodynamics, model);
  assert.ok(value.residual.every(Number.isFinite));
  assert.ok(value.edges.some(edge => edge.sides?.some(side => side.meanMachSquared > 1)));
  assert.ok(explicit.bl.kernel.station(value.layers.states[0]).machSquared > 1);
  // Actual outer density is retained independently of the common-isentrope
  // BL closure. This artificial state need not satisfy either set of rows.
  const physical = value.outer.sections[explicit.euler.layout.bodies[0].leadingIndex + 1][1][0];
  const gas = explicit.bl.kernel.station(value.layers.states[0]);
  assert.ok(Math.abs(physical.rho - gas.rho) > .01);
  assert.deepEqual(fixture, before);
  const bad = explicit.initial.slice(), k = explicit.ne;
  bad[k + 3] = Math.sqrt(2 * explicit.euler.conditions.h0) * 1.001;
  assert.equal(explicit.admissible(bad), false, 'Positive thermal energy remains mandatory.');
  bad.set(explicit.initial); bad[k + 2] = 1.01 * bad[k + 1];
  assert.equal(explicit.admissible(bad), false, 'Raw kinematic shape remains mandatory above Mach one.');
});

test('explicit historical option is an exact subsonic residual and Jacobian replay', () => {
  const fixture = frozenFixture(.6, .2), a = createCoupledStreamtubeBody(fixture.input, fixture.options);
  const b = createCoupledStreamtubeBody(fixture.input, { ...fixture.options, blThermodynamics: model });
  const av = a.evaluate(a.initial), bv = b.evaluate(b.initial);
  assert.deepEqual(a.initial, b.initial); assert.deepEqual(av.residual, bv.residual);
  assert.deepEqual(av.outer.nodes, bv.outer.nodes); assert.deepEqual(av.layers.states, bv.layers.states);
  assert.deepEqual(a.jacobian(a.initial), b.jacobian(b.initial));
  assert.equal(Object.hasOwn(a.conditions, 'blThermodynamics'), false);
  assert.equal(Object.hasOwn(coupledStreamtubeResult(a, a.initial), 'edgeThermodynamics'), false);
});

test('frozen supersonic coupled Jacobian matches all mixed fields without changing the state', t => {
  const fixture = frozenFixture(), system = createCoupledStreamtubeBody(fixture.input,
    { ...fixture.options, blThermodynamics: model }), x = system.initial;
  assert.ok(system.n < 300, 'Keep the transonic coupled derivative fixture small.');
  const matrix = system.jacobian(x), before = system.evaluate(x), checks = [];
  for (const direction of [0, 1, 2]) {
    const d = x.map((v, i) => (i < system.ne ? .001 : .01 * Math.max(.01, Math.abs(v)))
      * Math.sin((.73 + .17 * direction) * i + .2));
    const exact = sparseProduct(matrix, d);
    for (const h of [4e-6, 2e-6]) {
      const values = [-2, -1, 1, 2].map(f => system.residual(x.map((v, i) => v + f * h * d[i])));
      let error = 0, worst = null;
      for (let row = 0; row < system.n; row++) {
        const fd = (values[0][row] - 8 * values[1][row] + 8 * values[2][row] - values[3][row]) / (12 * h);
        const relative = Math.abs(exact[row] - fd) / Math.max(1, Math.abs(exact[row]), Math.abs(fd));
        if (relative > error) { error = relative; worst = { row, exact: exact[row], fd }; }
      }
      checks.push({ direction, h, error, worst });
      assert.ok(error < 5e-6, JSON.stringify(checks.at(-1)));
    }
  }
  assert.deepEqual(system.evaluate(x).residual, before.residual);
  t.diagnostic(JSON.stringify({ unknowns: system.n, checks }));
});

test('a supersonic result exports distinct physical bank gas and historical BL gas without another solve', () => {
  const fixture = frozenFixture(), system = createCoupledStreamtubeBody(fixture.input,
    { ...fixture.options, blThermodynamics: model });
  const before = system.initial.slice(), result = coupledStreamtubeResult(system, system.initial);
  assert.deepEqual(system.initial, before); assert.equal(result.edgeThermodynamics.model, model);
  assert.equal(result.edgeThermodynamics.stations.length, system.bl.stations.length);
  assert.ok(result.edgeThermodynamics.maxDensityMismatch > .01);
  assert.ok(result.edgeThermodynamics.stations.some(s => s.boundaryLayer.machSquared > 1));
  for (const s of result.edgeThermodynamics.stations) {
    assert.equal(s.banks.length, s.kind === 'surface' ? 1 : 2);
    assert.ok(s.banks.every(bank => bank.sections.length === 2 && bank.sections.every(gas =>
      [gas.rho, gas.q, gas.p, gas.enthalpy, gas.entropyOverR].every(Number.isFinite))));
  }
  assert.match(result.limitations, /historical common-isentrope/);
  assert.match(result.limitations, /unequal-entropy/);
  assert.match(result.limitations, /unvalidated/);
});

test('historical mode rejects unsupported Euler and edge-matching combinations before initialization', () => {
  const input = { streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 }, upwind };
  for (const edgeMatching of ['pressure', 'section-velocity-distance'])
    assert.throws(() => createCoupledStreamtubeBody(input, { blThermodynamics: model, edgeMatching }), /arithmetic section-velocity/);
  for (const change of [{ streamwiseMode: 'isentropic' }, { streamwiseMode: 'momentum' }, { upwind: undefined }, { flowModel: 'incompressible' }])
    assert.throws(() => createCoupledStreamtubeBody({ ...input, ...change }, { blThermodynamics: model, edgeMatching: 'section-velocity' }), /Historical transonic BL requires/);
  assert.throws(() => createCoupledStreamtubeBody(input, { blThermodynamics: 'variable-entropy' }), /Unknown coupled BL/);
  assert.throws(() => initializeStreamtubeBLStation({ allowSupersonicEdge: 'yes' }), /Invalid BL edge-domain/);
});

test('opting into the thermal-domain inverse initializer leaves its established subsonic solution unchanged', () => {
  const kernel = createIntegralKernel({ reynolds, mach: .2 });
  // Retained native inverse case already documented in
  // streamtube-bl-initial-state.test.js; only a scalar initializer is run.
  const input = { upstream: { s: 1.0167217420548411, aux: .046435413000181835,
    theta: .004884672338085486, deltaStar: .009117906695388468, ue: .8802286950132129 },
    s: 1.02658439442974, ue: .7855043882457006, regime: 'turbulent', tripS: .07058146487518047,
    reynolds, interval: data => kernel.interval(data), properties: state => kernel.station(state, 'turbulent') };
  const a = initializeStreamtubeBLStation(input), b = initializeStreamtubeBLStation({ ...input, allowSupersonicEdge: true });
  assert.equal(a.mode, 'inverse'); assert.deepEqual(a, b);
});
