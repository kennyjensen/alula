import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { initializeStreamtubeBLStation } from '../src/euler/streamtube-boundary-layers.js';

// Retained upper-surface transition station from the 320-contour-panel,
// 16-interval/7-tube NACA 0012 reference comparison at alpha 2, M .2, Re 1e6.
// The earlier initializer copied upstream N=1.999 into turbulent Ctau and
// approached a natural-transition branch before a Jacobian trial aborted.
const upstream = { s: .0390818242369273, aux: 1.9991613738624674,
  theta: .00009181911107807067, deltaStar: .00027697074290347585, ue: 1.3115899921511796 };
const input = { upstream, s: .08651547339850425, ue: 1.3424578017745805,
  regime: 'transition', tripS: .07058146487518047, reynolds: 1e6 };
const setup = () => {
  const kernel = createIntegralKernel({ reynolds: 1e6, mach: .2, ncrit: 9 });
  const interval = data => {
    const r = kernel.interval(data);
    if (data.regime === 'transition' && !r.transition.forced) throw new Error('Natural transition precedes fixed trip.');
    return r;
  };
  return { kernel, interval };
};

test('a retained transition interval initializes turbulent shear independently of laminar amplification', () => {
  const { kernel, interval } = setup(), before = structuredClone(input);
  const r = initializeStreamtubeBLStation({ ...input, interval });
  assert.equal(r.converged, true); assert.ok(r.state.aux > 0 && r.state.aux < .1);
  assert.deepEqual(input, before);
  const check = kernel.interval({ upstream, downstream: r.state, regime: input.regime, tripS: input.tripS });
  assert.equal(check.transition.forced, true);
  assert.ok(Math.max(...check.residual.map(Math.abs)) < 1e-10);
  assert.ok(Math.abs(check.transition.s - input.tripS) < 1e-14);
});

test('BL initialization retains true natural-transition rejection instead of forcing the selected trip', () => {
  const { interval } = setup();
  assert.throws(() => initializeStreamtubeBLStation({ ...input, upstream: { ...upstream, aux: 10 }, interval }), /Natural transition/);
});

test('the retained upper trailing-edge direct failure admits a native inverse initialization without prescribing the final edge velocity', () => {
  const kernel = createIntegralKernel({ reynolds: 1e6, mach: .2, ncrit: 9 });
  const data = { upstream: { s: 1.0167217420548411, aux: .046435413000181835,
    theta: .004884672338085486, deltaStar: .009117906695388468, ue: .8802286950132129 },
    s: 1.02658439442974, ue: .7855043882457006, regime: 'turbulent',
    tripS: .07058146487518047, reynolds: 1e6, interval: input => kernel.interval(input) };
  assert.throws(() => initializeStreamtubeBLStation(data), /line search failed/);
  const r = initializeStreamtubeBLStation({ ...data, properties: s => kernel.station(s, 'turbulent') });
  assert.equal(r.mode, 'inverse'); assert.equal(r.directReason, 'line search failed');
  const state = kernel.station(r.state, 'turbulent');
  assert.ok(Math.abs(state.hk - 2.5) < 1e-10); assert.ok(state.machSquared < 1);
  assert.ok(Math.abs(r.state.ue - data.ue) > .01, 'initializer must allow a pressure-matching residual for the global solve');
  const block = kernel.interval({ upstream: data.upstream, downstream: r.state, regime: data.regime });
  assert.ok(Math.max(...block.residual.map(Math.abs)) < 1e-10);
});

test('inverse initialization constrains raw Hk and escapes the retained flap station closure clamp', () => {
  const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/streamtube-laminar-inverse.json', import.meta.url)));
  for (const [path, hash] of Object.entries(fixture.native.sourceHashes))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${path}`, import.meta.url))).digest('hex'), hash);
  const kernel = createIntegralKernel({ ...fixture.parameters, exactJacobian: true });
  const before = structuredClone(fixture), p = kernel.station(fixture.downstream);
  assert.equal(p.hk, 1.05); assert.ok(p.rawHk > 1 && p.rawHk < 1.000001);
  const original = kernel.interval({ upstream: fixture.upstream, downstream: fixture.downstream, regime: 'laminar' });
  original.residual.forEach((v, i) => assert.ok(Math.abs(v - fixture.native.cases[0].nativeResidual[i]) < 1e-12));
  const result = initializeStreamtubeBLStation({ upstream: fixture.upstream, s: fixture.downstream.s, ue: fixture.downstream.ue,
    regime: 'laminar', reynolds: fixture.parameters.reynolds, interval: kernel.interval, properties: kernel.station });
  assert.equal(result.mode, 'inverse'); assert.equal(result.converged, true); assert.equal(result.directReason, 'line search failed');
  assert.ok(Math.abs(result.state.ue - .7234154740869885) < 1e-9);
  assert.ok(Math.abs(kernel.station(result.state).rawHk - result.targetHK) < 1e-10);
  assert.ok(Math.abs(result.state.ue - fixture.downstream.ue) > .4);
  const block = kernel.interval({ upstream: fixture.upstream, downstream: result.state, regime: 'laminar' });
  assert.ok(Math.max(...block.residual.map(Math.abs)) < 1e-10);
  for (const [path, hash] of Object.entries(fixture.inverse.native.sourceHashes))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${path}`, import.meta.url))).digest('hex'), hash);
  const nativeRoot = fixture.inverse.native.cases.find(c => c.name === 'existing inverse initializer');
  for (const key of ['aux', 'theta', 'deltaStar', 'ue']) assert.ok(Math.abs(result.state[key] - nativeRoot.downstream[key]) < 1e-12);
  block.residual.forEach((v, i) => assert.ok(Math.abs(v - nativeRoot.nativeResidual[i]) < 1e-12));
  assert.deepEqual(fixture, before);
});
