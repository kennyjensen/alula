import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { axset } from '../src/viscous/xfoil/xblsys.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';

const fixtures = JSON.parse(fs.readFileSync(new URL('fixtures/transition-root-failures.json', import.meta.url)));
const compatibility = JSON.parse(fs.readFileSync(new URL('fixtures/transition-root-critical-cycle-compatibility.json', import.meta.url)),
  (_key, value) => value && typeof value === 'object' && value.$negativeZero === true ? -0 : value);
const close = (a, b, tolerance) => assert.ok(Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= tolerance, `${a} versus ${b}`);
const independentRoots = {
  // Independent fraction bracketing: fraction-root-audit/report.json, source
  // hashes retained there. These are not outputs of the safeguarded Newton.
  0: { fraction: .6916370562120933, s: .4022233637256082, amplification: 4.031702957509561 },
  12: { fraction: .7919961886962379, s: .40537314721299994, amplification: 4.020398685597311 },
};

test('all thirteen captured critical cycles resolve the unchanged amplification equation at the original tolerance', () => {
  for (const [index, fixture] of fixtures.cases.entries()) {
    const { parameters, input } = fixture.payload, snapshot = JSON.stringify(input), k = createIntegralKernel(parameters);
    const result = k.transitionCheck(input);
    assert(result.transition && !result.forced);
    const fraction = (result.s - input.upstream.s) / (input.downstream.s - input.upstream.s);
    assert(fraction > 0 && fraction < 1);
    const at = { s: result.s, aux: parameters.ncrit };
    for (const key of ['theta', 'deltaStar', 'ue']) at[key] = input.upstream[key] * (1 - fraction) + input.downstream[key] * fraction;
    const upstream = k.station(input.upstream), transition = k.station(at);
    assert(transition.rawHk > 1 && transition.rho > 0 && transition.viscosity > 0);
    const ax = axset(upstream.hk, input.upstream.theta, upstream.reTheta, input.upstream.aux,
      transition.rawHk, at.theta, transition.reTheta, parameters.ncrit, parameters.ncrit, 0).ax;
    const residual = result.amplification - input.upstream.aux - ax * (input.downstream.s - input.upstream.s);
    close(residual, 0, parameters.transitionTolerance);
    const oracle = independentRoots[index];
    if (oracle) {
      close(fraction, oracle.fraction, 3e-12);
      close(result.s, oracle.s, 3e-12);
      close(result.amplification, oracle.amplification, 3e-12);
    }
    assert.equal(JSON.stringify(input), snapshot);
  }
});

test('the thirteen direct mixed intervals use the same resolved root and finite local equations', () => {
  for (const fixture of fixtures.cases) {
    const { parameters, input } = fixture.payload, k = createIntegralKernel(parameters);
    const root = k.transitionCheck(input);
    const interval = k.interval({ ...input, regime: 'transition' });
    assert.equal(interval.transition.s, root.s);
    assert.equal(interval.transition.amplification, root.amplification);
    assert([...interval.residual, ...interval.upstream.flat(), ...interval.downstream.flat()].every(Number.isFinite));
  }
});

test('all 26 successful native precision cases preserve complete pre-safeguard outputs exactly', () => {
  assert.equal(compatibility.successful.length, 26);
  for (const c of compatibility.successful) {
    const k = createIntegralKernel(c.parameters);
    const actual = c.input.regime === 'transition' ? evaluateTransitionInterval(k, c.input) : k.interval(c.input);
    assert.deepEqual(actual, c.output);
  }
});

test('default accuracy, zero Ncrit and earlier forced trips retain their original local behavior', () => {
  assert.equal(compatibility.exclusions.length, 3);
  for (const c of compatibility.exclusions) {
    const k = createIntegralKernel(c.parameters);
    if (c.error) assert.throws(() => k.transitionCheck(c.input), error => {
      assert.deepEqual({ message: error.message, code: error.code, transitionRootFailure: error.transitionRootFailure }, c.error);
      return true;
    });
    else assert.deepEqual(k.transitionCheck(c.input), c.output);
  }
});

test('the existing endpoint snap still permits a natural root before a terminal material trip', () => {
  const { parameters, input } = fixtures.cases[0].payload, k = createIntegralKernel(parameters);
  const natural = k.transitionCheck(input);
  const snapDistance = 1e-12 * Math.max(1, Math.abs(input.downstream.s));
  const snapped = k.transitionCheck({ ...input, tripS: input.downstream.s - .5 * snapDistance });
  assert.deepEqual(snapped, natural);
});
