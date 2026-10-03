import test from 'node:test';
import assert from 'node:assert/strict';
import { colderStreamtubeStartupMach } from '../src/euler/tests/streamtube-cold-startup.js';
import { evaluateStreamtubeCell, streamtubeCellGeometry } from '../src/euler/streamtube-cell.js';
import { isentropicSectionDensity } from '../src/euler/streamtube-initial-state.js';

test('only a pre-iteration gas-domain rejection permits bounded colder initialization', () => {
  const pressure = { code: 'streamtube-interface-pressure', stage: 'gas-initialization' };
  for (const [mach, next] of [[.2, .1], [.1, .05], [.05, .025], [.03, .025], [.025, null]])
    assert.equal(colderStreamtubeStartupMach(pressure, mach), next);
  const capacity = { code: 'streamtube-sonic-capacity', stage: 'gas-initialization', diagnostics: { stage: 'gas-initialization' } };
  assert.equal(colderStreamtubeStartupMach(capacity, .2), .1);
  capacity.diagnostics.stagnationDensityFallback = { code: 'streamtube-interface-pressure' };
  assert.equal(colderStreamtubeStartupMach(capacity, .2), .1);
  for (const error of [new Error('cancel'), {}, { ...pressure, stage: 'admissibility' },
    { ...capacity, stage: 'boundary-layer-initialization' }, { ...capacity, stage: 'transition-refinement' },
    { ...capacity, stage: undefined },
    { ...capacity, diagnostics: {} }, { ...capacity, diagnostics: { ...capacity.diagnostics,
      stagnationDensityFallback: { reason: 'Crossed, reversed or degenerate quadrilateral polygon.' } } }])
    assert.equal(colderStreamtubeStartupMach(error, .2), null);
});

test('typed static-enthalpy density rejections permit only the existing bounded colder baseline', () => {
  for (const stagnationCode of ['streamtube-interface-pressure', 'streamtube-static-enthalpy'])
    for (const densityCode of [undefined, 'streamtube-interface-pressure', 'streamtube-static-enthalpy']) {
      const error = { code: 'streamtube-sonic-capacity', stage: 'gas-initialization',
        diagnostics: { stage: 'gas-initialization', capacityRatio: 3.881,
          stagnationDensityFallback: { attempted: true, admissible: false, code: stagnationCode },
          ...(densityCode === undefined ? {} : { pressureDomainDensityFallback:
            { attempted: true, admissible: false, code: densityCode } }) } };
      const before = structuredClone(error);
      for (const [mach, next] of [[.2, .1], [.1, .05], [.05, .025], [.03, .025], [.025, null],
        [0, null], [NaN, null], [Infinity, null]])
        assert.equal(colderStreamtubeStartupMach(error, mach), next);
      assert.deepEqual(error, before);
      assert.equal(colderStreamtubeStartupMach({ ...error, stage: 'admissibility' }, .2), null);
      assert.equal(colderStreamtubeStartupMach({ ...error, diagnostics: { ...error.diagnostics, stage: 'euler' } }, .2), null);
    }
  // This change concerns a captured secondary failure. An unrelated primary
  // error does not become recoverable just because its message mentions h.
  assert.equal(colderStreamtubeStartupMach({ code: 'streamtube-static-enthalpy', stage: 'gas-initialization' }, .2), null);
});

test('unknown and geometric density fallback failures remain terminal at every captured level', () => {
  for (const code of ['streamtube-sonic-capacity', 'streamtube-interface-pressure']) {
    const error = { code, stage: 'gas-initialization', diagnostics: { stage: 'gas-initialization',
      stagnationDensityFallback: { code: 'streamtube-static-enthalpy' },
      pressureDomainDensityFallback: { code: 'streamtube-interface-pressure' } } };
    for (const level of ['stagnationDensityFallback', 'pressureDomainDensityFallback'])
      for (const failure of [{}, { code: 'unknown' }, { code: 'streamtube-pressure-domain-transport' },
        { reason: 'Nonpositive streamtube static enthalpy.' },
        { code: 'streamtube-grid-domain', reason: 'Crossed, reversed or degenerate quadrilateral polygon.' }]) {
        const rejected = structuredClone(error);
        rejected.diagnostics[level] = failure;
        const before = structuredClone(rejected);
        assert.equal(colderStreamtubeStartupMach(rejected, .2), null);
        assert.deepEqual(rejected, before);
      }
  }
});

test('captured RAE128 leading-edge cell retains its strict pressure rejection and admits the colder gas', () => {
  // Six physical points from the independent maximum-grid capture, i256/g1/j1.
  const lower = [{ x: .0014404228660905398, y: -.008014894452047748 },
    { x: .001006467704871092, y: -.00797881579269581 }, { x: .0016440374481485288, y: -.00795044954859053 }];
  const upper = [{ x: -.000027869790183397258, y: -.007999044758592784 },
    { x: -.00021069763184409802, y: -.007969000674521337 }, { x: .0001335837723471468, y: -.007942593819575548 }];
  const before = structuredClone({ lower, upper }), massFlow = .00027740425204852526;
  const geometry = streamtubeCellGeometry(lower, upper), gamma = 1.4;
  const h0 = mach => .5 + 1 / ((gamma - 1) * mach * mach);
  const rhoTotal = mach => (1 + .5 * (gamma - 1) * mach * mach) ** (1 / (gamma - 1));
  assert.throws(() => evaluateStreamtubeCell({ lower, upper, massFlow, gamma,
    stagnationEnthalpy: h0(.2), densities: [rhoTotal(.2), rhoTotal(.2)] }), error => {
    assert.equal(error.code, 'streamtube-interface-pressure');
    assert.ok(Math.abs(error.diagnostics.interfacePressure.lower + 6.464876699) < 1e-7);
    return true;
  });
  const densities = geometry.normalAreas.map(area => isentropicSectionDensity({ massFlux: massFlow / area, mach: .1, gamma }));
  const cold = evaluateStreamtubeCell({ lower, upper, massFlow, gamma, stagnationEnthalpy: h0(.1), densities });
  assert.ok(cold.interfacePressure.lower > 42 && cold.interfacePressure.upper > 92);
  assert.ok(cold.states.every(s => s.machSquared < 1 && s.p > 0));
  assert.deepEqual({ lower, upper }, before);
});
