import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { initializeStreamtubePressureDomain } from '../src/euler/streamtube-pressure-domain-initial-state.js';
import { initializeStreamtubeStartup } from '../src/euler/streamtube-startup.js';
import { initializeStreamtubeDensities, isentropicSonicMassFlux } from '../src/euler/streamtube-initial-state.js';
import { evaluateStreamtubeCell, streamtubeCellGeometry } from '../src/euler/streamtube-cell.js';
import { prepareStreamtubeTransportChain } from '../src/euler/streamtube-transport-chain.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/nlr7301-pressure-cell.json', import.meta.url)));
const close = (a, b, tolerance = 5e-13) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const bytes = state => Buffer.from(state.buffer, state.byteOffset, state.byteLength).toString('hex');

// Real finite-volume gas/momentum kernels on supplied physical banks. This
// keeps the failure geometry and equation checks without any meshing or LU.
function localSystem({ lower = fixture.lower, upper = fixture.upper, massFlow = fixture.massFlow,
  gamma = fixture.gamma, mach = fixture.mach, pressureCorrectionFactor = .1, upwind } = {}) {
  const nx = lower.length - 1, nodes = [lower.map((p, i) => [{ ...p }, { ...upper[i] }])];
  const allocation = { groups: [[{ massFlow }]] };
  const h0 = 1 / ((gamma - 1) * mach ** 2) + .5, pInf = 1 / (gamma * mach ** 2);
  const rhoTotal = (1 + .5 * (gamma - 1) * mach ** 2) ** (1 / (gamma - 1));
  const initial = Float64Array.from({ length: 2 * nx + 1 }, (_, i) => (i % 2 ? 0 : .123 * (i + 1)));
  const conditions = { flowModel: 'compressible', h0, pInf, rhoTotal, gamma, mach, pressureCorrectionFactor,
    geometryDomain: 'convex', ...(upwind ? { upwind: structuredClone(upwind) } : {}) };
  const system = { conditions, initial, layout: { nx, tubes: [1], densityIndex: i => 2 * i + 1 },
    evaluations: 0, decode: () => ({ nodes, allocation }),
    evaluate(state) {
      this.evaluations++;
      const densities = Array.from({ length: nx }, (_, i) => Math.exp(state[2 * i + 1]));
      const chain = upwind ? prepareStreamtubeTransportChain({ lower, upper, densities, massFlow,
        stagnationEnthalpy: h0, gamma, upwind }) : null;
      const cells = [], sections = [];
      for (let i = 1; i < nx; i++) {
        const cell = evaluateStreamtubeCell({ lower: lower.slice(i - 1, i + 2), upper: upper.slice(i - 1, i + 2),
          densities: densities.slice(i - 1, i + 1), massFlow, stagnationEnthalpy: h0, gamma, pressureCorrectionFactor,
          ...(chain ? { transportSpeeds: chain.transportSpeeds.slice(i - 1, i + 1) } : {}) });
        assert.ok(cell.states.every(s => s.machSquared < 1));
        cells.push([[cell]]);
        if (i === 1) sections.push([[cell.states[0]]]);
        sections.push([[cell.states[1]]]);
      }
      const residual = Float64Array.from(cells, row => row[0][0].isentropicResidual / pInf);
      return { nodes, allocation, cells, sections, residual, diagnostics: {
        residual: Math.max(...residual.map(Math.abs)), residualByFamily: { streamwise: Math.max(...residual.map(Math.abs)) },
        maxMach: Math.sqrt(Math.max(...sections.map(row => row[0][0].machSquared))),
        maxEntropyJump: Math.max(...cells.map(row => Math.abs(row[0][0].entropyJump))) } };
    } };
  return system;
}

test('a saved convex NLR cell has positive section gas but negative side pressure at both ordinary seeds', () => {
  const system = localSystem(), exact = initializeStreamtubeDensities(system, system.initial);
  assert.throws(() => system.evaluate(exact), error => {
    assert.equal(error.code, 'streamtube-interface-pressure');
    close(error.diagnostics.interfacePressure.lower, fixture.strictMinimumPressure);
    return true;
  });
  const stagnation = system.initial.slice();
  for (let i = 0; i < 2; i++) stagnation[2 * i + 1] = Math.log(system.conditions.rhoTotal);
  assert.throws(() => system.evaluate(stagnation), error => {
    assert.equal(error.code, 'streamtube-interface-pressure');
    close(error.diagnostics.interfacePressure.lower, fixture.stagnationMinimumPressure);
    return true;
  });
  const before = bytes(system.initial), controls = structuredClone(system.conditions), geometry = structuredClone(system.decode());
  system.evaluations = 0;
  const result = initializeStreamtubeStartup(system, system.initial);
  assert.equal(system.evaluations, 2);
  assert.equal(result.diagnostics.method, 'pressure-domain-density');
  assert.equal(result.diagnostics.isentropicFailure.code, 'streamtube-interface-pressure');
  close(result.diagnostics.isentropicFailure.diagnostics.interfacePressure.lower, fixture.strictMinimumPressure);
  assert.ok(result.diagnostics.minimumCertifiedPressure > 0);
  assert.ok(result.flow.cells[0][0][0].interfacePressure.lower >= result.diagnostics.minimumCertifiedPressure);
  assert.ok(result.diagnostics.entropy.maximumAbsoluteDeparture > .001);
  assert.equal(result.diagnostics.initialGuessOnly, true); assert.equal(result.diagnostics.targetEquationsUnchanged, true);
  assert.equal(result.diagnostics.converged, undefined);
  assert.equal(bytes(system.initial), before); assert.deepEqual(system.decode(), geometry); assert.deepEqual(system.conditions, controls);
  for (let i = 0; i < system.initial.length; i += 2) assert.equal(result.initial[i], system.initial[i]);
  for (const row of result.flow.sections) {
    const section = row[0][0];
    close(section.enthalpy + .5 * section.q ** 2, controls.h0);
    close(section.p, (controls.gamma - 1) / controls.gamma * section.rho * section.enthalpy);
  }
});

test('uniform-density certificate bounds the full signed-curvature kernel across gas conditions', () => {
  for (const gamma of [1.2, 1.4, 1.67]) for (const mach of [.185, .4, .7]) for (const bend of [-.2, .2]) {
    const lower = [{ x: 0, y: -.6 }, { x: 1, y: -.4 + bend }, { x: 2, y: -.55 }];
    const upper = [{ x: 0, y: .6 }, { x: 1.1, y: .45 + bend }, { x: 2.2, y: .5 }];
    const system = localSystem({ lower, upper, massFlow: 3, gamma, mach, pressureCorrectionFactor: .8 });
    const seed = initializeStreamtubePressureDomain(system, system.initial), { density, bounds } = seed.diagnostics;
    assert.equal(system.evaluations, 0);
    const cell = system.evaluate(seed.initial).cells[0][0][0], { A, B, C, D, K } = bounds.pressure;
    const correction = cell.pressureCorrection;
    assert.ok(D > 0 && K !== 0);
    assert.ok(Math.abs(correction) <= D / density);
    close(cell.states[0].p + cell.states[1].p, A * density - B / density);
    close(cell.interfacePressure.upper - cell.interfacePressure.lower, C / density + K * correction);
    const minimum = Math.min(cell.interfacePressure.lower, cell.interfacePressure.upper);
    assert.ok(minimum >= seed.diagnostics.minimumCertifiedPressure - 1e-13);
    assert.ok(cell.states.every(s => s.enthalpy > 0 && s.machSquared < 1));
    close(Math.max(...cell.states.map(s => s.machSquared)), seed.diagnostics.transportCheck.maximumMachSquared);
    assert.ok(bounds.sonic.squaredDensity > bounds.thermal.squaredDensity);
    close(density ** 2, seed.diagnostics.boundarySquaredDensity / .95);
    const geometry = streamtubeCellGeometry(lower, upper);
    cell.states.forEach((s, i) => close(s.rho * s.q * geometry.normalAreas[i], 3));
  }
});

test('pressure-controlled density scales with captured mass without changing geometry or H0', () => {
  const a = localSystem(), b = localSystem({ massFlow: 2 * fixture.massFlow });
  const sa = initializeStreamtubePressureDomain(a, a.initial), sb = initializeStreamtubePressureDomain(b, b.initial);
  assert.ok(sa.diagnostics.bounds.pressure.squaredDensity > a.conditions.rhoTotal ** 2);
  close(sb.diagnostics.bounds.pressure.squaredDensity, 4 * sa.diagnostics.bounds.pressure.squaredDensity);
  close(sb.diagnostics.density, 2 * sa.diagnostics.density);
  const ca = a.evaluate(sa.initial).cells[0][0][0], cb = b.evaluate(sb.initial).cells[0][0][0];
  ca.states.forEach((s, i) => { close(cb.states[i].q, s.q); close(cb.states[i].p, 2 * s.p); });
  close(cb.interfacePressure.lower, 2 * ca.interfacePressure.lower);
  close(cb.interfacePressure.upper, 2 * ca.interfacePressure.upper);
  assert.equal(a.conditions.h0, b.conditions.h0);
});

test('a capacity failure followed by real stagnation-seed pressure failure retains both diagnoses on success', () => {
  const system = localSystem({ massFlow: .4 }), before = bytes(system.initial);
  const result = initializeStreamtubeStartup(system, system.initial);
  assert.equal(result.diagnostics.method, 'pressure-domain-density');
  assert.equal(result.diagnostics.isentropicFailure.code, 'streamtube-sonic-capacity');
  assert.ok(result.diagnostics.isentropicFailure.diagnostics.capacityRatio > 1);
  const previous = result.diagnostics.stagnationDensityFallback;
  assert.equal(previous.code, 'streamtube-interface-pressure');
  assert.equal(previous.admissible, false); assert.ok(previous.diagnostics.interfacePressure.lower < 0);
  assert.ok(result.diagnostics.minimumCertifiedPressure > 0);
  assert.ok(result.flow.cells[0][0][0].interfacePressure.lower > 0);
  assert.equal(system.evaluations, 2); assert.equal(bytes(system.initial), before);
});

test('288.1% excess capacity and negative stagnation-seed enthalpy admit a certified fixed-grid density guess', () => {
  const mach = .2, gamma = 1.4, massFlow = 3.881 * isentropicSonicMassFlux({ mach, gamma });
  const lower = [0, 1, 2, 3].map(x => ({ x, y: 0 })), upper = lower.map(p => ({ ...p, y: 1 }));
  const upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
  const system = localSystem({ lower, upper, massFlow, mach, gamma, upwind });
  const before = bytes(system.initial), geometry = structuredClone(system.decode()), controls = structuredClone(system.conditions);
  const result = initializeStreamtubeStartup(system, system.initial);
  assert.equal(result.diagnostics.method, 'pressure-domain-density');
  assert.equal(result.diagnostics.isentropicFailure.code, 'streamtube-sonic-capacity');
  close(result.diagnostics.isentropicFailure.diagnostics.capacityRatio, 3.881);
  const thermal = result.diagnostics.stagnationDensityFallback;
  assert.equal(thermal.code, 'streamtube-static-enthalpy'); assert.equal(thermal.admissible, false);
  assert.deepEqual(thermal.diagnostics.section, { i: 0, bankStations: [0, 1] });
  assert.ok(thermal.diagnostics.enthalpy < 0);
  assert.ok(result.diagnostics.density > controls.rhoTotal);
  assert.equal(result.diagnostics.transportCheck.biasedSections, 0);
  assert.equal(result.diagnostics.targetEquationsUnchanged, true);
  assert.equal(result.diagnostics.initialGuessOnly, true);
  assert.ok(result.diagnostics.entropy.maximumAbsoluteDeparture > .1);
  assert.equal(result.diagnostics.converged, undefined);
  assert.equal(system.evaluations, 2);
  assert.equal(bytes(system.initial), before); assert.deepEqual(system.decode(), geometry); assert.deepEqual(system.conditions, controls);
  for (let i = 0; i < system.initial.length; i += 2) assert.equal(result.initial[i], system.initial[i]);
  for (const row of result.flow.sections) {
    const s = row[0][0];
    close(s.rho * s.q, massFlow); close(s.enthalpy + .5 * s.q ** 2, controls.h0);
    close(s.p, (gamma - 1) / gamma * s.rho * s.enthalpy);
    assert.ok(s.enthalpy > 0 && s.p > 0 && s.machSquared < 1);
  }
});

test('thermal recovery still evaluates and retains rejection from the complete selected equations', () => {
  const lower = [0, 1, 2].map(x => ({ x, y: 0 })), upper = lower.map(p => ({ ...p, y: 1 }));
  const system = localSystem({ lower, upper, massFlow: 12, mach: .2 });
  const evaluate = system.evaluate.bind(system), before = bytes(system.initial);
  const failure = Object.assign(new Error('Selected remote boundary rejected this density.'), { code: 'remote-boundary-state' });
  system.evaluate = state => { const flow = evaluate(state); assert.ok(flow.sections[0][0][0].p > 0); throw failure; };
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), error => {
    assert.equal(error.code, 'streamtube-sonic-capacity');
    assert.equal(error.diagnostics.stagnationDensityFallback.code, 'streamtube-static-enthalpy');
    const rejected = error.diagnostics.pressureDomainDensityFallback;
    assert.equal(rejected.code, failure.code); assert.equal(rejected.admissible, false);
    assert.ok(rejected.candidate.minimumCertifiedPressure > 0);
    return true;
  });
  assert.equal(system.evaluations, 2); assert.equal(bytes(system.initial), before);
});

test('a locally subsonic certificate never disables active selected speed upwinding', () => {
  const lower = [0, 1, 2, 3].map((x, i) => ({ x, y: -.5 * [1, .9, .7, .8][i] }));
  const upper = lower.map(p => ({ x: p.x, y: -p.y }));
  const upwind = { mucon: 1, mcrit: 0, boundary: { kind: 'unfiltered-first-two' } };
  const system = localSystem({ lower, upper, massFlow: 2, mach: .4, upwind });
  const before = bytes(system.initial), controls = structuredClone(system.conditions);
  assert.throws(() => initializeStreamtubePressureDomain(system, system.initial), error => {
    assert.equal(error.code, 'streamtube-pressure-domain-transport');
    assert.ok(error.diagnostics.transportCheck.biasedSections > 0);
    assert.ok(error.diagnostics.transportCheck.maximumMachSquared < 1);
    assert.ok(error.diagnostics.minimumCertifiedPressure > 0);
    return true;
  });
  assert.equal(system.evaluations, 0); assert.equal(bytes(system.initial), before); assert.deepEqual(system.conditions, controls);
  const physical = localSystem({ lower, upper, massFlow: 2, mach: .4, upwind: { ...upwind, mcrit: 1 } });
  const seed = initializeStreamtubePressureDomain(physical, physical.initial);
  assert.equal(seed.diagnostics.transportCheck.biasedSections, 0);
  assert.ok(physical.evaluate(seed.initial).diagnostics.maxMach < 1);
});

test('full selected-equation failure preserves the original pressure error and candidate evidence', () => {
  const system = localSystem(), evaluate = system.evaluate.bind(system), initial = bytes(system.initial);
  let original;
  const rejected = Object.assign(new Error('Selected boundary row rejects this gas state.'), {
    code: 'remote-boundary-state', diagnostics: { endpoint: { side: 'exit', mach: 1.1 } } });
  system.evaluate = state => {
    if (!original) {
      try { return evaluate(state); } catch (error) { original = error; throw error; }
    }
    evaluate(state); // The candidate's real cell passes; a separate full-flow condition fails.
    throw rejected;
  };
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), error => {
    assert.equal(error, original); assert.equal(error.code, 'streamtube-interface-pressure');
    close(error.diagnostics.interfacePressure.lower, fixture.strictMinimumPressure);
    const fallback = error.diagnostics.pressureDomainDensityFallback;
    assert.equal(fallback.code, rejected.code); assert.equal(fallback.admissible, false);
    assert.deepEqual(fallback.diagnostics, rejected.diagnostics);
    assert.ok(fallback.candidate.minimumCertifiedPressure > 0);
    rejected.diagnostics.endpoint.side = 'mutated';
    assert.equal(fallback.diagnostics.endpoint.side, 'exit');
    return true;
  });
  assert.equal(system.evaluations, 2); assert.equal(bytes(system.initial), initial);
});

test('invalid geometry and gas parameters are not repaired by the analytical helper', () => {
  const invalid = localSystem(); invalid.decode().nodes[0][2][1] = { x: -1, y: -1 };
  assert.throws(() => initializeStreamtubePressureDomain(invalid, invalid.initial), /Folded|degenerate/);
  assert.equal(invalid.evaluations, 0);
  for (const [key, value] of [['gamma', 1], ['h0', 0], ['rhoTotal', NaN], ['pressureCorrectionFactor', -1]]) {
    const system = localSystem(); system.conditions[key] = value;
    assert.throws(() => initializeStreamtubePressureDomain(system, system.initial), /initialization conditions/);
    assert.equal(system.evaluations, 0);
  }
});
