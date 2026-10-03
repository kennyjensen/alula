import test from 'node:test';
import assert from 'node:assert/strict';
import { initializeStreamtubeStartup } from '../src/euler/streamtube-startup.js';
import { initializeStreamtubeDensities } from '../src/euler/streamtube-initial-state.js';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const bytes = state => Buffer.from(state.buffer, state.byteOffset, state.byteLength).toString('hex');
const close = (a, b, tolerance = 3e-14) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

// Two straight-centerline finite-volume sections with analytic widths. The
// complete local gas/momentum kernel is real; the small layout avoids meshing,
// panel flow, Jacobians and Newton solves. Nondensity state slots are sentinels.
function channel({ gamma = 1.4, mach = .2, widths = [1.2, 1, .8], massFlow } = {}) {
  const nodes = [widths.map((width, x) => [{ x, y: -.5 * width }, { x, y: .5 * width }])];
  const rhoTotal = (1 + .5 * (gamma - 1) * mach ** 2) ** (1 / (gamma - 1));
  const h0 = 1 / ((gamma - 1) * mach ** 2) + .5, pInf = 1 / (gamma * mach ** 2);
  const criticalTemperature = (1 + .5 * (gamma - 1) * mach ** 2) / (1 + .5 * (gamma - 1));
  const criticalFlux = criticalTemperature ** (1 / (gamma - 1) + .5) / mach;
  massFlow ??= criticalFlux * .9 * 1.005;
  const initial = new Float64Array([.137, .01, -.294, .02, .827]);
  const allocation = { groups: [[{ massFlow }]] };
  const system = { conditions: { flowModel: 'compressible', gamma, mach, rhoTotal, h0, pInf },
    initial, layout: { nx: 2, tubes: [1], densityIndex: i => 2 * i + 1 },
    decode: () => ({ nodes, allocation }), evaluations: 0,
    evaluate(state) {
      this.evaluations++;
      const cell = evaluateStreamtubeCell({ lower: nodes[0].map(row => row[0]), upper: nodes[0].map(row => row[1]),
        densities: [Math.exp(state[1]), Math.exp(state[3])], massFlow, stagnationEnthalpy: h0, gamma });
      if (cell.states.some(s => s.machSquared >= 1)) throw new Error('Local sonic flow is outside this channel verification system.');
      const sections = cell.states.map(s => [[s]]), residual = new Float64Array([cell.isentropicResidual / pInf]);
      const entropy = cell.states.map(s => Math.log(s.p / pInf) - gamma * Math.log(s.rho));
      const maxStagnationPressureError = Math.max(...entropy.map(e => Math.abs(Math.expm1(-e / (gamma - 1)))));
      return { nodes, allocation, sections, residual, cells: [[[cell]]], diagnostics: {
        residual: Math.max(...residual.map(Math.abs)), residualByFamily: { streamwise: Math.abs(residual[0]) },
        maxMach: Math.sqrt(Math.max(...cell.states.map(s => s.machSquared))),
        maxEntropyJump: Math.abs(cell.entropyJump), maxStagnationPressureError } };
    } };
  return { system, widths, massFlow, criticalFlux };
}

test('successful exact body initialization remains byte-identical and evaluates only once', () => {
  const body = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, alpha: 2, tubes: 5, tubeGrowth: 3 }));
  const inputBytes = bytes(body.initial), conditions = structuredClone(body.conditions);
  const decoded = structuredClone(body.decode(body.initial));
  const exact = initializeStreamtubeDensities(body, body.initial);
  let calls = 0;
  const system = { ...body, evaluate: state => { calls++; return body.evaluate(state); } };
  const result = initializeStreamtubeStartup(system, body.initial);
  assert.equal(bytes(result.initial), bytes(exact)); assert.equal(calls, 1);
  assert.equal(result.diagnostics.method, 'isentropic');
  assert.equal(result.diagnostics.isentropicFailure, undefined);
  assert.ok(result.diagnostics.entropy.maximumAbsoluteDeparture < 3e-14);
  assert.equal(bytes(body.initial), inputBytes); assert.deepEqual(body.conditions, conditions);
  assert.deepEqual(body.decode(body.initial), decoded);
  assert.deepEqual(result.flow.residual, body.evaluate(exact).residual);
});

test('stagnation-density seed preserves mass and enthalpy while explicitly retaining entropy residuals', () => {
  for (const gamma of [1.2, 1.4, 1.67]) for (const mach of [.2, .4, .6]) {
    const { system, widths, massFlow } = channel({ gamma, mach });
    const inputBytes = bytes(system.initial), conditions = structuredClone(system.conditions), geometry = structuredClone(system.decode());
    const result = initializeStreamtubeStartup(system, system.initial);
    assert.equal(result.diagnostics.method, 'stagnation-density'); assert.equal(system.evaluations, 1);
    const d = result.diagnostics.isentropicFailure;
    assert.equal(d.code, 'streamtube-sonic-capacity'); close(d.diagnostics.capacityRatio, 1.005);
    assert.deepEqual(d.diagnostics.section, { i: 1, group: 0, tube: 0 });
    const entropy = [];
    for (let i = 0; i < 2; i++) {
      const s = result.flow.sections[i][0][0], area = .5 * (widths[i] + widths[i + 1]);
      const q = massFlow / (conditions.rhoTotal * area), h = conditions.h0 - .5 * q ** 2;
      const p = (gamma - 1) / gamma * conditions.rhoTotal * h;
      close(s.rho, conditions.rhoTotal); close(s.q, q); close(s.p, p);
      close(s.rho * s.q * area, massFlow); close(s.enthalpy + .5 * s.q ** 2, conditions.h0);
      assert.ok(s.p > 0 && s.enthalpy > 0 && s.machSquared < 1);
      entropy.push(Math.log(p / conditions.pInf) - gamma * Math.log(conditions.rhoTotal));
    }
    close(result.diagnostics.entropy.maximumAbsoluteDeparture, Math.max(...entropy.map(Math.abs)));
    close(result.diagnostics.maxEntropyJump, Math.abs(entropy[1] - entropy[0]) / (gamma - 1));
    assert.ok(result.diagnostics.entropy.maximumAbsoluteDeparture > .01);
    assert.ok(result.diagnostics.residual > .001);
    assert.equal(result.diagnostics.targetEquationsUnchanged, true);
    for (const i of [0, 2, 4]) assert.equal(result.initial[i], system.initial[i]);
    assert.equal(bytes(system.initial), inputBytes); assert.deepEqual(system.conditions, conditions);
    assert.deepEqual(system.decode(), geometry);
  }
});

test('a truly sonic stagnation-density seed keeps the original capacity error and rejection reason', () => {
  const { system } = channel({ widths: [1, 1, 1], massFlow: 6 });
  const inputBytes = bytes(system.initial); let original;
  try { initializeStreamtubeDensities(system, system.initial); } catch (error) { original = error; }
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), error => {
    assert.equal(error.code, original.code); assert.equal(error.message, original.message);
    const { stagnationDensityFallback, fallbackAttempted, ...strictDiagnostics } = error.diagnostics;
    assert.equal(fallbackAttempted, true);
    assert.deepEqual(strictDiagnostics, original.diagnostics);
    assert.deepEqual(stagnationDensityFallback, { attempted: true, method: 'stagnation-density', admissible: false,
      reason: 'Local sonic flow is outside this channel verification system.' });
    return true;
  });
  assert.equal(system.evaluations, 1); assert.equal(bytes(system.initial), inputBytes);
});

test('geometry failure is never masked by an alternative density guess', () => {
  const { system } = channel({ massFlow: 10 });
  system.decode().nodes[0][2][1] = { x: .8, y: .2 };
  const inputBytes = bytes(system.initial);
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), /Initial grid i=1, group=0, tube=0: Folded/);
  assert.equal(system.evaluations, 0); assert.equal(bytes(system.initial), inputBytes);
});

test('non-pressure full evaluation errors after successful exact inversion never invoke fallback', () => {
  const { system } = channel({ massFlow: .5 });
  const failure = Object.assign(new Error('Evaluation rejected the exact seed.'), { code: 'streamtube-sonic-capacity' });
  system.evaluate = () => { system.evaluations++; throw failure; };
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), error => error === failure);
  assert.equal(system.evaluations, 1); assert.equal(failure.diagnostics, undefined);
});

test('fallback full evaluation failure is retained without changing the strict capacity diagnosis', () => {
  const { system } = channel();
  system.evaluate = () => { system.evaluations++; throw new Error('Nonpositive interface pressure.'); };
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), error => {
    assert.equal(error.code, 'streamtube-sonic-capacity'); close(error.diagnostics.capacityRatio, 1.005);
    assert.deepEqual(error.diagnostics.stagnationDensityFallback,
      { attempted: true, method: 'stagnation-density', admissible: false, reason: 'Nonpositive interface pressure.' });
    return true;
  });
  assert.equal(system.evaluations, 1);
});

test('typed fallback pressure details are detached while the original capacity error and state remain unchanged', () => {
  const { system } = channel(), before = bytes(system.initial), geometry = structuredClone(system.decode());
  let original;
  try { initializeStreamtubeDensities(system, system.initial); } catch (error) { original = error; }
  assert.equal(original.code, 'streamtube-sonic-capacity');
  const pressure = Object.assign(new Error('Nonpositive or nonfinite streamline interface pressure.'), {
    code: 'streamtube-interface-pressure', diagnostics: { interfacePressure: { lower: -.2, upper: Infinity },
      pressureSum: Infinity, pressureDifference: NaN, normalInertia: -4, cell: { i: 124, group: 0, tube: 2 } },
  });
  const expected = structuredClone(pressure.diagnostics);
  system.evaluate = () => { system.evaluations++; throw pressure; };
  let caught;
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), error => {
    caught = error;
    assert.equal(error.code, original.code); assert.equal(error.message, original.message);
    const { fallbackAttempted, stagnationDensityFallback, pressureDomainDensityFallback, ...strict } = error.diagnostics;
    assert.equal(fallbackAttempted, true); assert.deepEqual(strict, original.diagnostics);
    assert.deepEqual(stagnationDensityFallback, { attempted: true, method: 'stagnation-density', admissible: false,
      reason: pressure.message, code: pressure.code, diagnostics: expected });
    assert.equal(pressureDomainDensityFallback.method, 'pressure-domain-density');
    assert.equal(pressureDomainDensityFallback.admissible, false);
    assert.deepEqual(pressureDomainDensityFallback.diagnostics, expected);
    assert.ok(pressureDomainDensityFallback.candidate.minimumCertifiedPressure > 0);
    return true;
  });
  pressure.diagnostics.cell.i = 999; pressure.diagnostics.interfacePressure.lower = 1;
  assert.deepEqual(caught.diagnostics.stagnationDensityFallback.diagnostics, expected);
  assert.deepEqual(caught.diagnostics.pressureDomainDensityFallback.diagnostics, expected);
  assert.equal(system.evaluations, 2); assert.equal(bytes(system.initial), before);
  assert.deepEqual(system.decode(), geometry);
});

test('missing physical stagnation density rejects fallback before flow evaluation', () => {
  const { system } = channel(); delete system.conditions.rhoTotal;
  assert.throws(() => initializeStreamtubeStartup(system, system.initial), error => {
    assert.equal(error.code, 'streamtube-sonic-capacity');
    assert.equal(error.diagnostics.stagnationDensityFallback.reason, 'Invalid physical stagnation density.');
    return true;
  });
  assert.equal(system.evaluations, 0);
});
