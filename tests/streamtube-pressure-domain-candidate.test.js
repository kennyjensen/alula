// SPDX-License-Identifier: GPL-2.0-or-later
// Small isolated policy/domain fixtures, no real mesh/flow/J/LU solve.
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { StreamtubePressureDomainTransportError as Typed,
  evaluateStreamtubePressureDomainCandidate as evaluate } from '../src/euler/streamtube-pressure-domain-initial-state.js';
const currentStartup = new URL('../src/euler/streamtube-startup.js', import.meta.url);
const originalStartup = new URL('../docs/solver-reliability/gui-defaults-current/30p30n/pressure-domain-evaluated-draft/before-promotion/streamtube-startup.js', import.meta.url);
const fixture = () => {
  const conditions = { flowModel: 'compressible', upwind: { mucon: 1, mcrit: .99 }, hybrid: { ismom: 4 }, gamma: 1.4, pInf: 1 };
  const initial = new Float64Array([Math.log(2), Math.log(2), .125]), state = new Float64Array([0, 0, .125]);
  const diagnostics = { density: 2, densityUnknowns: 2, minimumCertifiedPressure: .3,
    transportCheck: { biasedSections: 1, maxTransportDeparture: .02 } };
  const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }], [{ x: 2, y: 0 }, { x: 2, y: 1 }]]];
  const flow = { nodes: structuredClone(nodes), sections: [0, 1].map(() => [[{ p: 2, rho: 2, enthalpy: 3, machSquared: .5 }]]),
    cells: [[[{ interfacePressure: { lower: .5, upper: .25 } }]]], residual: new Float64Array([1, 2, 3]),
    diagnostics: { residual: 3, residualByFamily: { streamwise: 1 }, maxMach: Math.sqrt(.5), maxEntropyJump: .1, maxStagnationPressureError: .2 } };
  let calls = 0;
  const chart = [{ column: 2, offset: { x: 0, y: 0 }, normal: { x: 0, y: 1 } }];
  const system = { conditions, decode: () => ({ nodes }), geometryChart: () => structuredClone(chart), layout: { nx: 2, tubes: [1], densityIndex: i => i }, evaluate: x => {
    calls++; assert.deepEqual(x, initial); return flow; } };
  return { initial, state, diagnostics, flow, system, chart, originNodes: nodes, get calls() { return calls; }, error: new Typed(initial, diagnostics, system, nodes) };
};
test('actual biased candidate admission records failed certificate separately and preserves controls/state', () => {
  const f = fixture(), before = structuredClone({ state: f.state, conditions: f.system.conditions, candidate: f.error.candidate });
  const result = evaluate(f.system, f.state, f.error);
  assert.equal(f.calls, 1); assert.equal(result.flow, f.flow); assert.notEqual(result.initial, f.error.candidate.initial);
  assert.equal(result.diagnostics.pressureCertificateApplicable, false);
  assert.equal(result.diagnostics.minimumCertifiedPressure, undefined);
  assert.equal(result.diagnostics.actualEvaluation.minInterfacePressure, .25);
  assert.equal(result.diagnostics.certificateInapplicable.diagnostics.minimumCertifiedPressure, .3);
  assert.deepEqual({ state: f.state, conditions: f.system.conditions, candidate: f.error.candidate }, before);
});
test('wrong model, changed selected controls and changed nondensity state are rejected before evaluation', () => {
  for (const change of [f => { f.system.conditions.flowModel = 'incompressible'; },
    f => { f.system.conditions.upwind.mucon = 0; }, f => { f.state[2] += .01; },
    f => { f.error.candidate.initial[0] = NaN; }, f => { f.error.candidate.initial[0] = Math.log(3); },
    f => { f.error.candidate.initial.fill(Math.log(3), 0, 2); f.error.diagnostics.density = 3; }]) {
    const f = fixture(); change(f); assert.throws(() => evaluate(f.system, f.state, f.error)); assert.equal(f.calls, 0);
  }
});
test('a different solver, adopted physical geometry or refreshed chart cannot reuse the candidate', () => {
  const noChart = fixture(); delete noChart.system.geometryChart;
  const rejected = new Typed(noChart.initial, noChart.diagnostics, noChart.system, noChart.originNodes);
  assert.equal(rejected.code, 'streamtube-pressure-domain-transport');
  assert.throws(() => evaluate(noChart.system, noChart.state, rejected), /originating/);
  assert.equal(noChart.calls, 0);
  const other = fixture(), clone = { ...other.system };
  assert.throws(() => evaluate(clone, other.state, other.error), /originating system or geometry chart/);
  assert.equal(other.calls, 0);
  for (const change of [f => { f.chart[0].normal.x = .2; },
    f => { f.chart[0].offset.y = .01; }, f => { f.originNodes[0][1][0].x += .01; }]) {
    const f = fixture(); change(f); assert.throws(() => evaluate(f.system, f.state, f.error), /originating/);
    assert.equal(f.calls, 0);
  }
});
test('a code-only or unrelated/observer exception does not authorize evaluation', () => {
  for (const error of [new Error('observer cancelled'), Object.assign(new Error('not the generated candidate'),
    { code: 'streamtube-pressure-domain-transport' }), Object.assign(new Error('negative pressure'),
    { code: 'streamtube-interface-pressure' })]) {
    const f = fixture(); assert.throws(() => evaluate(f.system, f.state, error), e => e === error); assert.equal(f.calls, 0);
  }
});
test('negative actual pressure, enthalpy, nonfinite residuals and a concave grid remain rejected', () => {
  for (const change of [f => { f.flow.cells[0][0][0].interfacePressure.upper = -.1; },
    f => { f.flow.sections[0][0][0].enthalpy = -.1; }, f => { f.flow.sections[0][0][0].p = -.1; },
    f => { f.flow.residual[0] = NaN; }, f => { f.flow.nodes[0][1][0].y = 1.5; }]) {
    const f = fixture(); change(f); assert.throws(() => evaluate(f.system, f.state, f.error)); assert.equal(f.calls, 1);
  }
});
test('an exception from actual evaluation propagates unchanged without retry or caller mutation', () => {
  const f = fixture(), error = Object.assign(new Error('observer cancelled'), { code: 'observer-cancellation' });
  const before = f.state.slice(); let calls = 0; f.system.evaluate = () => { calls++; throw error; };
  assert.throws(() => evaluate(f.system, f.state, f.error), e => e === error);
  assert.equal(calls, 1); assert.deepEqual(f.state, before);
});
let serial = 0;
async function startup(draft, pressureOutcome) {
  const f = fixture(), original = Object.assign(new Error('bad side pressure'), { code: 'streamtube-interface-pressure', diagnostics: { original: true } });
  let calls = 0, pressureCalls = 0;
  f.system.evaluate = x => { calls++; if (calls === 1) throw original; assert.deepEqual(x, f.initial); return f.flow; };
  const key = `__startup${serial++}`;
  globalThis[key] = { initializeStreamtubeDensities: () => f.state.slice(),
    initializeStreamtubePressureDomain: () => { pressureCalls++; if (pressureOutcome === 'typed') throw f.error;
      if (pressureOutcome instanceof Error) throw pressureOutcome; return { initial: f.initial.slice(), diagnostics: { untouched: true } }; },
    StreamtubePressureDomainTransportError: Typed, evaluateStreamtubePressureDomainCandidate: evaluate };
  let source = fs.readFileSync(draft ? currentStartup : originalStartup, 'utf8');
  source = source.replace(/^import \{ ([^}]+) \} from .*;$/gm, (_, names) => `const { ${names} } = globalThis[${JSON.stringify(key)}];`);
  const module = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`); delete globalThis[key];
  return { f, original, run: () => module.initializeStreamtubeStartup(f.system, f.state), counts: () => ({ calls, pressureCalls }) };
}
test('healthy unfiltered-certificate startup output and evaluation counts equal the original exactly', async () => {
  const old = await startup(false), draft = await startup(true);
  assert.deepEqual(draft.run(), old.run()); assert.deepEqual(draft.counts(), old.counts());
});
test('startup routes only the typed certificate failure to the same candidate once', async () => {
  const h = await startup(true, 'typed'), result = h.run();
  assert.deepEqual(result.initial, h.f.initial); assert.equal(result.diagnostics.method, 'pressure-domain-density-evaluated');
  assert.equal(result.diagnostics.pressureCertificateApplicable, false);
  assert.equal(result.diagnostics.isentropicFailure.code, 'streamtube-interface-pressure');
  assert.deepEqual(h.counts(), { calls: 2, pressureCalls: 1 });
});
test('startup preserves other failures without another gas evaluation', async () => {
  for (const failure of [Object.assign(new Error('cancelled'), { code: 'observer-cancellation' }),
    Object.assign(new Error('inapplicable spoof'), { code: 'streamtube-pressure-domain-transport' }),
    new Error('invalid geometry')]) {
    const h = await startup(true, failure);
    assert.throws(h.run, e => e === h.original);
    assert.deepEqual(h.counts(), { calls: 1, pressureCalls: 1 });
    assert.equal(h.original.diagnostics.pressureDomainDensityFallback.reason, failure.message);
  }
});
