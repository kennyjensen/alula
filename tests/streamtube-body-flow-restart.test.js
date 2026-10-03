import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { initializeStreamtubeBodyFromFlow } from '../src/euler/tests/streamtube-body-flow-restart.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const entropy = (s, c) => Math.log(c.rhoTotal / s.rho) + Math.log(s.enthalpy / c.h0) / (c.gamma - 1);
function fixture() {
  const input = { ...intrinsicBodyFixture({ elements: 1, bodySegments: 4, tubes: 2, mach: .3 }),
    pressureCorrectionFactor: 0, upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const system = createStreamtubeBodySystem(input), initial = system.initial.slice();
  const jump = Math.floor(system.layout.nx / 2);
  for (let i = jump; i < system.layout.nx; i++) for (let g = 0; g < system.layout.tubes.length; g++)
    for (let j = 0; j < system.layout.tubes[g]; j++) initial[system.layout.densityIndex(i, g, j)] = Math.log(1.13);
  // This is an admissible synthetic density/entropy discontinuity, not a
  // converged shock. No Newton iterations or solves are needed to test transfer.
  const flow = system.evaluate(initial);
  assert.ok(Math.abs(entropy(flow.sections[jump][0][0], system.conditions)
    - entropy(flow.sections[jump - 1][0][0], system.conditions)) > .08);
  return { input, system, initial, flow, jump };
}

test('same-Mach physical-flow restart retains the complete density and entropy discontinuity', () => {
  const f = fixture(), before = { input: structuredClone(f.input), initial: f.initial.slice(), chart: f.system.geometryChart(),
    sourceInitial: f.system.initial.slice(), conditions: structuredClone(f.system.conditions) };
  const r = initializeStreamtubeBodyFromFlow(f.input, f.system, { initial: f.initial });
  assert.deepEqual(r.flow.nodes, f.flow.nodes);
  assert.deepEqual(r.flow.sections, f.flow.sections);
  assert.deepEqual(r.flow.residual, f.flow.residual);
  assert.deepEqual(r.flow.captured, f.flow.captured);
  assert.deepEqual(r.flow.stagnation, f.flow.stagnation);
  assert.deepEqual(r.flow.strengths, f.flow.strengths);
  assert.deepEqual(r.initial.subarray(0, f.system.layout.densityCount), f.initial.subarray(0, f.system.layout.densityCount));
  assert.equal(r.diagnostics.maximumEntropyChange, 0);
  assert.equal(r.diagnostics.maximumGeometryChange, 0);
  assert.equal(r.diagnostics.sameMach, true);
  assert.equal(r.diagnostics.converged, false);
  assert.match(r.diagnostics.densityInitialization, /no isentropic inversion/);
  assert.deepEqual(f.input, before.input); assert.deepEqual(f.initial, before.initial);
  assert.deepEqual(f.system.initial, before.sourceInitial); assert.deepEqual(f.system.geometryChart(), before.chart);
  assert.deepEqual(f.system.conditions, before.conditions);
  r.initialEuler.x[0] = 100;
  assert.notEqual(r.initial[0], 100, 'Exported state must not alias the solver state.');
  r.initialEuler.nodes[0][0][0].x += 1;
  assert.deepEqual(f.system.evaluate(f.initial).nodes, f.flow.nodes, 'Exported mesh must not alias the source.');
});

test('Mach continuation copies density and mass, recomputes target gas, and reports the resulting entropy change', () => {
  const f = fixture(), input = { ...f.input, mach: .34 };
  const r = initializeStreamtubeBodyFromFlow(input, f.system, { initial: f.initial });
  assert.deepEqual(r.flow.nodes, f.flow.nodes);
  assert.deepEqual(r.initial.subarray(0, f.system.layout.densityCount), f.initial.subarray(0, f.system.layout.densityCount));
  assert.equal(r.diagnostics.sameMach, false);
  assert.ok(r.diagnostics.maximumEntropyChange > 1e-5);
  assert.match(r.diagnostics.entropyInterpretation, /relative entropy changes/);
  assert.notEqual(r.diagnostics.sourceH0, r.diagnostics.targetH0);
  let expectedMaximum = 0;
  for (let i = 0; i < f.system.layout.nx; i++) for (let g = 0; g < f.system.layout.tubes.length; g++)
    for (let j = 0; j < f.system.layout.tubes[g]; j++) {
      const old = f.flow.sections[i][g][j], next = r.flow.sections[i][g][j], a = f.system.conditions, b = r.system.conditions;
      assert.equal(next.rho, old.rho); assert.equal(next.q, old.q);
      assert.equal(r.flow.allocation.groups[g][j].massFlow, f.flow.allocation.groups[g][j].massFlow);
      assert.ok(next.enthalpy > 0 && next.p > 0);
      assert.equal(next.enthalpy, b.h0 - .5 * old.q * old.q);
      // Independent normalization identity at fixed rho and q.
      const change = Math.log(b.rhoTotal / a.rhoTotal)
        + Math.log((next.enthalpy / b.h0) / (old.enthalpy / a.h0)) / (a.gamma - 1);
      assert.ok(Math.abs(entropy(next, b) - entropy(old, a) - change) < 2e-15);
      expectedMaximum = Math.max(expectedMaximum, Math.abs(change));
    }
  assert.ok(Math.abs(r.diagnostics.maximumEntropyChange - expectedMaximum) < 2e-15);
  assert.notEqual(r.flow.sections[f.jump][0][0].rho, 1);
});

test('moved two-body coordinates, free captured mass, stagnation and multipoles survive rebasing', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .3 }),
    pressureCorrectionFactor: 0, upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const source = createStreamtubeBodySystem(input), initial = source.initial.slice();
  source.layout.positions.forEach(({ column }) => { initial[column] = 1e-5 * Math.sin(column); });
  const g = source.layout.globals;
  for (const col of g.stagnation) if (col !== null) initial[col] = 1e-6;
  for (const col of g.capture) if (col !== null) initial[col] = 1e-5;
  for (const key of ['circulation', 'source', 'doubletX', 'doubletY']) initial[g[key]] = 1e-5;
  const before = source.evaluate(initial), chart = source.geometryChart();
  const r = initializeStreamtubeBodyFromFlow(input, source, { initial });
  assert.deepEqual(r.flow.nodes, before.nodes);
  assert.deepEqual(r.flow.captured, before.captured);
  assert.deepEqual(r.flow.stagnation, before.stagnation);
  assert.deepEqual(r.flow.strengths, before.strengths);
  assert.deepEqual(r.flow.sections, before.sections);
  assert.equal(r.diagnostics.maximumEntropyChange, 0);
  assert.equal(r.diagnostics.maximumGeometryChange, 0);
  r.system.layout.positions.forEach(({ column }) => assert.equal(r.initial[column], 0));
  assert.deepEqual(source.geometryChart(), chart);
  assert.deepEqual(source.evaluate(initial).residual, before.residual);
});

test('matching prescribed displacement preserves both fluid nodes and the underlying solid chart', () => {
  const f = fixture(), input = { ...f.input, displacement: {
    surfaces: f.input.bodies.map(b => ({
      upper: Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => i === 0 ? 0 : .001),
      lower: Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => i === 0 ? 0 : .001),
    })), wakes: f.input.bodies.map(b => Array(f.system.layout.nx - b.trailingIndex).fill(.002)),
  } };
  const source = createStreamtubeBodySystem(input), before = source.evaluate(source.initial);
  const r = initializeStreamtubeBodyFromFlow(input, source);
  assert.deepEqual(r.flow.nodes, before.nodes);
  assert.deepEqual(r.flow.undisplacedNodes, before.undisplacedNodes);
  assert.deepEqual(r.flow.sections, before.sections);
  assert.deepEqual(r.initialEuler.undisplacedNodes, before.undisplacedNodes);
  assert.equal(r.diagnostics.maximumEntropyChange, 0);
});

test('incomplete states, different topology and changed material/reference data are rejected', () => {
  const f = fixture();
  assert.throws(() => initializeStreamtubeBodyFromFlow(f.input, {}), /complete compressible source/);
  assert.throws(() => initializeStreamtubeBodyFromFlow(f.input, f.system, { initial: [1] }), /encoded source state/);
  const incomplete = f.initial.slice(1);
  assert.throws(() => initializeStreamtubeBodyFromFlow(f.input, f.system, { initial: incomplete }), /encoded source state/);
  const mismatch = { ...f.input, weights: f.input.weights.map(w => [...w, 1]) };
  assert.throws(() => initializeStreamtubeBodyFromFlow(mismatch, f.system, { initial: f.initial }), /identical topology/);
  assert.throws(() => initializeStreamtubeBodyFromFlow({ ...f.input, gamma: 1.5 }, f.system, { initial: f.initial }), /changed gamma/);
  assert.throws(() => initializeStreamtubeBodyFromFlow({ ...f.input, alpha: 1 }, f.system, { initial: f.initial }), /changed alpha/);
  const shifted = structuredClone(f.input); shifted.cutPaths[0][0].y += .001;
  assert.throws(() => initializeStreamtubeBodyFromFlow(shifted, f.system, { initial: f.initial }), /identical topology/);
  const displaced = { ...f.input, displacement: {
    surfaces: f.input.bodies.map(b => ({ upper: Array(b.trailingIndex - b.leadingIndex + 1).fill(.001),
      lower: Array(b.trailingIndex - b.leadingIndex + 1).fill(.001) })),
    wakes: f.input.bodies.map(b => Array(f.system.layout.nx - b.trailingIndex).fill(.002)) } };
  assert.throws(() => initializeStreamtubeBodyFromFlow(displaced, f.system, { initial: f.initial }), /displacement/);
  assert.deepEqual(f.system.evaluate(f.initial).residual, f.flow.residual, 'Rejected transfers must not alter the source.');
});
