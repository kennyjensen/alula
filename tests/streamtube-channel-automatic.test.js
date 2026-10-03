// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { solveUpwindStreamtubeChannelAutomatic } from '../src/euler/tests/streamtube-channel-automatic.js';

const max = a => Math.max(...Array.from(a, Math.abs));
const close = (a, b, tolerance = 2e-9) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
function uniform(referenceMass = .03) {
  const gamma = 1.4, stagnationEnthalpy = 5, stagnationDensity = 1;
  const totalPressure = (gamma - 1) / gamma * stagnationDensity * stagnationEnthalpy;
  return { x: [0, .4, 1, 2], lower: [0, 0, 0, 0], upper: [.2, .2, .2, .2],
    massFlows: [referenceMass], gamma, stagnationEnthalpy, stagnationDensity,
    referencePressure: totalPressure, outletPressure: .92 * totalPressure,
    upwind: { mucon: 1, mcrit: .6, boundary: { kind: 'unfiltered-first-two' } } };
}
const controls = { tolerance: 1e-11, linearBackend: 'dense', directMaxIterations: 12, stageMaxIterations: 20 };

test('automatic geometry/reservoir startup reaches the requested uniform root independently of reference mass', () => {
  const results = [];
  for (const referenceMass of [.03, 3]) {
    const input = uniform(referenceMass), before = structuredClone(input);
    const result = solveUpwindStreamtubeChannelAutomatic(input, controls);
    assert.ok(result.system.n <= 10);
    assert.equal(result.converged, true, result.reason);
    assert.equal(result.continuation.reachedTarget, true);
    assert.equal(result.automaticInitialization.initialGuessOnly, true);
    assert.equal(result.continuation.analyticShockInformationUsed, false);
    assert.deepEqual(result.system.conditions.outletPressure, [input.outletPressure]);
    assert.deepEqual(result.continuation.targetPressure, [input.outletPressure]);
    assert.ok(max(result.residual) < controls.tolerance);
    const rho = input.stagnationDensity * .92 ** (1 / input.gamma);
    const h = input.stagnationEnthalpy * .92 ** ((input.gamma - 1) / input.gamma);
    const q = Math.sqrt(2 * (input.stagnationEnthalpy - h));
    for (const row of result.sections) {
      close(row[0].rho, rho); close(row[0].q, q); close(row[0].p, input.outletPressure);
    }
    close(result.massFlows[0], .2 * rho * q);
    assert.deepEqual(input, before);
    results.push(result);
  }
  close(results[0].massFlows[0], results[1].massFlows[0]);
  results[0].sections.forEach((row, i) => {
    for (const key of ['rho', 'q', 'p', 'enthalpy', 'machSquared']) close(row[0][key], results[1].sections[i][0][key]);
  });
});

test('zero update budgets terminate without accepting a failed startup or intermediate pressure', () => {
  const input = uniform();
  const zero = solveUpwindStreamtubeChannelAutomatic(input,
    { ...controls, directMaxIterations: 0, stageMaxIterations: 0, maxStages: 2 });
  assert.equal(zero.converged, false);
  assert.equal(zero.continuation.reachedTarget, false);
  assert.equal(zero.continuation.attempts.length, 2);
  assert.ok(zero.continuation.attempts.every(a => !a.converged && a.iterations === 0));
  assert.match(zero.reason, /startup failed/i);

  const capped = solveUpwindStreamtubeChannelAutomatic(input,
    { ...controls, directMaxIterations: 0, maxStages: 2 });
  assert.equal(capped.converged, false);
  assert.equal(capped.continuation.reachedTarget, false);
  assert.equal(capped.continuation.attempts.length, 2);
  assert.equal(capped.continuation.attempts.at(-1).converged, true);
  assert.ok(max(capped.residual) < controls.tolerance,
    'An intermediate converged residual must not be called a converged target solution.');
  assert.notDeepEqual(capped.system.conditions.outletPressure, [input.outletPressure]);
  assert.deepEqual(capped.continuation.targetPressure, [input.outletPressure]);
  assert.match(capped.reason, /stage limit/i);
});

test('pressure continuation copies accepted states and isolates observers and caller input edits', () => {
  const input = uniform(), before = structuredClone(input), stages = [], snapshots = new Map();
  let oracleCalls = 0;
  input.unusedOracle = () => { oracleCalls++; throw new Error('An analytic oracle must not be used.'); };
  const result = solveUpwindStreamtubeChannelAutomatic(input, { ...controls, directMaxIterations: 0,
    initialFractionStep: .5, maxStages: 8,
    onStage: event => {
      stages.push(structuredClone(event));
      event.outletPressure.fill(NaN); event.label = 'observer mutation'; event.fraction = -1;
      input.massFlows.fill(1e6); input.x[1] = NaN; input.stagnationEnthalpy = 0;
      input.outletPressure = 100; input.upwind.mucon = NaN;
    },
    onIteration: event => { event.outletPressure.fill(NaN); event.residual = NaN; },
    onState: event => {
      const key = `${event.stage}:${event.continuationFraction}`;
      if (!snapshots.has(key)) snapshots.set(key, []);
      snapshots.get(key).push(structuredClone(event));
      event.x.fill(NaN); event.iteration.residual = NaN; event.outletPressure.fill(NaN);
    },
  });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.continuation.reachedTarget, true);
  assert.deepEqual(result.system.conditions.outletPressure, [before.outletPressure]);
  assert.deepEqual(result.continuation.targetPressure, [before.outletPressure]);
  assert.deepEqual(result.system.conditions.massFlows, before.massFlows);
  assert.deepEqual(result.system.conditions.stagnationEnthalpy, [before.stagnationEnthalpy]);
  assert.deepEqual(result.nodes.map(row => row[0].x), before.x);
  assert.equal(oracleCalls, 0);
  const attempts = result.continuation.attempts;
  assert.equal(attempts[0].converged, false);
  assert.equal(attempts[1].converged, true);
  assert.equal(attempts.at(-1).fraction, 1);
  assert.deepEqual(attempts.at(-1).outletPressure, [before.outletPressure]);
  assert.deepEqual(attempts.map(a => ({ label: a.label, fraction: a.fraction, outletPressure: a.outletPressure })), stages);
  assert.ok(attempts.every(a => a.outletPressure.every(Number.isFinite)
    && a.history.every(h => Number.isFinite(h.residual))));
  for (let i = 2; i < attempts.length; i++) {
    assert.equal(attempts[i].converged, true);
    const previous = attempts[i - 1], next = attempts[i];
    const last = snapshots.get(`${previous.label}:${previous.fraction}`).at(-1);
    const first = snapshots.get(`${next.label}:${next.fraction}`)[0];
    assert.equal(first.iteration.iteration, 0);
    assert.deepEqual(first.x, last.x, 'Continuation must preserve all encoded physical densities and free masses.');
  }
  assert.ok(result.x.every(Number.isFinite));
  assert.ok(max(result.residual) < controls.tolerance);
});

test('positive reservoir flow rejects back pressure at or above total pressure', () => {
  const input = uniform(), p0 = (input.gamma - 1) / input.gamma * input.stagnationDensity * input.stagnationEnthalpy;
  for (const outletPressure of [p0, 1.01 * p0])
    assert.throws(() => solveUpwindStreamtubeChannelAutomatic({ ...input, outletPressure }, controls), /below inlet total pressure/);
});
