// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { solveUpwindStreamtubeChannelAutomatic } from '../src/euler/tests/streamtube-channel-automatic.js';
import { solveStreamtubeBodyAutomatic } from '../src/euler/tests/streamtube-body-automatic.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { initializeStreamtubeBodyFromFlow } from '../src/euler/tests/streamtube-body-flow-restart.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const controls = { tolerance: 1e-10, directMaxIterations: 0, stageMaxIterations: 15,
  initialFractionStep: .5, maxStages: 6 };
const epsilonP = .001;
const upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
const maximum = a => Math.max(...Array.from(a, Math.abs));
const bodyInput = () => ({ ...intrinsicBodyFixture({ elements: 1, bodySegments: 4, tubes: 2, mach: .2 }),
  streamwiseMode: 'hybrid', hybrid: { epsilonP }, upwind: structuredClone(upwind) });
function channelInput() {
  const gamma = 1.4, stagnationEnthalpy = 5, stagnationDensity = 1;
  const totalPressure = (gamma - 1) / gamma * stagnationDensity * stagnationEnthalpy;
  return { x: [0, .4, 1, 2], lower: [0, 0, 0, 0], upper: [.2, .2, .2, .2],
    massFlows: [.03], gamma, stagnationEnthalpy, stagnationDensity,
    referencePressure: totalPressure, outletPressure: .92 * totalPressure,
    streamwiseMode: 'hybrid', hybrid: { epsilonP }, upwind: structuredClone(upwind) };
}

test('automatic channel continuation retains explicit hybrid rows and epsilon through its case whitelist', () => {
  const input = channelInput(), before = structuredClone(input), stages = [], states = new Map();
  const result = solveUpwindStreamtubeChannelAutomatic(input, { ...controls, linearBackend: 'dense',
    onStage: event => {
      stages.push(structuredClone(event));
      input.hybrid.epsilonP = NaN; input.streamwiseMode = 'momentum';
      event.outletPressure.fill(NaN);
    },
    onState: event => {
      const key = `${event.stage}/${event.continuationFraction}`;
      if (!states.has(key)) states.set(key, []);
      states.get(key).push(event.x.slice());
      event.x.fill(NaN);
    },
  });
  assert.equal(result.system.n, 4);
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.continuation.reachedTarget, true);
  assert.equal(result.streamwiseMode, 'hybrid');
  assert.equal(result.system.conditions.streamwiseMode, 'hybrid');
  assert.deepEqual(result.system.conditions.hybrid, before.hybrid);
  assert.deepEqual(result.system.conditions.outletPressure, [before.outletPressure]);
  assert.ok(maximum(result.residual) < controls.tolerance);
  assert.ok(result.hybrid.flat().every(cell => cell.fraction === 0), 'Uniform subcritical flow selects isentropy.');
  const attempts = result.continuation.attempts;
  assert.equal(attempts.length, 4);
  assert.equal(attempts[0].converged, false);
  assert.ok(attempts.slice(1).every(a => a.converged));
  assert.deepEqual(attempts.map(({ label, fraction, outletPressure }) => ({ label, fraction, outletPressure })), stages);
  for (let i = 2; i < attempts.length; i++) {
    const previous = attempts[i - 1], next = attempts[i];
    assert.deepEqual(states.get(`${next.label}/${next.fraction}`)[0],
      states.get(`${previous.label}/${previous.fraction}`).at(-1), 'A hybrid stage must start from the accepted complete physical state.');
  }
});

test('automatic body Mach continuation keeps the selected hybrid equations and scale after caller mutation', () => {
  const input = bodyInput(), before = structuredClone(input), stages = [];
  const result = solveStreamtubeBodyAutomatic(input, { ...controls, onStage: event => {
    stages.push({ ...event });
    input.hybrid.epsilonP = NaN; input.streamwiseMode = 'momentum';
    event.mach = NaN; event.targetMach = NaN;
  } });
  assert.ok(result.system.layout.n < 150);
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.continuation.reachedTarget, true);
  assert.equal(result.system.conditions.mach, before.mach);
  assert.equal(result.streamwiseMode, 'hybrid');
  assert.deepEqual(result.system.conditions.hybrid, before.hybrid);
  assert.equal(result.diagnostics.hybrid.epsilonP, before.hybrid.epsilonP);
  assert.ok(result.hybridCells.flat(2).every(cell => cell.fraction === 0));
  assert.ok(maximum(result.residual) < controls.tolerance);
  const attempts = result.continuation.attempts;
  assert.equal(attempts.length, 4);
  assert.equal(attempts[0].converged, false);
  assert.ok(attempts.slice(1).every(a => a.converged));
  assert.deepEqual(attempts.map(({ label, mach, fraction }) => ({ label, mach, targetMach: before.mach, fraction })), stages);
  for (const stage of attempts.slice(2)) {
    assert.equal(stage.startup.physicalDensityPreserved, true);
    assert.equal(stage.startup.physicalMassPreserved, true);
    assert.equal(stage.startup.targetMach, stage.mach);
  }
  assert.equal(result.physicalAcceptance, false);
  assert.equal(result.fullSolver, false);
});

test('omitting hybrid mode preserves automatic-body momentum and both wrappers reject a missing hybrid scale', () => {
  const input = intrinsicBodyFixture({ elements: 1, bodySegments: 4, tubes: 2, mach: .2 });
  const result = solveStreamtubeBodyAutomatic(input, { ...controls, stageMaxIterations: 0, maxStages: 2 });
  assert.equal(result.converged, false);
  assert.equal(result.system.conditions.streamwiseMode, 'momentum');
  assert.equal(Object.hasOwn(result.system.conditions, 'hybrid'), false);
  assert.equal(Object.hasOwn(result, 'hybridCells'), false);
  const channel = channelInput(); delete channel.hybrid;
  assert.throws(() => solveUpwindStreamtubeChannelAutomatic(channel, controls), /epsilonP/);
  const body = bodyInput(); delete body.hybrid;
  assert.throws(() => solveStreamtubeBodyAutomatic(body, controls), /hybrid|epsilonP/i);
});

test('physical-flow restart permits Mach changes with the same hybrid scale and rejects changed equations', () => {
  const input = bodyInput(), source = createStreamtubeBodySystem(input);
  const state = source.initial.map((v, i) => v + (i < source.layout.densityCount ? .002 : 1e-6) * Math.sin(i + 1));
  const before = source.evaluate(state), chart = source.geometryChart(), savedState = state.slice();
  const target = { ...input, mach: .3, hybrid: { epsilonP } };
  const prepared = initializeStreamtubeBodyFromFlow(target, source, { initial: state });
  assert.equal(prepared.system.conditions.streamwiseMode, 'hybrid');
  assert.deepEqual(prepared.system.conditions.hybrid, input.hybrid);
  assert.equal(prepared.diagnostics.sourceMach, input.mach);
  assert.equal(prepared.diagnostics.targetMach, target.mach);
  assert.equal(prepared.diagnostics.physicalDensityPreserved, true);
  assert.equal(prepared.diagnostics.physicalMassPreserved, true);
  for (let i = 0; i < source.layout.nx; i++) for (let g = 0; g < source.layout.tubes.length; g++)
    for (let j = 0; j < source.layout.tubes[g]; j++)
      assert.equal(prepared.flow.sections[i][g][j].rho, before.sections[i][g][j].rho);
  assert.ok(prepared.diagnostics.maximumGeometryChange < 1e-12);
  assert.deepEqual(prepared.flow.allocation, before.allocation);
  assert.ok(prepared.flow.residual.every(Number.isFinite));
  assert.throws(() => initializeStreamtubeBodyFromFlow({ ...target, hybrid: { epsilonP: 2 * epsilonP } }, source,
    { initial: state }), /changed hybrid/);
  assert.deepEqual(state, savedState);
  assert.deepEqual(source.geometryChart(), chart);
  assert.deepEqual(source.evaluate(state).residual, before.residual);
});
