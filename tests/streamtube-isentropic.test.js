import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { createStreamtubeChannel, solveStreamtubeChannel } from '../src/euler/tests/streamtube-channel.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { vortexChannel, vortexErrors, directChannelConservation, nozzleChannel } from './oracles/streamtube.js';
import { directBodyConservation } from './oracles/streamtube-body.js';
import { stagnationPressureError } from './oracles/streamtube-entropy.js';

const max = a => Math.max(...Array.from(a, Math.abs));

test('isentropic cell row equals the independent stagnation-pressure change and rejects a shock entropy jump', () => {
  const gamma = 1.4, lower = [{ x: 0, y: 0 }, { x: .4, y: 0 }, { x: 1, y: 0 }];
  const upper = lower.map(p => ({ ...p, y: .2 }));
  for (const mach of [.2, 2]) {
    const rho = 1, p = 1, q = mach * Math.sqrt(gamma * p / rho);
    const ratio = mach > 1 ? (gamma + 1) * mach ** 2 / ((gamma - 1) * mach ** 2 + 2) : 1.02;
    const c = evaluateStreamtubeCell({ lower, upper, densities: [rho, ratio], massFlow: rho * q * .2,
      stagnationEnthalpy: gamma / (gamma - 1) * p / rho + .5 * q * q, gamma });
    const p0 = c.states.map(s => s.p * (1 + .5 * (gamma - 1) * s.rho * s.q ** 2 / (gamma * s.p)) ** (gamma / (gamma - 1)));
    assert.ok(Math.abs(c.entropyJump + Math.log(p0[1] / p0[0])) < 2e-15);
    assert.ok(Math.abs(c.isentropicResidual - .5 * (c.states[0].p + c.states[1].p) * Math.log(p0[1] / p0[0])) < 5e-15);
    if (mach > 1) { assert.ok(Math.abs(c.streamwiseResidual) < 2e-14); assert.ok(Math.abs(c.isentropicResidual) > .1); }
  }
});

test('isentropic moving-grid vortex preserves stagnation pressure and refines independent pressure, geometry and momentum defects', () => {
  const studies = [];
  for (const [nx, nt] of [[16, 4], [32, 8], [64, 16]]) {
    const input = vortexChannel(nx, nt), system = createStreamtubeChannel({ ...input, streamwiseMode: 'isentropic' });
    const r = solveStreamtubeChannel(system, { tolerance: 1e-11 });
    assert.equal(r.converged, true, r.reason); assert.equal(r.streamwiseMode, 'isentropic');
    const conservation = directChannelConservation(r), pressure = stagnationPressureError(r.sections, { ...input, freestreamMach: .3 });
    assert.ok(pressure.maxRelativeError < 2e-10, JSON.stringify(pressure));
    for (const key of ['maxLocal', 'total', 'external']) for (const k of [0, 3]) assert.ok(Math.abs(conservation[key][k]) < 2e-9);
    assert.ok(max(conservation.internalCancellation) < 2e-9);
    // Isentropy does not impose exact vector momentum in each discrete cell.
    // Measure that truncation error instead of declaring it conserved.
    studies.push({ ...vortexErrors(r, input), momentum: max(conservation.maxLocal.slice(1, 3)) });
  }
  for (let i = 1; i < studies.length; i++) for (const key of ['pressureMax', 'positionMax', 'momentum'])
    assert.ok(studies[i][key] < .3 * studies[i - 1][key], JSON.stringify(studies));
  assert.ok(studies.at(-1).pressureMax < 4e-5); assert.ok(studies.at(-1).positionMax < 2e-6);
});

test('isentropic dense and sparse body/channel solutions agree with independent finite-difference Newton', () => {
  for (const elements of [1, 2]) {
    const input = { ...intrinsicBodyFixture({ elements }), streamwiseMode: 'isentropic' };
    const reference = solveStreamtubeBody(createStreamtubeBodySystem(input), { jacobianBackend: 'finite-difference', tolerance: 1e-11 });
    assert.equal(reference.converged, true, reference.reason);
    for (const linearBackend of ['dense', 'klu']) {
      const system = createStreamtubeBodySystem(input), r = solveStreamtubeBody(system, { linearBackend, tolerance: 1e-11 });
      assert.equal(r.converged, true, r.reason); assert.equal(r.streamwiseMode, 'isentropic');
      assert.ok(max(r.x.map((v, i) => v - reference.x[i])) < 2e-9);
      for (let g = 0; g < r.nodes.length; g++) for (let i = 0; i < r.nodes[g].length; i++) for (let j = 0; j < r.nodes[g][i].length; j++)
        assert.ok(Math.hypot(r.nodes[g][i][j].x - reference.nodes[g][i][j].x, r.nodes[g][i][j].y - reference.nodes[g][i][j].y) < 2e-9);
      const conservation = directBodyConservation(r, input.bodies, system.conditions);
      for (const k of [0, 3]) assert.ok(Math.abs(conservation.balance[k]) < 2e-9);
      assert.ok(max(conservation.cutTraction) < 2e-9);
      const pressure = stagnationPressureError(r.sections.map(row => row.flat()), {
        gamma: input.gamma ?? 1.4, referencePressure: system.conditions.pInf, freestreamMach: input.mach });
      assert.ok(pressure.maxRelativeError < 2e-10, JSON.stringify(pressure));
      assert.ok(Math.abs(pressure.maxRelativeError - r.diagnostics.maxStagnationPressureError) < 3e-15);
      assert.match(r.forceStatus, /Unvalidated/);
    }
  }
  const system = createStreamtubeChannel({ ...nozzleChannel(12), streamwiseMode: 'isentropic' });
  const reference = solveStreamtubeChannel(system, { jacobianBackend: 'finite-difference', tolerance: 1e-12 });
  const sparse = solveStreamtubeChannel(system, { tolerance: 1e-12 });
  assert.equal(reference.converged, true); assert.equal(sparse.converged, true);
  assert.ok(max(sparse.x.map((v, i) => v - reference.x[i])) < 2e-9);
  assert.throws(() => createStreamtubeChannel({ ...nozzleChannel(8), streamwiseMode: 'hybrid' }), /mode/);
  assert.throws(() => createStreamtubeBodySystem({ ...intrinsicBodyFixture(), streamwiseMode: 'hybrid' }), /mode/);
});
