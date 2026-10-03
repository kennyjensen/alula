import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeEdgeState, streamtubeEdgePressure } from '../src/euler/streamtube-edge-velocity.js';
import { isentropicState } from '../src/potential/isentropic.js';
import { createIntegralKernel } from '../src/viscous/integral.js';

const close = (a, b, tolerance = 2e-10) => assert.ok(Math.abs(a - b) < tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
test('pressure-matched BL edges recover independent perfect-gas states and their pressure derivatives', () => {
  for (const mach of [.02, .2, .6]) for (const gamma of [1.4, 1.67]) {
    const conditions = { flowModel: 'compressible', mach, gamma, pInf: 1 / (gamma * mach ** 2) };
    const native = createIntegralKernel({ reynolds: 1e6, mach, gamma });
    for (const ue of [.2, .8, 1, 1.2]) {
      const gas = isentropicState(ue, 0, conditions), pressure = conditions.pInf + .5 * gas.cp;
      const actual = streamtubeEdgeState(pressure, conditions);
      const direct = streamtubeEdgePressure(ue, conditions);
      close(direct.pressure, pressure, 2e-14);
      close(direct.derivative, -gas.rho * ue, 2e-14);
      close(actual.ue, ue); close(actual.rho, gas.rho); close(actual.machSquared, gas.machSquared);
      const bl = native.station({ s: .2, ue: actual.ue, theta: .0006, deltaStar: .0016, aux: 0 });
      close(bl.rho, actual.rho); close(bl.machSquared, actual.machSquared);
      const h = 1e-6, fd = (streamtubeEdgeState(pressure + h, conditions).ue - streamtubeEdgeState(pressure - h, conditions).ue) / (2 * h);
      close(actual.velocityPressureDerivative, fd, 2e-7);
    }
    const gm = gamma - 1, qSonic = Math.sqrt((1 + .5 * gm * mach ** 2) / (.5 * (gamma + 1) * mach ** 2));
    const p = q => conditions.pInf * (1 + .5 * gm * mach ** 2 * (1 - q * q)) ** (gamma / gm);
    assert.throws(() => streamtubeEdgeState(p(1.01 * qSonic), conditions), /sonic/);
    assert.throws(() => streamtubeEdgePressure(1.01 * qSonic, conditions), /sonic/);
    assert.throws(() => streamtubeEdgeState(p(0) * (1 + 1e-10), conditions), /stagnating/);
    assert.throws(() => streamtubeEdgeState(0, conditions), /Nonpositive/);
  }
});

test('incompressible edge matching uses the Bernoulli pressure gauge and rejects unresolved stagnation', () => {
  const conditions = { flowModel: 'incompressible' };
  for (const ue of [.1, .5, 1, 2]) {
    const state = streamtubeEdgeState(-.5 * ue * ue, conditions);
    close(state.ue, ue); close(state.rho, 1); close(state.velocityPressureDerivative, -1 / ue);
    assert.deepEqual(streamtubeEdgePressure(ue, conditions), { pressure: -.5 * ue * ue, derivative: -ue });
  }
  for (const pressure of [0, 1, NaN]) assert.throws(() => streamtubeEdgeState(pressure, conditions));
});

test('direct pressure residual remains finite above Euler stagnation pressure', () => {
  const conditions = { flowModel: 'compressible', mach: .01, gamma: 1.4, pInf: 1 / (1.4 * .01 ** 2) };
  const impossibleInversion = conditions.pInf + 1;
  assert.throws(() => streamtubeEdgeState(impossibleInversion, conditions), /stagnating/);
  const matching = streamtubeEdgePressure(.9, conditions);
  assert.ok(Number.isFinite(matching.pressure - impossibleInversion));
  const h = 1e-4, fd = (streamtubeEdgePressure(.9 + h, conditions).pressure
    - streamtubeEdgePressure(.9 - h, conditions).pressure) / (2 * h);
  close(matching.derivative, fd, 2e-8);
  for (const ue of [0, -1, NaN, Infinity]) assert.throws(() => streamtubeEdgePressure(ue, conditions));
});
