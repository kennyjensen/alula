import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeChannel, solveStreamtubeChannel } from '../src/euler/tests/streamtube-channel.js';
import { channelGas, nozzleChannel, vortexChannel, vortexErrors, directChannelConservation } from './oracles/streamtube.js';

const max = a => Math.max(...Array.from(a, Math.abs));
const close = (a, b, tol = 2e-10) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);
const straight = (nt = 3, slope = 0) => {
  const x = [0, .1, .3, .55, .8, 1.2, 1.5, 1.8, 2], height = .2;
  return { ...channelGas(), x, lower: x.map(x => slope * x), upper: x.map(x => slope * x + height),
    massFlows: Array.from({ length: nt }, (_, j) => height * 2 * (j + 1) / (nt * (nt + 1)) / Math.hypot(1, slope)),
    inletSlopes: Array(nt - 1).fill(slope), outletSlopes: Array(nt - 1).fill(slope) };
};

test('intrinsic channel closes a square system and recovers uniform flow from perturbed density and grid guesses', () => {
  for (const nt of [1, 3]) for (const slope of [0, .25]) {
    const input = straight(nt, slope), system = createStreamtubeChannel(input);
    assert.equal(system.n, 8 * nt + 9 * (nt - 1));
    assert.equal(system.residual(system.initial).length, system.n);
    assert.ok(max(system.residual(system.initial)) < 1e-14);
    const initial = system.initial.map((_, i) => (i < system.densityCount ? .005 : .003) * Math.sin(1.7 * (i + 1)));
    assert.ok(max(system.residual(initial)) > 1e-3);
    const r = solveStreamtubeChannel(system, { initial, tolerance: 1e-11 });
    assert.equal(r.converged, true, r.reason);
    assert.ok(max(r.residual) < 1e-11);
    assert.ok(max(r.x) < 1e-9);
    r.sections.flat().forEach(s => { close(s.rho, 1); close(s.q, 1); close(s.p, input.referencePressure); });
  }
});

test('intrinsic channel rejects folds, invalid inputs and sonic states and reports an iteration-limited solve', () => {
  const input = straight(), system = createStreamtubeChannel(input);
  const folded = system.initial.slice(); folded[system.positionIndex(3, 1)] = 2;
  assert.equal(system.admissible(folded), false);
  assert.throws(() => system.evaluate(folded), /Folded|area/);
  assert.throws(() => createStreamtubeChannel({ ...input, x: [0, 0, 1] }), /Invalid channel/);
  assert.throws(() => createStreamtubeChannel({ ...input, massFlows: [0] }), /positive mass/);
  assert.throws(() => createStreamtubeChannel({ ...input, stagnationEnthalpy: -1 }), /stagnation enthalpy/);
  assert.throws(() => createStreamtubeChannel({ ...input, inletSlopes: [NaN, 0] }), /slopes/);
  const sonic = createStreamtubeChannel({ ...input, massFlows: input.massFlows.map(m => 3.5 * m) });
  assert.equal(sonic.admissible(sonic.initial), false);
  assert.throws(() => sonic.evaluate(sonic.initial), /subsonic/);
  const initial = system.initial.map((_, i) => .002 * Math.sin(i + 1));
  const limited = solveStreamtubeChannel(system, { initial, maxIterations: 0 });
  assert.equal(limited.converged, false); assert.equal(limited.reason, 'iteration limit');
  assert.ok(max(limited.residual) > 1e-4);
});

test('both flow/grid coupling blocks are active, including inlet thermodynamics and end geometry', () => {
  const system = createStreamtubeChannel(nozzleChannel(8)), state = system.initial.map((_, i) => .001 * Math.sin(i + 1));
  const h = 2e-6, n = system.n, nd = system.densityCount, j = new Float64Array(n * n);
  let flowFromGrid = 0, gridFromFlow = 0;
  for (let col = 0; col < n; col++) {
    const p = state.slice(), m = state.slice(); p[col] += h; m[col] -= h;
    const rp = system.residual(p), rm = system.residual(m);
    for (let row = 0; row < n; row++) {
      const v = (rp[row] - rm[row]) / (2 * h); j[row * n + col] = v;
      if (row < nd && col >= nd) flowFromGrid = Math.max(flowFromGrid, Math.abs(v));
      if (row >= nd && col < nd) gridFromFlow = Math.max(gridFromFlow, Math.abs(v));
    }
  }
  assert.ok(flowFromGrid > .01); assert.ok(gridFromFlow > .01);
  // A separate whole-state perturbation catches missing end rows or a wrong
  // layout; this verifies assembly, not an analytic-Jacobian claim.
  const direction = state.map((_, i) => Math.cos(.7 * (i + 1))), step = 1e-5;
  const rp = system.residual(state.map((v, i) => v + step * direction[i]));
  const rm = system.residual(state.map((v, i) => v - step * direction[i]));
  for (let row = 0; row < n; row++) {
    let product = 0; for (let col = 0; col < n; col++) product += j[row * n + col] * direction[col];
    close(product, (rp[row] - rm[row]) / (2 * step), 2e-7);
  }
});

test('the assembled curved channel conserves independently integrated vector fluxes and shared interfaces', () => {
  const input = nozzleChannel(12), system = createStreamtubeChannel(input), r = solveStreamtubeChannel(system, { tolerance: 1e-12 });
  assert.equal(r.converged, true, r.reason);
  assert.ok(max(r.x.subarray(system.densityCount)) > .001, 'Streamlines must move with the solved flow.');
  const check = directChannelConservation(r, input.gamma);
  for (const key of ['maxLocal', 'total', 'external', 'internalCancellation'])
    assert.ok(max(check[key]) < 2e-11, JSON.stringify(check));
});

test('intrinsic nozzle refinement approaches the independent slender isentropic comparison', () => {
  const errors = [];
  for (const nx of [8, 16, 32]) {
    const input = nozzleChannel(nx), r = solveStreamtubeChannel(createStreamtubeChannel(input), { tolerance: 1e-11 });
    assert.equal(r.converged, true, r.reason);
    const dp = r.sections.flatMap((row, i) => row.map(s => Math.abs(s.p - input.exact((input.x[i] + input.x[i + 1]) / 2).p)));
    errors.push(dp.reduce((s, v) => s + v, 0) / dp.length);
  }
  assert.ok(errors[1] < .4 * errors[0], String(errors));
  assert.ok(errors[2] < .5 * errors[1], String(errors));
  assert.ok(errors[2] < .0006, String(errors));
  // A quasi-1D nozzle is a slender-channel comparison, not an exact 2D
  // solution. The vortex below supplies the exact continuum refinement gate.
});

test('simultaneous intrinsic flow/grid refinement converges to the exact 2D compressible vortex', () => {
  const errors = [];
  for (const [nx, nt] of [[8, 2], [16, 4], [32, 8]]) {
    const input = vortexChannel(nx, nt), r = solveStreamtubeChannel(createStreamtubeChannel(input), { tolerance: 1e-11 });
    assert.equal(r.converged, true, r.reason);
    errors.push(vortexErrors(r, input));
    const flux = directChannelConservation(r, input.gamma);
    assert.ok(max(flux.external) < 2e-9, JSON.stringify(flux));
  }
  for (let i = 1; i < errors.length; i++) for (const key of ['pressureMax', 'positionMax'])
    assert.ok(errors[i][key] < .3 * errors[i - 1][key], JSON.stringify(errors));
  assert.ok(errors.at(-1).pressureMax < .00015, JSON.stringify(errors));
  assert.ok(errors.at(-1).positionMax < 8e-6, JSON.stringify(errors));
});
