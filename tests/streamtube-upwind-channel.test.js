// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createUpwindStreamtubeChannel } from '../src/euler/tests/streamtube-upwind-channel.js';
import { solveStreamtubeChannel } from '../src/euler/tests/streamtube-channel.js';
import { sparseDense, sparseProduct } from '../src/numerics/sparse.js';
import { solveLinear, linearResidual } from '../src/numerics/linear.js';
import { channelGas, nozzleChannel } from './oracles/streamtube.js';
import { directUpwindChannelConservation } from './oracles/upwind-streamtube.js';

const close = (a, b, tolerance = 3e-9, label = '') => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${label}: ${a} != ${b}`);
const max = a => Math.max(...Array.from(a, Math.abs));
const input = () => {
  const shape = nozzleChannel(7, 3), gas = channelGas(.75);
  return { ...shape, ...gas, outletPressure: gas.referencePressure,
    upwind: { mucon: 1, mcrit: .6, boundary: { kind: 'unfiltered-first-two' } } };
};
const skewState = system => system.initial.map((v, i) => v + .003 * Math.sin(.7 * (i + 1)));
const fourth = (values, h) => (values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (12 * h);

test('complete upwind-channel sparse and dense Jacobians match every fourth-order residual column', t => {
  const p = input(), system = createUpwindStreamtubeChannel(p), state = skewState(system), { n, nx, nt } = system;
  assert.ok(n < 100);
  assert.equal(n, system.densityCount + system.positionCount + nt);
  const value = system.evaluate(state), dense = system.jacobian(state), sparse = system.jacobian(state, { sparse: true });
  assert.deepEqual(sparseDense(sparse), dense);
  assert.ok(value.sections.every(row => row.every(s => Number.isFinite(s.machSquared))));
  assert.ok(max(value.transportSpeeds.flatMap((row, i) => row.map((q, j) => q - value.sections[i][j].q))) > 1e-5);
  let maximumNormalizedError = 0, comparisons = 0;
  for (let col = 0; col < n; col++) for (const h of [2e-4, 1e-4]) {
    const values = [-2, -1, 1, 2].map(k => { const x = state.slice(); x[col] += k * h; return system.evaluate(x); });
    for (let row = 0; row < n; row++) {
      const fd = fourth(values.map(v => v.residual[row]), h), analytic = dense[row * n + col];
      maximumNormalizedError = Math.max(maximumNormalizedError, Math.abs(fd - analytic) / Math.max(1, Math.abs(fd), Math.abs(analytic)));
      comparisons++; close(analytic, fd, 3e-9, `row ${row}, column ${col}, h ${h}`);
    }
    if (col === system.massIndex(1)) {
      const i = 4, rows = system.densityCount + (i - 1) * (nt - 1);
      const lower = fourth(values.map(v => v.cells[i - 1][1].interfacePressure.lower), h) / p.referencePressure;
      const upper = fourth(values.map(v => v.cells[i - 1][1].interfacePressure.upper), h) / p.referencePressure;
      assert.ok(Math.abs(lower) > .01 && Math.abs(upper) > .01);
      close(dense[rows * n + col], -lower, 3e-9, 'N momentum subtracts the upper passage lower-bank pressure');
      close(dense[(rows + 1) * n + col], upper, 3e-9, 'N momentum adds the lower passage upper-bank pressure');
    }
  }
  // A section two upstream supplies qtilde and its arc geometry. Neither
  // variable belongs to the ordinary three-node cell stencil at i=4.
  const row = 3 * nt + 1, remoteDensity = nt + 1, remoteGeometry = system.positionIndex(1, 1);
  assert.ok(Math.abs(dense[row * n + remoteDensity]) > 1e-7);
  assert.ok(Math.abs(dense[row * n + remoteGeometry]) > 1e-7);
  const unfiltered = createUpwindStreamtubeChannel({ ...p, upwind: { ...p.upwind, mucon: 0 } }).jacobian(state);
  assert.equal(unfiltered[row * n + remoteDensity], 0);
  assert.equal(unfiltered[row * n + remoteGeometry], 0);
  const direction = state.map((_, i) => .3 * Math.cos(i + .2)), h = 1e-4;
  const samples = [-2, -1, 1, 2].map(k => system.residual(state.map((v, i) => v + k * h * direction[i])));
  sparseProduct(sparse, direction).forEach((v, row) => close(v, fourth(samples.map(r => r[row]), h)));
  system.jacobian(state.map(v => .9 * v), { sparse: true });
  assert.deepEqual(sparseDense(sparse), dense, 'Later assembly must not overwrite an existing sparse matrix.');
  t.diagnostic(JSON.stringify({ unknowns: n, comparisons, maximumNormalizedError, streamwiseIntervals: nx }));
});

test('log mass-flow unknowns drive physical mass, reservoir entropy and per-tube back pressure', () => {
  const p = input(), targets = p.massFlows.map((_, j) => p.referencePressure * (1 + .01 * j));
  const system = createUpwindStreamtubeChannel({ ...p, outletPressure: targets }), x = skewState(system);
  const baseCount = system.densityCount + system.positionCount;
  for (let j = 0; j < system.nt; j++) {
    assert.equal(system.massIndex(j), baseCount + j);
    assert.equal(system.initial[system.massIndex(j)], 0);
    x[system.massIndex(j)] = .025 * (j + 1);
  }
  const original = structuredClone(p.massFlows), value = system.evaluate(x);
  for (let j = 0; j < system.nt; j++) {
    assert.equal(value.massFlows[j], p.massFlows[j] * Math.exp(x[system.massIndex(j)]));
    const inlet = value.sections[0][j], outlet = value.sections.at(-1)[j];
    const reservoir = Math.log(inlet.rho / p.stagnationDensity)
      - Math.log(inlet.enthalpy / p.stagnationEnthalpy) / (p.gamma - 1);
    close(value.residual[(system.nx - 1) * system.nt + j], reservoir, 2e-13);
    close(value.residual[baseCount + j], (outlet.p - targets[j]) / p.referencePressure, 2e-13);
  }
  assert.deepEqual(p.massFlows, original);
  assert.ok(max(value.massFlows.map((m, j) => m - p.massFlows[j])) > .001);
  // Inlet slopes plus anchored inlet streamline labels replace outlet slopes.
  const anchorStart = baseCount - (system.nt - 1), dense = system.jacobian(x);
  for (let j = 1; j < system.nt; j++) {
    close(value.residual[anchorStart + j - 1], x[system.positionIndex(0, j)], 2e-13);
    for (let col = 0; col < system.n; col++)
      assert.equal(dense[(anchorStart + j - 1) * system.n + col], col === system.positionIndex(0, j) ? 1 : 0);
  }
  assert.throws(() => createUpwindStreamtubeChannel({ ...p, outletSlopes: [0, 0] }));
});

test('independent channel flux integration separates biased momentum from physical mass and enthalpy', () => {
  const p = input(), system = createUpwindStreamtubeChannel(p), state = skewState(system), value = system.evaluate(state);
  const check = directUpwindChannelConservation(value, p.stagnationEnthalpy, p.gamma);
  close(check.maximumMassMismatch, 0, 2e-12);
  close(check.maximumEnthalpyMismatch, 0, 2e-12);
  assert.ok(check.maximumTransportMomentumDifference > 1e-6);
  for (const key of ['biased', 'physical']) for (const field of ['total', 'external', 'maxLocal'])
    for (const k of [0, 3]) close(check[key][field][k], 0, 2e-12);
  for (const cell of check.local) {
    const r = value.residual[(cell.i - 1) * system.nt + cell.j] * p.referencePressure;
    close(cell.biased[1], r * cell.transverse.y, 2e-12);
    close(cell.biased[2], -r * cell.transverse.x, 2e-12);
  }
  // On this unconverged state, interface-pressure mismatches need not vanish.
  // Shared cross-section transport fluxes still cancel exactly; the physical
  // and biased totals differ solely by their exterior transport correction.
  for (let k = 0; k < 4; k++) close(check.physicalMinusBiased.total[k], check.physicalMinusBiased.external[k], 2e-12);
  const [first, second] = value.transportSpeeds;
  for (let j = 0; j < system.nt; j++) {
    assert.equal(first[j], value.sections[0][j].q);
    assert.equal(second[j], value.sections[1][j].q);
    const centers = value.nodes.map(row => ({ x: .5 * (row[j].x + row[j + 1].x), y: .5 * (row[j].y + row[j + 1].y) }));
    const lengths = centers.slice(1).map((p, i) => Math.hypot(p.x - centers[i].x, p.y - centers[i].y));
    for (let i = 1; i < system.nx; i++) close(value.sectionArcs[j][i] - value.sectionArcs[j][i - 1],
      .5 * (lengths[i - 1] + lengths[i]), 2e-12, 'Spacing follows the bent centerline polyline.');
  }
});

test('a uniform subsonic reservoir/back-pressure channel is an exact root with free tube masses', () => {
  const gas = channelGas(.3), x = [0, .25, .5, .75, 1], height = .2;
  const p = { ...gas, x, lower: x.map(() => 0), upper: x.map(() => height), massFlows: [.04, .07, .09],
    outletPressure: gas.referencePressure, upwind: { mucon: 1, mcrit: .6, boundary: { kind: 'unfiltered-first-two' } } };
  const system = createUpwindStreamtubeChannel(p);
  assert.ok(max(system.residual(system.initial)) < 1e-12);
  const result = solveStreamtubeChannel(system, { tolerance: 1e-12, maxIterations: 1 });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.history.length, 1, 'The exact initial root needs no Newton update.');
  const jacobian = system.jacobian(system.initial), direction = system.initial.map((_, i) => Math.sin(i + .4));
  const rhs = Float64Array.from({ length: system.n }, (_, row) => direction.reduce((sum, v, col) => sum + jacobian[row * system.n + col] * v, 0));
  const solved = solveLinear(jacobian, rhs);
  assert.ok(max(linearResidual(jacobian, solved, rhs)) < 1e-12);
  solved.forEach((v, i) => close(v, direction[i], 3e-9, 'Uniform channel must have no streamline-label nullspace.'));
  const conservation = directUpwindChannelConservation(result, p.stagnationEnthalpy, p.gamma);
  for (const kind of ['biased', 'physical']) for (const key of ['total', 'external', 'maxLocal', 'internalCancellation'])
    assert.ok(max(conservation[kind][key]) < 2e-12);
});

test('supersonic interior sections are admissible while reservoir and outlet remain subsonic', () => {
  const p = input(), system = createUpwindStreamtubeChannel(p), state = system.initial.slice();
  for (let j = 0; j < system.nt; j++) state[3 * system.nt + j] = -.45;
  const value = system.evaluate(state);
  assert.ok(value.sections[3].every(s => s.machSquared > 1));
  assert.ok(value.sections[0].every(s => s.machSquared < 1));
  assert.ok(value.sections.at(-1).every(s => s.machSquared < 1));
  for (const i of [0, system.nx - 1]) {
    const bad = state.slice();
    for (let j = 0; j < system.nt; j++) bad[i * system.nt + j] = -.45;
    assert.throws(() => system.evaluate(bad));
  }
});
