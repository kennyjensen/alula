// SPDX-License-Identifier: GPL-2.0-or-later
// A nonlinear shock discriminator, independent of the RAE startup topology.
// The channel uses the production cell, upwind and entropy/momentum kernels.
import test from 'node:test';
import assert from 'node:assert/strict';
import { transonicNozzleReference } from './oracles/transonic-nozzle.js';
import { createUpwindStreamtubeChannel } from '../src/euler/tests/streamtube-upwind-channel.js';
import { solveStreamtubeChannel } from '../src/euler/tests/streamtube-channel.js';

const maximum = values => values.reduce((m, v) => Math.max(m, Math.abs(v)), 0);

function nozzle(nx, nt, shear = 0) {
  const exact = transonicNozzleReference({ shockMach: 1.5, areaRatio: 2, throatArea: .02 });
  const x = Array.from({ length: nx + 1 }, (_, i) =>
    exact.domain[0] + (exact.domain[1] - exact.domain[0]) * i / nx);
  const input = {
    x, lower: x.map(z => shear * z - .5 * exact.area(z)),
    upper: x.map(z => shear * z + .5 * exact.area(z)),
    inletSlopes: Array(nt - 1).fill(shear),
    massFlows: Array(nt).fill(exact.chokedMass / nt / Math.hypot(1, shear)),
    stagnationEnthalpy: exact.stagnationEnthalpy, stagnationDensity: exact.stagnationDensity,
    referenceDensity: 1, referencePressure: exact.totalPressure,
    outletPressure: exact.backPressure, streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } },
  };
  const system = createUpwindStreamtubeChannel(input);
  const analytic = system.initial.slice();
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++)
    analytic[i * nt + j] = Math.log(exact.stateAt(.5 * (x[i] + x[i + 1])).rho);
  return { exact, system, analytic, x };
}

for (const [nx, nt, shear] of [[32, 3, 0], [64, 3, 0], [128, 3, 0], [64, 12, 1]]) {
  test(`hybrid Euler captures a nozzle shock from uniform gas: ${nx} intervals, ${nt} tubes, shear ${shear}`, t => {
    const { exact, system, analytic, x } = nozzle(nx, nt, shear);
    const initial = system.initial.slice();
    assert.ok(system.evaluate(initial).sections.every(row => row.every(s => s.machSquared < 1)),
      'the cold seed must not contain a pre-positioned supersonic region or shock');
    const cold = solveStreamtubeChannel(system, { maxIterations: 30, tolerance: 1e-10 });
    const warm = solveStreamtubeChannel(system, { initial: analytic, maxIterations: 15, tolerance: 1e-10 });
    for (const result of [cold, warm]) {
      assert.equal(result.converged, true, result.reason);
      assert.ok(maximum(result.residual) <= 1e-10);
      assert.ok(result.sections.every(row => row.every(s => s.rho > 0 && s.p > 0 && Number.isFinite(s.q))));
      assert.ok(result.sections.some(row => row.some(s => s.machSquared > 1.5)), 'a resolved supersonic region is required');
      assert.ok(result.sections.at(-1).every(s => s.machSquared < 1));
    }
    assert.ok(maximum(cold.x.map((v, i) => v - warm.x[i])) < 1e-7, 'both starts must reach the same discrete shock solution');
    assert.deepEqual(system.initial, initial, 'solving must not mutate the next cold start');
    const mach = cold.sections.map(row => Math.sqrt(row[0].machSquared));
    let shockIndex = 1;
    for (let i = 2; i < nx; i++)
      if (mach[i - 1] - mach[i] > mach[shockIndex - 1] - mach[shockIndex]) shockIndex = i;
    const dx = x[1] - x[0], shockError = Math.abs(x[shockIndex] - exact.shockX);
    assert.ok(shockError < 1.5 * dx, `shock location error ${shockError} exceeds cell spacing ${dx}`);
    const massError = Math.abs(cold.massFlows.reduce((a, b) => a + b, 0)
      * Math.hypot(1, shear) / exact.chokedMass - 1);
    // Quasi-1D is an independent slender-nozzle reference, not an exact
    // moving-grid solution. This bound shrinks with streamwise refinement.
    assert.ok(massError < 6 / nx ** 2, `relative choked-mass error ${massError}`);
    t.diagnostic(JSON.stringify({ coldUpdates: cold.history.length - 1, warmUpdates: warm.history.length - 1,
      residual: maximum(cold.residual), shockError, massError }));
  });
}
