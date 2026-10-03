// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { runCoupled } from '../src/viscous/context.js';
import { solveCoupled } from '../src/viscous/tests/solve.js';
import { assembleInterval } from '../src/viscous/tests/interval.js';
import { cft } from '../src/viscous/xfoil/xblsys.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const flow = JSON.parse(readFileSync(new URL('./fixtures/fortran/coupled.json', import.meta.url)));
const kernels = JSON.parse(readFileSync(new URL('./fixtures/fortran/kernels.json', import.meta.url)));
const close = (actual, expected, tolerance, label) => assert.ok(Number.isFinite(actual) && Math.abs(actual - expected) <= tolerance,
  `${label}: ${actual} vs ${expected}, tolerance ${tolerance}`);

test('native fixtures identify the original Fortran and the exact headless harness', () => {
  for (const fixture of [flow, kernels]) for (const [path, hash] of Object.entries(fixture.provenance.sha256)) {
    assert.equal(createHash('sha256').update(readFileSync(path)).digest('hex'), hash, `${path}: regenerate native fixtures after source/harness changes`);
  }
});

for (const fixture of flow.cases) test(`coupled native Fortran parity: ${fixture.name}`, () => {
  const { points, options, expected: f } = fixture; const j = runCoupled(points, options);
  assert.equal(j.converged, f.converged);
  close(j.bl.CL, f.cl, 2e-8, 'CL'); close(j.bl.CM, f.cm, 2e-8, 'Cm');
  close(j.bl.CD, f.cd, 2e-9, 'CD'); close(j.bl.CDF, f.cdf, 1e-8, 'CDf');
  for (let side = 1; side <= 2; side++) close(j.bl.XOCTR[side], f.transition[side - 1], 2e-7, 'transition');
  for (const p of f.cp) {
    close(1 - j.qvis[p.index] ** 2, p.cp, 2e-7, 'viscous Cp');
    close(1 - j.qinv[p.index] ** 2, p.cpInviscid, 2e-9, 'inviscid Cp');
  }
  for (const p of f.bl) {
    assert.equal(j.bl.IPAN[p.station][p.side], p.index);
    for (const [field, key, tolerance] of [['UEDG', 'ue', 2e-7], ['THET', 'theta', 2e-9], ['DSTR', 'deltaStar', 2e-9], ['CTAU', 'ctau', 2e-6]]) {
      close(j.bl[field][p.station][p.side], p[key], tolerance, field);
    }
  }
});

test('laminar, turbulent and wake residuals and Jacobian entries match native Fortran', () => {
  for (const { input, expected } of kernels.intervals) {
    const actual = assembleInterval(input);
    for (const key in expected) expected[key].flat().forEach((value, i) =>
      close(actual[key].flat()[i], value, 2e-12 * (1 + Math.abs(value)), `${input.name}/${key}/${i}`));
  }
});

test('independent finite differences verify the exact parts of the native interval linearization', () => {
  for (const { input: c } of kernels.intervals) {
    const jac = assembleInterval(c);
    for (let side = 0; side < 2; side++) for (let col = 0; col < 5; col++) {
      const index = [c.type === 1 ? 1 : 2, 3, 4, 6, 0][col];
      const h = Math.max(1e-8, Math.abs(c.stations[side][index]) * 1e-5);
      const lo = structuredClone(c); const hi = structuredClone(c);
      lo.stations[side][index] -= h; hi.stations[side][index] += h;
      const l = assembleInterval(lo); const u = assembleInterval(hi);
      for (let row = 0; row < 3; row++) {
        // Original BLDIF omits the explicit UQ_RTA contribution in the turbulent
        // lag row. Preserve native parity; document this approximate Jacobian.
        if (c.type === 2 && row === 0) continue;
        const fd = -(u.residual[row] - l.residual[row]) / (2 * h);
        const a = jac[side ? 'downstream' : 'upstream'][row][col];
        close(fd, a, 1e-6 * (1 + Math.abs(a)), `${c.name}/${side}/${row}/${col}`);
      }
    }
  }
});

test('turbulent skin-friction clamp derivatives match original Fortran and finite differences', () => {
  for (const { input, expected } of kernels.friction) {
    const c = cft(...input); const values = [c.cf, c.cfHk, c.cfRt, c.cfMsq];
    values.forEach((v, i) => close(v, expected[i], 1e-13, 'CFT native'));
    input.forEach((v, i) => {
      const h = 1e-5 * Math.max(Math.abs(v), 1); const lo = [...input]; const hi = [...input];
      lo[i] -= h; hi[i] += h;
      close((cft(...hi).cf - cft(...lo).cf) / (2 * h), values[i + 1], 1e-9, 'CFT derivative');
    });
  }
  assert.equal(Math.abs(cft(1.5, 10, 0).cfRt), 0);
});

test('coupled solution changes pressure, satisfies interaction, resolves transition and wake drag', () => {
  const input = { elements: [{ points: naca4('0012', 160) }], alpha: 4 };
  const r = solveCoupled(input);
  assert.equal(r.status, 'solved'); assert.ok(r.cd > 0.005 && r.cd < 0.01);
  assert.ok(r.diagnostics.couplingResidual < 1e-6);
  assert.ok(Math.max(...r.elements[0].cp.map(p => Math.abs(p.cp - p.cpInviscid))) > 0.1);
  assert.ok(r.history.length > 1 && r.history.at(-1).rmsUpdate < 1e-4);
  const [top, bottom] = r.boundaryLayer.surfaces;
  const firstWake = r.boundaryLayer.wake[0];
  close(firstWake.theta, top.stations.at(-1).theta + bottom.stations.at(-1).theta, 1e-7, 'TE momentum-thickness matching');
  close(firstWake.deltaStar, top.stations.at(-1).deltaStar + bottom.stations.at(-1).deltaStar, 1e-7, 'TE displacement matching');
  const refined = solveCoupled({ ...input, elements: [{ points: naca4('0012', 320) }] });
  close(refined.cl, r.cl, 0.001, 'lift refinement'); close(refined.cd, r.cd, 2e-5, 'drag refinement');
  const forced = solveCoupled({ ...input, trips: [0.05, 0.1] });
  assert.ok(forced.cd > r.cd); close(forced.diagnostics.transition[0], 0.05, 1e-6, 'upper trip');
  close(forced.diagnostics.transition[1], 0.1, 1e-6, 'lower trip');
});

test('reference length and moment origin are explicit in coupled results', () => {
  const points = naca4('2412', 160);
  const r = solveCoupled({ elements: [{ points }], alpha: 4 });
  const shifted = solveCoupled({ elements: [{ points: transform(points, { chord: 2, x: 1.6, y: -0.6 }) }],
    alpha: 4, referenceChord: 2, momentReference: { x: 2.1, y: -0.6 } });
  close(shifted.cl, r.cl, 2e-7, 'CL frame invariance'); close(shifted.cd, r.cd, 2e-8, 'CD frame invariance');
  close(shifted.cm, r.cm, 2e-7, 'Cm frame invariance');
});

test('unconverged and unsupported coupled cases cannot report successful coefficients', () => {
  const input = { elements: [{ points: naca4('0012', 160) }], alpha: 4 };
  const failed = solveCoupled({ ...input, maxIterations: 1 });
  assert.equal(failed.status, 'unconverged'); assert.equal(failed.cl, null); assert.equal(failed.cm, null); assert.equal(failed.cd, null);
  assert.ok(failed.diagnostics.couplingResidual > failed.diagnostics.couplingTolerance);
  for (const change of [{ elements: [...input.elements, ...input.elements] }, { reynolds: NaN }, { mach: 0.3 },
    { trips: [0, 1] }, { ncrit: -1 }, { maxIterations: 0 }, { referenceChord: 0 }]) assert.throws(() => solveCoupled({ ...input, ...change }));
});
