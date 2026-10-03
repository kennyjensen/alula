// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile, readdir } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { runCoupled } from '../src/viscous/context.js';
import { solveCoupled } from '../src/viscous/tests/solve.js';
import { assembleInterval } from '../src/viscous/tests/interval.js';
import { naca4 } from '../src/geometry/airfoil.js';

const flow = JSON.parse(await readFile('tests/fixtures/fortran/coupled.json', 'utf8'));
const kernels = JSON.parse(await readFile('tests/fixtures/fortran/kernels.json', 'utf8'));
const report = { generatedAt: new Date().toISOString(), runtime: process.version, status: 'passed',
  scope: 'Executed original-Fortran parity for the single-element incompressible coupled XFOIL model; not experimental or multielement MSES validation.',
  nativeProvenance: flow.provenance, cases: [], intervals: [], refinement: [], sha256: {} };
for (const [path, hash] of Object.entries(flow.provenance.sha256)) {
  assert.equal(createHash('sha256').update(await readFile(path)).digest('hex'), hash, path);
}
for (const f of flow.cases) {
  const j = runCoupled(f.points, f.options); const n = f.expected;
  const errors = { cl: Math.abs(j.bl.CL - n.cl), cm: Math.abs(j.bl.CM - n.cm), cd: Math.abs(j.bl.CD - n.cd),
    cdf: Math.abs(j.bl.CDF - n.cdf), cpMax: Math.max(...n.cp.map(p => Math.abs(1 - j.qvis[p.index] ** 2 - p.cp))),
    transitionMax: Math.max(...n.transition.map((v, side) => Math.abs(v - j.bl.XOCTR[side + 1]))) };
  for (const [key, name] of [['UEDG', 'ue'], ['THET', 'theta'], ['DSTR', 'deltaStar'], ['CTAU', 'ctau']]) {
    errors[`${name}Max`] = Math.max(...n.bl.map(s => Math.abs(j.bl[key][s.station][s.side] - s[name])));
  }
  assert.equal(j.converged, n.converged); assert.ok(errors.cl < 2e-8 && errors.cm < 2e-8 && errors.cd < 2e-9 && errors.cdf < 1e-8);
  assert.ok(errors.cpMax < 2e-7 && errors.transitionMax < 2e-7 && errors.ueMax < 2e-7);
  assert.ok(errors.thetaMax < 2e-9 && errors.deltaStarMax < 2e-9 && errors.ctauMax < 2e-6);
  report.cases.push({ name: f.name, panels: f.points.length - 1, conditions: f.options, converged: j.converged,
    iterations: j.history.length, native: { cl: n.cl, cm: n.cm, cd: n.cd, transition: n.transition },
    javascript: { cl: j.bl.CL, cm: j.bl.CM, cd: j.bl.CD }, errors });
}
for (const { input, expected } of kernels.intervals) {
  const actual = assembleInterval(input);
  const nativeScaledMax = Math.max(...Object.keys(expected).flatMap(key => expected[key].flat().map((v, i) => Math.abs(actual[key].flat()[i] - v) / (1 + Math.abs(v)))));
  assert.ok(nativeScaledMax < 2e-12);
  const fdRowScaledMax = [0, 0, 0];
  for (let side = 0; side < 2; side++) for (let col = 0; col < 5; col++) {
    const index = [input.type === 1 ? 1 : 2, 3, 4, 6, 0][col]; const h = Math.max(1e-8, Math.abs(input.stations[side][index]) * 1e-5);
    const lo = structuredClone(input); const hi = structuredClone(input); lo.stations[side][index] -= h; hi.stations[side][index] += h;
    const l = assembleInterval(lo); const u = assembleInterval(hi);
    for (let row = 0; row < 3; row++) {
      const a = actual[side ? 'downstream' : 'upstream'][row][col]; const fd = -(u.residual[row] - l.residual[row]) / (2 * h);
      fdRowScaledMax[row] = Math.max(fdRowScaledMax[row], Math.abs(fd - a) / (1 + Math.abs(a)));
    }
  }
  fdRowScaledMax.forEach((v, row) => { if (!(input.type === 2 && row === 0)) assert.ok(v < 1e-6); });
  report.intervals.push({ name: input.name, nativeScaledMax, finiteDifferenceRowScaledMax: fdRowScaledMax,
    note: input.type === 2 ? 'Original BLDIF turbulent lag Jacobian omits the explicit UQ_RTA term. It is an approximate linearization; native parity is preserved.' : 'All active rows pass independent finite differences.' });
}
for (const panels of [80, 160, 240, 320]) {
  const r = solveCoupled({ elements: [{ points: naca4('0012', panels) }], alpha: 4 });
  assert.equal(r.status, 'solved'); assert.ok(r.diagnostics.couplingResidual < 1e-5);
  report.refinement.push({ panels, cl: r.cl, cd: r.cd, cm: r.cm, ...r.diagnostics });
}
for (let i = 2; i < report.refinement.length; i++) for (const k of ['cl', 'cd']) {
  assert.ok(Math.abs(report.refinement[i][k] - report.refinement[i - 1][k]) < Math.abs(report.refinement[i - 1][k] - report.refinement[i - 2][k]));
}
report.refinementNote = 'CL and CD successive changes decrease in this four-resolution study. Transition positions are reported without asserting monotonicity; moving transition intervals can be nonmonotone.';
for (const path of ['tests/fixtures/fortran/coupled.json', 'tests/fixtures/fortran/kernels.json',
  ...(await readdir('src/viscous', { recursive: true })).filter(p => p.endsWith('.js')).map(p => `src/viscous/${p}`)]) {
  report.sha256[path] = createHash('sha256').update(await readFile(path)).digest('hex');
}
await mkdir('docs', { recursive: true });
await writeFile('docs/coupled-validation-results.json', JSON.stringify(report, null, 2) + '\n');
console.log(`Coupled validation passed: ${report.cases.length} native cases, ${report.intervals.length} interval blocks, ${report.refinement.length} resolutions. See docs/coupled-validation-results.json.`);
