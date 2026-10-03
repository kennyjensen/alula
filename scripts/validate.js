// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { solveInviscid } from '../src/inviscid/linear-vortex.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { solveFlatPlate } from '../src/boundary-layer/tests/flat-plate.js';
import { joukowski } from '../tests/fixtures/joukowski.js';
import { blasius } from '../tests/fixtures/blasius.js';

const report = { generatedAt: new Date().toISOString(), runtime: process.version, status: 'passed',
  scope: 'Inviscid initializer and laminar flat-plate laboratory only; not MSES validation.', joukowski: [], multielement: [], flatPlate: [] };
for (const n of [80, 160, 320]) {
  const exact = joukowski(n); const r = solveInviscid({ elements: [exact], alpha: 4 });
  let sum = 0; let weight = 0;
  r.elements[0].cp.forEach((p, i) => { sum += (p.cp - exact.cpAtAngle(2 * Math.PI * (i + 0.5) / n)) ** 2 * p.length; weight += p.length; });
  report.joukowski.push({ panels: n, cl: r.cl, exactCl: exact.cl, clError: Math.abs(r.cl - exact.cl),
    cpRms: Math.sqrt(sum / weight), dragError: Math.abs(r.diagnostics.pressureDrag), ...r.diagnostics });
  assert.ok(r.diagnostics.normalVelocityResidual < 1e-10); assert.ok(r.diagnostics.kuttaResidual < 1e-10);
}
for (let i = 1; i < report.joukowski.length; i++) {
  for (const key of ['clError', 'cpRms', 'dragError']) assert.ok(report.joukowski[i][key] < report.joukowski[i - 1][key]);
}
assert.ok(report.joukowski.at(-1).clError < 0.002); assert.ok(report.joukowski.at(-1).cpRms < 0.01); assert.ok(report.joukowski.at(-1).dragError < 0.0011);
for (const n of [80, 160, 320]) {
  const r = solveInviscid({ elements: [{ name: 'Main', points: naca4('2412', n) },
    { name: 'Flap', points: transform(naca4('0012', n), { chord: 0.3, angle: -15, x: 0.94, y: -0.08 }) }], alpha: 4 });
  report.multielement.push({ panels: r.panelCount, cl: r.cl, cm: r.cm, ...r.diagnostics });
  assert.ok(r.diagnostics.normalVelocityResidual < 1e-10);
}
for (let i = 1; i < report.multielement.length; i++) {
  assert.ok(Math.abs(report.multielement[i].pressureDrag) < Math.abs(report.multielement[i - 1].pressureDrag));
  assert.ok(report.multielement[i].liftMismatch < report.multielement[i - 1].liftMismatch);
}
const exactBL = blasius(); report.blasiusOracle = exactBL;
for (const reynolds of [1e5, 1e6, 1e7]) {
  const r = solveFlatPlate({ reynolds }); assert.ok(r.converged);
  const p = r.stations.at(-1);
  const errors = { thetaRelative: Math.abs(p.theta * Math.sqrt(reynolds) / exactBL.theta - 1),
    hRelative: Math.abs(p.h / exactBL.h - 1), cfRelative: Math.abs(p.cf * Math.sqrt(reynolds) / exactBL.cf - 1) };
  assert.ok(Object.values(errors).every(v => v < 0.01));
  report.flatPlate.push({ reynolds, stationCount: r.stations.length, thetaAtX1: p.theta, h: p.h, cfAtX1: p.cf,
    errors, history: r.history });
}
const sources = ['src/inviscid/linear-vortex.js', 'src/inviscid/panel.js', 'src/numerics/linear.js', 'src/numerics/newton.js',
  'src/geometry/airfoil.js', 'src/boundary-layer/tests/closures.js', 'src/boundary-layer/tests/flat-plate.js',
  'tests/fixtures/blasius.js', 'tests/fixtures/joukowski.js', 'third_party/Xfoil/src/xblsys.f', 'third_party/mses/mses.pdf'];
report.sha256 = {};
for (const file of sources) report.sha256[file] = createHash('sha256').update(await readFile(file)).digest('hex');
await mkdir('docs', { recursive: true });
await writeFile('docs/validation-results.json', JSON.stringify(report, null, 2) + '\n');
console.table(report.joukowski.map(({ panels, clError, cpRms, dragError }) => ({ panels, clError, cpRms, dragError })));
console.log(`Blasius: theta error ${(report.flatPlate[1].errors.thetaRelative * 100).toFixed(3)}%, H error ${(report.flatPlate[1].errors.hRelative * 100).toFixed(3)}%.`);
console.log('All report gates passed. Wrote docs/validation-results.json.');
