// SPDX-License-Identifier: GPL-2.0-or-later
import { mkdir, writeFile } from 'node:fs/promises';
import { channelMesh } from '../src/euler/mesh.js';
import { multielementMesh } from '../src/euler/multielement-mesh.js';
import { solveEuler } from '../src/euler/solve.js';
import { solveStreamlineChannel } from '../src/euler/streamline.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { annularVortex, isentropicNozzle, weightedError } from '../tests/oracles/euler.js';

const record = r => ({ converged: r.converged, reason: r.reason, cells: r.mesh.cells.length,
  iterations: r.history.length - 1, diagnostics: r.diagnostics, grid: r.mesh.diagnostics });
const report = { formulation: 'First-order conservative finite-volume Euler verification backend; HLLC; simultaneous damped Newton',
  scope: 'Independent reference, not the MSES discretization. No BL coupling, force accuracy release, or validated shock capture.',
  normalization: 'c_ref, rho_inf=1, U_inf=1, gamma=1.4; entropy error is ln(p/rho^gamma) minus inlet value.',
  vortex: [], nozzle: [], airfoilRefinement: [], multielement: null, movingStreamlines: null };
for (const [nx, ny] of [[6, 2], [12, 4], [24, 8]]) {
  const { mesh, exact } = annularVortex(nx, ny);
  const r = solveEuler(mesh, { initial: mesh.cells.map(exact), alpha: 120 });
  report.vortex.push({ nx, ny, mach: .3, ...record(r), pressureL1: weightedError(r, exact, 'p'), velocityXL1: weightedError(r, exact, 'u') });
}
for (const nx of [12, 24, 48]) {
  const r = solveEuler(channelMesh({ nx, ny: 2, upper: x => isentropicNozzle(x).area }));
  report.nozzle.push({ nx, ny: 2, mach: .3, ...record(r), pressureL1: weightedError(r, c => isentropicNozzle(c.x), 'p'), velocityXL1: weightedError(r, c => isentropicNozzle(c.x), 'u'),
    caveat: 'Quasi-1D comparison includes finite-height model difference: height=.08, length=2, contraction=.12.' });
}
const moving = solveStreamlineChannel(channelMesh({ nx: 10, ny: 3, upper: x => 1 - .08 * Math.sin(Math.PI * x / 2) ** 2 }));
report.movingStreamlines = { ...record(moving), mach: .3, flowUnknowns: 120, gridUnknowns: 20,
  maxGridDisplacement: Math.max(...moving.displacements.map(Math.abs)) };
for (const [panels, rows] of [[20, 2], [40, 3], [80, 4]]) {
  const r = solveEuler(multielementMesh([naca4('0012', panels)], { rows, padding: 2 }), { mach: .2 });
  report.airfoilRefinement.push({ panels, rows, padding: 2, mach: .2, alpha: 0, ...record(r),
    pressureDragError: 2 * r.wallForces[0].x, liftSymmetryError: 2 * Math.abs(r.wallForces[0].y),
    caveat: 'Exact smooth subcritical inviscid drag is zero. This mesh family does not yet pass an airfoil-force accuracy gate.' });
}
const contours = [naca4('0012', 20), transform(naca4('0012', 20), { chord: .4, x: 1.05, y: -.2, angle: -10 })];
const multi = solveEuler(multielementMesh(contours, { rows: 2 }));
report.multielement = { ...record(multi), mach: .3, alpha: 0, contours,
  rawWallForces: multi.wallForces, connectedWakeCuts: multi.mesh.cuts.filter(c => c.type === 'wake-cut').length,
  caveat: 'Conservation/topology demonstration; finite-domain and first-order errors remain. Forces are diagnostics, not released predictions.' };
const cases = [...report.vortex, ...report.nozzle, ...report.airfoilRefinement, report.movingStreamlines, report.multielement];
if (cases.some(r => !r.converged || r.diagnostics.residual > 1e-9 || r.diagnostics.relativeMassImbalance > 1e-8)) throw new Error('Euler verification solve failed; report not updated.');
for (let i = 1; i < report.vortex.length; i++) if (report.vortex[i].pressureL1 >= .6 * report.vortex[i - 1].pressureL1) throw new Error('Vortex refinement gate failed.');
for (let i = 1; i < report.nozzle.length; i++) if (report.nozzle[i].velocityXL1 >= .8 * report.nozzle[i - 1].velocityXL1) throw new Error('Nozzle refinement gate failed.');
await mkdir('docs', { recursive: true });
await writeFile('docs/euler-validation-results.json', JSON.stringify(report, null, 2) + '\n');
console.log(`Euler verification: ${cases.length} converged solves; exact vortex pressure L1 ${report.vortex.map(r => r.pressureL1.toExponential(3)).join(' → ')}.`);
console.log(`Airfoil pressure-drag error (not released CD): ${report.airfoilRefinement.map(r => r.pressureDragError.toFixed(5)).join(' → ')}.`);
