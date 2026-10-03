// SPDX-License-Identifier: GPL-2.0-or-later
// Compare the app stopping policy with continued strict equation convergence.
// Usage: node scripts/validate-coupled-stopping.js [single|flap|rae2822-mses ...]
import assert from 'node:assert/strict';
import { buildReliabilityCase } from './validation/solver-reliability-cases.js';
import { solveCoupledStreamtubeAssembly } from '../src/euler/streamtube-coupled-assembly.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { coupledConvergenceSatisfied } from '../src/euler/streamtube-coupled-convergence.js';
import { quadCoupledResultForDisplay } from '../src/ui/quad-coupled-result.js';

const cases = {
  single: {}, flap: {},
  'rae2822-mses': { alpha: 2.68, mach: .74, reynolds: 1e6, ncrit: 9, gridIntervals: 16,
    gridTubes: 7, gridInletIntervals: 16, gridOutletIntervals: 16 },
};
for (const preset of process.argv.length > 2 ? process.argv.slice(2) : Object.keys(cases)) {
  assert.ok(Object.hasOwn(cases, preset), `Unknown stopping calibration case: ${preset}`);
  const { caseData } = buildReliabilityCase({ preset, mode: 'streamtube-bl', changes: cases[preset] });
  const start = performance.now();
  const practical = solveCoupledStreamtubeAssembly(caseData, { direct: true, convergence: 'mses',
    maxIterations: 40, eulerMaxIterations: 40 });
  assert.equal(practical.converged, true, practical.reason);
  assert.equal(coupledConvergenceSatisfied(practical, 1e-10), true);
  const strict = solveCoupledStreamtubeIses(undefined, { ...practical.checkpoint.continuation,
    resume: practical.checkpoint, convergence: 'residual', maxIterations: 15, tolerance: 1e-10 });
  assert.equal(strict.converged, true, strict.reason);
  assert.equal(strict.residualConverged, true);
  assert.equal(practical.mesh.quality.valid && strict.mesh.quality.valid, true);
  const a = quadCoupledResultForDisplay(practical, caseData);
  const b = quadCoupledResultForDisplay({ ...practical, ...strict, convergence: undefined }, caseData);
  const differences = Object.fromEntries(['cl', 'cd', 'cm'].map(k => [k, Math.abs(a[k] - b[k])]));
  differences.cp = Math.max(...a.elements.flatMap((e, body) => e.cp.map((p, i) => Math.abs(p.cp - b.elements[body].cp[i].cp))));
  for (const k of ['cl', 'cd', 'cm']) assert.ok(differences[k] < 1e-5, `${preset} ${k}: ${differences[k]}`);
  assert.ok(differences.cp < 1e-4, `${preset} Cp: ${differences.cp}`);
  console.log(JSON.stringify({ preset, conditions: { alpha: caseData.alpha, mach: caseData.mach,
    reynolds: caseData.reynolds, ncrit: caseData.ncrit, gridIntervals: caseData.gridIntervals, gridTubes: caseData.gridTubes },
    reason: practical.reason, iterations: practical.history.length - 1, strictAdditionalIterations: strict.history.length - 1,
    families: practical.families, strictFamilies: strict.families, differences,
    milliseconds: performance.now() - start }));
}
