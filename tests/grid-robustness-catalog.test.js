// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { gridRobustnessCases, gridRobustnessControls, gridRobustnessPresets,
  buildGridRobustnessCase, planGridRobustnessCase, assertGridRobustnessGuiCoverage } from '../scripts/validation/grid-robustness-cases.js';
import { parseGridRobustnessArgs, runGridRobustness, assessGridRobustnessSolve } from '../scripts/check-grid-robustness.js';

test('both benchmark geometries cover every visible surface/tube pair with SLOR on and off', () => {
  const ids = gridRobustnessCases.map(c => c.id); assert.equal(ids.length, new Set(ids).size);
  for (const preset of gridRobustnessPresets) for (const n of gridRobustnessControls.surfaceIntervals)
    for (const tubes of gridRobustnessControls.tubes) for (const smoothing of [false, true]) {
      const matching = gridRobustnessCases.filter(c => c.tier === 'surface-tube-matrix' && c.preset === preset
        && c.changes.gridIntervals === n && c.changes.gridTubes === tubes && c.changes.gridEllipticSmoothing === smoothing);
      assert.equal(matching.length, 1);
    }
  assert.equal(gridRobustnessCases.filter(c => c.tier === 'surface-tube-matrix').length, 80);
});

test('coverage checks track current visible values and reject an uncovered addition', () => {
  const html = fs.readFileSync('index.html', 'utf8'); assert.equal(assertGridRobustnessGuiCoverage(html), true);
  assert.throws(() => assertGridRobustnessGuiCoverage(html.replace(/(<select\b[^>]*\bid="grid-intervals"[^>]*>)/,
    (_, tag) => `${tag}<option value="256">256</option>`)), /grid-intervals/);
  assert.throws(() => assertGridRobustnessGuiCoverage(html.replace(/(<input\b[^>]*\bid="grid-upper-streamlines"[^>]*\bmax=")32"/,
    (_, prefix) => `${prefix}40"`)), /grid-upper-streamlines/);
});

test('every inventory case plans without producing mesh, gas or convergence claims', () => {
  for (const spec of gridRobustnessCases) {
    const input = buildGridRobustnessCase(spec), before = structuredClone(input), plan = planGridRobustnessCase(input);
    assert.deepEqual(input, before); assert.equal(plan.actualMesh, null);
    assert.equal(plan.geometricAdmissibility, 'not-run'); assert.equal(plan.gasAdmissibility, 'not-run'); assert.equal(plan.convergence, 'not-run');
    assert.equal(plan.baseTubes.length, input.elements.length + 1);
    assert.equal(plan.inletIntervals, input.gridInletIntervals ?? 2 * input.gridIntervals);
    assert.equal(plan.outletIntervals, input.gridOutletIntervals ?? 2 * input.gridIntervals);
  }
});

test('all supported explicit streamline counts convert to tubes and advanced endpoints remain covered', () => {
  const source = buildGridRobustnessCase(gridRobustnessCases.find(c => c.preset === 'nlr7301'));
  for (let n = 4; n <= 32; n++) {
    const p = planGridRobustnessCase({ ...source, gridLowerTubes: n - 1, gridUpperTubes: n - 1, gridGapTubes: n - 1 });
    assert.deepEqual(p.baseStreamlines, [n, n, n]);
  }
  for (const preset of gridRobustnessPresets) for (const key of ['gridInletIntervals', 'gridOutletIntervals'])
    for (const n of [16, 32, 64, 128]) assert(gridRobustnessCases.some(c => c.preset === preset && c.changes[key] === n));
  for (const key of ['gridUpperTubes', 'gridLowerTubes', 'gridGapTubes']) for (const n of [3, 31])
    assert(gridRobustnessCases.some(c => c.changes[key] === n));
  assert.throws(() => planGridRobustnessCase({ ...source, gridUpperTubes: 32 }), /region tube/);
  assert.throws(() => planGridRobustnessCase({ ...source, gridInletIntervals: 3 }), /inlet/);
  assert.throws(() => planGridRobustnessCase({ ...source, gridIntervals: 140 }), /case controls/);
});

test('case generation preserves reported physical conditions and owns independent geometry', () => {
  const rae = gridRobustnessCases.find(c => c.id === 'rae2822-mses-128x7-slor-on');
  const a = buildGridRobustnessCase(rae), b = buildGridRobustnessCase(rae);
  assert.equal(a.mach, .74); assert.equal(a.alpha, 2.68); assert.equal(a.reynolds, 2510000); assert.equal(a.ncrit, 4);
  assert.equal(a.eulerIsmom, 4); assert.equal(a.transitionMode, 'automatic'); assert.deepEqual(a.materialTrips, [[1, 1]]);
  a.elements[0].points[0].y += 1; assert.notDeepEqual(a.elements, b.elements);
  const nlr = buildGridRobustnessCase(gridRobustnessCases.find(c => c.id === 'nlr7301-32x9-slor-on'));
  assert.equal(nlr.mach, .185); assert.equal(nlr.reynolds, 2700000); assert.equal(nlr.ncrit, 9); assert.equal(nlr.alpha, 6);
  assert.equal(nlr.elements.length, 2); assert(nlr.elements.every(e => e.trailingEdge.kind === 'finite-base'));
});

test('default command is cheap and costly stages require explicit case selection', () => {
  assert.equal(parseGridRobustnessArgs([]).stage, 'plan');
  for (const stage of ['topology', 'geometry', 'startup', 'solve'])
    assert.throws(() => parseGridRobustnessArgs([`--stage=${stage}`]), /Name --case/);
  assert.equal(parseGridRobustnessArgs(['--stage=geometry', '--all']).stage, 'geometry');
  assert.throws(() => parseGridRobustnessArgs(['--seconds=601']), /seconds/);
  assert.throws(() => parseGridRobustnessArgs(['--case=missing']), /Unknown/);
  assert.throws(() => parseGridRobustnessArgs(['--max-iterations=1']), /only to solve/);
});

test('retained lower Ncrit roots cannot pass the requested nonlinear-convergence gate', () => {
  const input = buildGridRobustnessCase(gridRobustnessCases.find(c => c.preset === 'nlr7301'));
  const stations = [{ theta: .001, deltaStar: .002, ue: 1 }, { theta: .001, deltaStar: .002, ue: 1 }];
  const result = { converged: true, mach: input.mach, alpha: input.alpha, cl: 1, cm: -.1, cd: .01,
    families: { euler: 1e-12, boundaryLayer: 1e-12, edgeMatching: 1e-12 }, mesh: { quality: { valid: true } },
    checkpoint: { restart: { options: { ncrit: 8 } } }, boundaryLayer: {
      surfaces: input.elements.flatMap((_, element) => ['upper', 'lower'].map(side => ({ element, side, stations }))),
      wakes: input.elements.map((_, element) => ({ element, stations })),
    } };
  assert.equal(assessGridRobustnessSolve(input, result).passed, false);
  result.checkpoint.restart.options.ncrit = 9;
  assert.equal(assessGridRobustnessSolve(input, result).passed, true);
  result.mach = .1; assert.equal(assessGridRobustnessSolve(input, result).passed, false);
});

test('bounded plan command writes immutable stage-scoped evidence and no numerical stage files', async () => {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'grid-robustness-plan-')), output = path.join(directory, 'receipt');
  const args = parseGridRobustnessArgs(['--case=rae2822-mses-128x7-slor-on', `--out=${output}`]);
  try {
    const report = await runGridRobustness(args);
    assert.equal(report.passed, true); assert.equal(report.currentConvergenceEnvelopeCertified, false);
    assert.equal(report.cases.length, 1); assert.equal(report.cases[0].gates.plan.status, 'pass');
    for (const [name, gate] of Object.entries(report.cases[0].gates)) if (name !== 'plan') assert.equal(gate.status, 'not-run');
    assert.deepEqual(fs.readdirSync(path.join(output, 'rae2822-mses-128x7-slor-on')), ['input.json']);
    await assert.rejects(runGridRobustness(args), /Preserve existing evidence/);
  } finally { fs.rmSync(directory, { recursive: true, force: true }); }
});

test('bounded raw-topology stage runs without constructing a panel/SLOR mesh or gas seed', async () => {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), 'grid-robustness-topology-')), output = path.join(directory, 'receipt');
  const args = parseGridRobustnessArgs(['--stage=topology', '--case=rae2822-mses-8x7-slor-off', `--out=${output}`, '--seconds=15']);
  try {
    const report = await runGridRobustness(args), entry = report.cases[0];
    assert.equal(report.passed, true); assert.equal(entry.gates.plan.status, 'pass'); assert.equal(entry.gates.topology.status, 'pass');
    assert.equal(entry.gates.geometry.status, 'not-run'); assert.equal(entry.gates.gasStart.status, 'not-run');
    assert.equal(entry.gates.nonlinearConvergence.status, 'not-run');
    assert.deepEqual(fs.readdirSync(path.join(output, entry.id)).sort(), ['input.json', 'outcome.json', 'topology.json']);
  } finally { fs.rmSync(directory, { recursive: true, force: true }); }
});
