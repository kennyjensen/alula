// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { naca4, naca4Standard, transform } from '../src/geometry/airfoil.js';
import { prepareAirfoilElement } from '../src/geometry/airfoil-element.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { builtinAirfoils } from '../src/geometry/builtin-airfoils.js';

const binary = fileURLToPath(new URL('./alula-euler', import.meta.url));
const sharp = { points: naca4('2412', 160) };
const openBase = naca4Standard('2412', 160);
const finite = {
  points: [...openBase, openBase[0]],
  trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: openBase.length - 1 },
};
const stacked = [-.65, .65].map(y => ({ points: transform(naca4('0012', 72), { x: 0, y }) }));
const encode = elements => [elements.length, ...elements.flatMap(({ points, trailingEdge: te }) => [
  `${points.length} ${te ? `F ${te.upperIndex} ${te.lowerIndex}` : 'S'}`,
  ...points.map(({ x, y }) => `${x} ${y}`),
])].join('\n') + '\n';

// Independent processes can select different OpenBLAS kernels on hybrid CPUs.
// Compare every field, allowing only roundoff in the floating-point results.
function compareRuns(actual, expected, path = 'report') {
  if (typeof expected === 'number') {
    assert.ok(typeof actual === 'number' && Math.abs(actual - expected) <= 1e-12 * Math.max(1, Math.abs(expected)),
      `${path}: ${actual} differs from ${expected}`);
  } else if (expected && typeof expected === 'object') {
    assert.ok(actual && typeof actual === 'object' && Array.isArray(actual) === Array.isArray(expected), path);
    assert.deepEqual(Object.keys(actual).sort(), Object.keys(expected).sort(), path);
    for (const key of Object.keys(expected)) compareRuns(actual[key], expected[key], `${path}.${key}`);
  } else {
    assert.equal(actual, expected, path);
  }
}

const cases = [
  { label: 'sharp NACA 2412', elements: [sharp], alpha: 0, mach: .1, surface: 8, bl: true },
  { label: 'finite-base NACA 2412', elements: [finite], alpha: 0, mach: .1, surface: 8, bl: true, wakes: true },
  { label: 'stacked NACA 0012', elements: stacked, alpha: 0, mach: .1, surface: 6, bl: true, symmetric: true },
  { label: 'coarse transonic RAE 2822', elements: getBenchmarkAirfoil('rae2822-mses').elements,
    alpha: 2.68, mach: .74, surface: 6, transonic: true },
];
for (const c of cases) {
  const args = ['--topology', '--cp', ...(c.bl ? ['--bl',1e6,9] : []), c.alpha, c.mach, c.surface, 4].map(String);
  const output = spawnSync(binary, args, { input: encode(c.elements), encoding: 'utf8' });
  assert.equal(output.status, 0, `${c.label}: ${output.stderr}${output.stdout}`);
  const r = JSON.parse(output.stdout);
  assert.equal(r.invalidCells, 0);
  assert.equal(r.flowSeedStatus, 0);
  assert.ok(r.minimumCornerSine > 0 && r.nodes > 100 && r.unknowns > r.nodes);
  assert.equal(r.wakeWidths > 0, Boolean(c.wakes));
  assert.equal(r.converged, true);
  assert.ok(r.iterations > 1 && r.residual <= 1e-9);
  if (r.backtracks) {
    const { chart, flowDomain, residual } = r.stepRejections;
    assert.equal(chart + flowDomain + residual, r.backtracks);
  }
  assert.ok(r.maxMach > r.mach && (c.transonic ? r.maxMach > 1.1 : r.maxMach < .2));
  assert.equal(r.cpProvisional, undefined);
  assert.ok(r.cp.every(p => [p.x, p.y, p.cp].every(Number.isFinite)));
  let cx = 0;
  let cy = 0;
  for (let e = 0; e < c.elements.length; ++e) {
    const [upper, lower, base] = ['upper', 'lower', 'base'].map(side =>
      r.cp.filter(p => p.contour === e && p.surface === side));
    assert.ok(upper.length >= c.surface + 1 && upper.length === lower.length);
    assert.equal(base.length > 0, Boolean(c.wakes));
    for (const key of ['x', 'y', 'cp']) assert.equal(upper[0][key], lower[0][key]);
    if (base.length) {
      assert.equal(base[0].cp, lower.at(-1).cp);
      assert.equal(base.at(-1).cp, upper.at(-1).cp);
    }
    for (const path of [upper.reverse(), lower, base]) for (let i = 1; i < path.length; ++i) {
      const a = path[i - 1], b = path[i], average = (a.cp + b.cp) / 2;
      cx -= average * (b.y - a.y);
      cy += average * (b.x - a.x);
    }
  }
  const alpha = c.alpha * Math.PI / 180;
  assert.ok(Math.abs(r.cl - (cy * Math.cos(alpha) - cx * Math.sin(alpha))) < 1e-12);
  assert.ok(Math.abs(r.pressureDrag - (cx * Math.cos(alpha) + cy * Math.sin(alpha))) < 1e-12);
  if (c.symmetric) assert.ok(Math.abs(r.circulation) < 1e-5);
  if (c.bl) {
    assert.equal(r.boundaryLayer.edgeSource,'converged-euler');
    assert.equal(r.boundaryLayer.complete,true);
    assert.equal(r.boundaryLayer.profiles.length,3*c.elements.length);
    for (const p of r.boundaryLayer.profiles) {
      assert.equal(p.status,0);
      assert.equal(p.completed,p.stations);
      assert.ok(p.maximumResidual <= 1e-9);
      assert.ok(p.displacementThickness > p.momentumThickness && p.momentumThickness > 0);
      if (p.side !== 'wake') assert.ok(p.transitionS > 0 && p.transitionS < 1.2);
    }

    const coupledArgs = ['--topology', '--coupled', 1e6, 9, c.alpha, c.mach, c.surface, 4].map(String);
    const output = spawnSync(binary, coupledArgs, { input: encode(c.elements), encoding: 'utf8' });
    assert.equal(output.status, 0, `${c.label} coupled: ${output.stderr}${output.stdout}`);
    const solution = JSON.parse(output.stdout);
    if (c === cases[0]) {
      const live = spawnSync(binary, ['--progress', ...coupledArgs], {
        input: encode(c.elements), encoding: 'utf8',
      });
      assert.equal(live.status, 0);
      const liveSolution = JSON.parse(live.stdout);
      compareRuns(liveSolution, solution);
      const events = live.stderr.trim().split('\n').map(line => JSON.parse(line));
      for (const stage of ['Euler', 'Euler/BL']) {
        const frames = events.filter(e => e.stage === stage && e.iteration > 0);
        assert.ok(frames.length > 0);
        assert.ok(frames.every(e => Number.isFinite(e.residual) && Number.isFinite(e.cl) &&
          Number.isFinite(e.cm) && Number.isFinite(e.pressureDrag) && e.cp.length > 0 &&
          e.streamlines.some(line => line.length > 1)));
      }
      assert.equal(events.at(-1).type, 'retained');
      assert.ok(events.filter(e=>e.stage==='Euler').every(e=>!e.boundaryLayerCurves));
      const blFrames=events.filter(e=>e.stage==='Euler/BL' && e.boundaryLayerCurves?.length);
      assert.ok(blFrames.length>0);
      const finalCurves=events.at(-1).boundaryLayerCurves;
      assert.equal(finalCurves.length,2);
      for (const curve of finalCurves) {
        const profile=liveSolution.boundaryLayer.profiles.find(p=>p.side===curve.surface);
        assert.equal(curve.contour,0);
        assert.equal(curve.samples.length,profile.completed);
        assert.ok(curve.samples.every(([s,d],i)=>s>0 && d>0 && (i===0 || s>curve.samples[i-1][0])));
        assert.ok(Math.abs(curve.samples.at(-1)[1]/profile.displacementThickness-1)<1e-8);
      }
      assert.equal(events.at(-2).cl, liveSolution.cl);
      assert.equal(events.at(-2).viscousWakeDrag, liveSolution.viscousWakeDrag);
      const desktop = spawnSync(binary, ['--live-progress', ...coupledArgs], {
        input: encode(c.elements), encoding: 'utf8',
      });
      assert.equal(desktop.status, 0);
      compareRuns(JSON.parse(desktop.stdout), liveSolution);
      const samples = desktop.stderr.trim().split('\n').map(line => JSON.parse(line));
      assert.equal(samples.length, events.length);
      for (const [i, sample] of samples.entries()) {
        if (sample.plotsDeferred) {
          const {cp, streamlines, boundaryLayerCurves, ...scalars} = events[i];
          compareRuns(sample, {...scalars, plotsDeferred: true}, `step ${i}`);
        } else compareRuns(sample, events[i], `step ${i}`);
      }
    }
    assert.equal(solution.converged, true);
    assert.equal(solution.coupled.status, 0);
    assert.ok(solution.coupled.updates > 0 && solution.coupled.eulerUpdates > 0);
    assert.ok(solution.residual <= 1e-9 && solution.coupled.blResidual <= 1e-9);
    assert.ok(solution.coupled.displacementResidual <= 1e-8);
    assert.equal(solution.boundaryLayer.edgeSource, 'displaced-euler');
    assert.equal(solution.boundaryLayer.complete, true);
    assert.ok([solution.cl, solution.cm, solution.pressureDrag].every(Number.isFinite));
    assert.ok(solution.viscousWakeDrag > 0);
    if (c.symmetric) assert.ok(Math.abs(solution.circulation) < 1e-5);
    else assert.ok(Math.abs(solution.circulation - r.circulation) > 1e-4);
  }
}

const raeInput = encode(cases.find(c => c.transonic).elements);
const refined = spawnSync(binary, ['--topology', '--potential-grid', '--coarse-start', '8',
  '2.68', '.74', '16', '4'], { input: raeInput, encoding: 'utf8' });
assert.equal(refined.status, 0, refined.stderr + refined.stdout);
const refinedEuler = JSON.parse(refined.stdout);
assert.equal(refinedEuler.coarseStartApplied, true);
assert.equal(refinedEuler.coarseStagesSolved, 1);
assert.equal(refinedEuler.intervals, 40);
assert.equal(refinedEuler.converged, true);
assert.equal(refinedEuler.invalidCells, 0);
assert.ok(refinedEuler.residual <= 1e-9 && refinedEuler.maxMach > 1.1);

const staged = spawnSync(binary, ['--topology', '--potential-grid', '--coarse-start', '8,16x4,20x12',
  '2.68', '.74', '24', '12'], { input: raeInput, encoding: 'utf8' });
assert.equal(staged.status, 0, staged.stderr + staged.stdout);
const stagedEuler = JSON.parse(staged.stdout);
assert.equal(stagedEuler.coarseStartApplied, true);
assert.equal(stagedEuler.coarseStagesSolved, 3);
assert.equal(stagedEuler.intervals, 56);
assert.equal(stagedEuler.converged, true);
assert.equal(stagedEuler.invalidCells, 0);
assert.ok(stagedEuler.residual <= 1e-9 && stagedEuler.maxMach > 1.1);

for (const stages of ['8,', '16x0', '16,8', '8,24x12']) {
  const invalid = spawnSync(binary, ['--coarse-start', stages, '2.68', '.74', '24', '12'], {
    input: raeInput, encoding: 'utf8',
  });
  assert.equal(invalid.status, 2, stages);
  assert.match(invalid.stderr, /^Usage:/);
}

const refinedViscous = spawnSync(binary, ['--topology', '--potential-grid', '--coarse-start', '8', '--cp',
  '--coupled', '2700000', '4', '--trip', '.03', '.07', '2.68', '.74', '16', '4'], {
  input: raeInput, encoding: 'utf8',
});
const refinedCoupled = JSON.parse(refinedViscous.stdout);
assert.equal(refinedViscous.status, 0, refinedViscous.stderr + refinedViscous.stdout);
assert.equal(refinedCoupled.coarseStartApplied, true);
assert.equal(refinedCoupled.converged, true);
assert.equal(refinedCoupled.solveStatus, 0);
assert.equal(refinedCoupled.invalidCells, 0);
assert.ok(refinedCoupled.residual <= 1e-9 && refinedCoupled.maxMach > 1);
assert.ok(refinedCoupled.cp.length > 0 && !refinedCoupled.cpProvisional);
assert.equal(refinedCoupled.boundaryLayer.complete, true);
assert.equal(refinedCoupled.boundaryLayer.edgeSource, 'displaced-euler');
const coupling = refinedCoupled.coupled;
assert.ok(coupling.updates > 0 && coupling.blResidual <= 1e-9);
assert.ok(coupling.displacementResidual <= 1e-8 && coupling.edgeResidual <= 1e-8);
assert.equal(Object.values(coupling.rejections).reduce((sum, count) => sum + count, 0), coupling.backtracks);
assert.equal(coupling.blResidual, Math.max(...refinedCoupled.boundaryLayer.profiles.map(p => p.maximumResidual)));
assert.ok(Math.abs(refinedCoupled.cl - .65903664) < 2e-5);
assert.ok(Math.abs(refinedCoupled.cm + .11248166) < 2e-5);
assert.ok(Math.abs(refinedCoupled.viscousWakeDrag - .01841677) < 2e-6);

for (const [intervals, crosslines] of [['6', 'algebraic'], ['8', 'potential']]) {
  const rae = spawnSync(binary, ['--topology', '--coupled', '2700000', '4', '--trip', '.03', '.07',
    '2.68', '.74', intervals, '4'], { input: raeInput, encoding: 'utf8' });
  assert.equal(rae.status, 0, rae.stderr + rae.stdout);
  const transonic = JSON.parse(rae.stdout);
  assert.equal(transonic.converged, true);
  assert.equal(transonic.crosslines, crosslines);
  assert.equal(transonic.coupled.usedInteraction, true);
  assert.ok(transonic.maxMach > 1.1 && transonic.residual <= 1e-9);
  assert.ok(transonic.coupled.blResidual <= 1e-9);
  assert.ok(transonic.coupled.displacementResidual <= 1e-8);
  assert.ok(transonic.coupled.edgeResidual <= 1e-8);
  assert.equal(transonic.boundaryLayer.complete, true);
  assert.ok(transonic.boundaryLayer.profiles[0].transitionS < transonic.boundaryLayer.profiles[1].transitionS);
  assert.ok(transonic.viscousWakeDrag > 0);
}

for (const tubes of ['4', '24']) {
  const fine = spawnSync(binary, ['--topology', '--potential-grid', '--mesh-only', '2.68', '.74', '64', tubes], {
    input: raeInput, encoding: 'utf8',
  });
  assert.equal(fine.status, 0, fine.stderr + fine.stdout);
  const report = JSON.parse(fine.stdout);
  assert.equal(report.crosslines, 'potential');
  assert.equal(report.invalidCells, 0);
  assert.equal(report.flowSeedStatus, 0);
  assert.ok(report.intervals > 64 && report.minimumCornerSine > .001);
}

const referenceDomain = spawnSync(binary, ['--topology', '--potential-grid', '--mesh-only',
  '--domain', '-2', '3', '-3', '3.5', '--inlet-intervals', '16', '--wake-intervals', '18',
  '2.68', '.74', '64', '24'], { input: raeInput, encoding: 'utf8' });
assert.equal(referenceDomain.status, 0, referenceDomain.stderr + referenceDomain.stdout);
const referenceMesh = JSON.parse(referenceDomain.stdout);
assert.equal(referenceMesh.intervals, 98);
assert.equal(referenceMesh.invalidCells, 0);
assert.equal(referenceMesh.flowSeedStatus, 0);

const fineEuler = spawnSync(binary, ['--topology', '--potential-grid',
  '--domain', '-2', '3', '-3', '3.5', '--inlet-intervals', '16', '--wake-intervals', '16',
  '2.68', '.74', '64', '24'], { input: raeInput, encoding: 'utf8' });
assert.equal(fineEuler.status, 0, fineEuler.stderr + fineEuler.stdout);
const referenceFlow = JSON.parse(fineEuler.stdout);
assert.equal(referenceFlow.converged, true);
assert.equal(referenceFlow.intervals, 96);
assert.equal(referenceFlow.invalidCells, 0);
assert.ok(referenceFlow.maxMach > 1 && referenceFlow.residual <= 1e-9);

function machContinuation(start, intervals, tubes, gridOptions = [], target = '.74') {
  const run = spawnSync(binary, ['--topology', '--potential-grid', ...gridOptions,
    '--mach-start', start, '--coupled', '2700000', '4', '--trip', '.03', '.07',
    '2.68', target, intervals, tubes], { input: raeInput, encoding: 'utf8' });
  return { status: run.status, result: JSON.parse(run.stdout) };
}

const ramp = machContinuation('.70', '8', '4');
assert.equal(ramp.status, 0);
assert.equal(ramp.result.startMach, .70);
assert.equal(ramp.result.actualMach, .74);
assert.equal(ramp.result.converged, true);
assert.ok(ramp.result.coupled.machStages > 1);
assert.ok(ramp.result.residual <= 1e-9 && ramp.result.coupled.blResidual <= 1e-9);

const higherMach = machContinuation('.74', '8', '4', ['--cp'], '.77');
assert.equal(higherMach.result.mach, .77);
assert.ok(higherMach.result.coupled.machStages > 0);
if (higherMach.result.converged) {
  assert.equal(higherMach.status, 0);
  assert.equal(higherMach.result.actualMach, .77);
} else {
  assert.equal(higherMach.status, 1);
  assert.ok(higherMach.result.actualMach > .74 && higherMach.result.actualMach < .77);
  assert.ok(higherMach.result.coupled.machRetries > 0);
  assert.equal(higherMach.result.finalGrid.invalidCells, 0);
  assert.equal(higherMach.result.cpProvisional, true);
}

if (process.env.ALULA_FINE_COUPLED === '1') {
  const fine = machContinuation('.64', '64', '24', ['--cp', '--coarse-start', '8', '--domain', '-2', '3', '-3', '3.5',
    '--inlet-intervals', '16', '--wake-intervals', '16', '--surface-uniformity', '.7']);
  const r = fine.result;
  assert.equal(fine.status, 0);
  assert.equal(r.converged, true);
  assert.equal(r.actualMach, .74);
  assert.equal(r.invalidCells, 0);
  assert.equal(r.boundaryLayer.complete, true);
  assert.ok(r.residual <= 1e-9 && r.coupled.blResidual <= 1e-9);
  assert.ok(r.coupled.displacementResidual <= 1e-8 && r.coupled.edgeResidual <= 1e-8);
  assert.ok(r.cp.length > 0 && r.cp.every(p => Number.isFinite(p.cp)) && !r.cpProvisional);
  assert.ok(Math.abs(r.cl - .7630) < .002 && Math.abs(r.cm + .09327) < .001);
  assert.ok(Math.abs(r.viscousWakeDrag - .01614) < .0002);
}

const clippedDomain = spawnSync(binary, ['--topology', '--mesh-only',
  '--domain', '-2', '.9', '-3', '3.5', '2.68', '.74'], { input: raeInput, encoding: 'utf8' });
assert.equal(clippedDomain.status, 2);
assert.match(clippedDomain.stderr, /strictly enclose/);

const flapInput = encode(builtinAirfoils.flap.map(p =>
  prepareAirfoilElement({ points: transform(naca4Standard(p.code, 160), p) })));
for (const coarseArgs of [[], ['--potential-grid', '--coarse-start', '7']]) {
  const flapped = spawnSync(binary, ['--topology', ...coarseArgs, '--coupled', '1000000', '9',
    '--trip', '.2', '.2', '--trip', '.2', '.2', '4', '.2', '8', '4'], { input: flapInput, encoding: 'utf8' });
  assert.equal(flapped.status, 0, flapped.stderr + flapped.stdout);
  const flapResult = JSON.parse(flapped.stdout);
  assert.equal(flapResult.converged, true);
  if (coarseArgs.length) assert.equal(flapResult.coarseStartApplied, true);
  assert.ok(flapResult.wakeWidths > 0 && flapResult.viscousWakeDrag > 0);
  assert.ok(flapResult.residual <= 1e-9 && flapResult.coupled.blResidual <= 1e-9);
  assert.ok(flapResult.coupled.displacementResidual <= 1e-8 && flapResult.coupled.edgeResidual <= 1e-8);
  assert.equal(flapResult.boundaryLayer.complete, true);
  assert.equal(flapResult.boundaryLayer.profiles.length, 6);
}

const trippedAssembly = spawnSync(binary, ['--topology', '--bl', '1000000', '20',
  '--trip', '.2', '.3', '--trip', '.5', '.6', '0', '.1', '6', '4'], {
  input: encode(stacked.slice().reverse()), encoding: 'utf8',
});
assert.equal(trippedAssembly.status, 0, trippedAssembly.stderr + trippedAssembly.stdout);
const [lowerUpper, lowerLower, , upperUpper, upperLower] = JSON.parse(trippedAssembly.stdout).boundaryLayer.profiles;
assert.ok(lowerUpper.transitionS > upperUpper.transitionS + .1);
assert.ok(lowerLower.transitionS > upperLower.transitionS + .1);

const nlrInput = encode(getBenchmarkAirfoil('nlr7301').elements);
const mesh = spawnSync(binary, ['--topology', '--mesh-only', '6', '.185', '8', '4'], {
  input: nlrInput, encoding: 'utf8',
});
assert.equal(mesh.status, 0, mesh.stderr + mesh.stdout);
const report = JSON.parse(mesh.stdout);
assert.equal(report.invalidCells, 0);
assert.equal(report.flowSeedStatus, 0);
assert.ok(report.wakeWidths > 0 && report.minimumCornerSine > 0);
assert.equal(report.solveStatus, undefined);

const nlrEulerGrids = [
  ['6', '4'], ['8', '4'],
  ['12', '6', ['--surface-uniformity', '.7']],
  ['32', '24'],
];
for (const [surface, tubes, spacing = []] of nlrEulerGrids) {
  const nlr = spawnSync(binary, ['--topology', '--potential-grid', ...spacing, '6', '.185', surface, tubes], {
    input: nlrInput, encoding: 'utf8',
  });
  assert.equal(nlr.status, 0, nlr.stderr + nlr.stdout);
  const solution = JSON.parse(nlr.stdout);
  assert.equal(solution.solveStatus, 0);
  assert.equal(solution.converged, true);
  assert.equal(solution.startupRecovery, 'wall-following');
  assert.equal(solution.startupGridRepaired, true);
  assert.equal(solution.invalidCells, 0);
  assert.ok(solution.residual <= 1e-9 && solution.wakeWidths > 0);
  assert.ok(Number.isFinite(solution.cl) && Number.isFinite(solution.cm));
  if (surface === '32') {
    assert.equal(solution.intervals, 105);
    assert.equal(solution.nodes, 7950);
    assert.ok(Math.abs(solution.cl - 2.66638) < .001 && Math.abs(solution.cm + .54545) < .001);
  }
}

const nlrViscous = spawnSync(binary, ['--topology', '--potential-grid', '--coupled', '2510000', '9',
  '6', '.185', '8', '4'], { input: nlrInput, encoding: 'utf8' });
assert.equal(nlrViscous.status, 0, nlrViscous.stderr + nlrViscous.stdout);
const nlrCoupled = JSON.parse(nlrViscous.stdout);
assert.equal(nlrCoupled.solveStatus, 0);
assert.equal(nlrCoupled.converged, true);
assert.equal(nlrCoupled.coupled.status, 0);
assert.ok(nlrCoupled.coupled.eulerGridRepairs > 0 && nlrCoupled.coupled.usedInteraction);
assert.ok(nlrCoupled.residual <= 1e-9 && nlrCoupled.coupled.blResidual <= 1e-9);
assert.ok(nlrCoupled.coupled.displacementResidual <= 1e-8 && nlrCoupled.coupled.edgeResidual <= 1e-8);
assert.equal(nlrCoupled.boundaryLayer.complete, true);
assert.equal(nlrCoupled.boundaryLayer.profiles.length, 6);
for (const profile of nlrCoupled.boundaryLayer.profiles) {
  assert.equal(profile.completed, profile.stations);
  if (profile.side !== 'wake') assert.ok(profile.transitionS > 0);
}
assert.ok(nlrCoupled.viscousWakeDrag > 0);

const nlrRefined = spawnSync(binary, ['--topology', '--potential-grid', '--coarse-start', '8x4',
  '6', '.185', '10', '4'], { input: nlrInput, encoding: 'utf8' });
assert.equal(nlrRefined.status, 0, nlrRefined.stderr + nlrRefined.stdout);
const nlrRefinedEuler = JSON.parse(nlrRefined.stdout);
assert.equal(nlrRefinedEuler.coarseStagesSolved, 1);
assert.equal(nlrRefinedEuler.coarseStartApplied, true);
assert.equal(nlrRefinedEuler.startupGridRepaired, true);
assert.equal(nlrRefinedEuler.solveStatus, 0);
assert.equal(nlrRefinedEuler.converged, true);
assert.ok(nlrRefinedEuler.residual <= 1e-9 && nlrRefinedEuler.wakeWidths > 0);

const nlrRefinedViscous = spawnSync(binary, ['--topology', '--potential-grid', '--coarse-start', '8x4',
  '--refine-coupled', '--coupled', '2510000', '9', '6', '.185', '10', '4'],
{ input: nlrInput, encoding: 'utf8' });
assert.equal(nlrRefinedViscous.status, 0, nlrRefinedViscous.stderr + nlrRefinedViscous.stdout);
const refinedNlr = JSON.parse(nlrRefinedViscous.stdout);
assert.equal(refinedNlr.converged, true);
assert.ok(refinedNlr.coupledRefinementFraction > .5 && refinedNlr.coupledRefinementFraction <= 1);
assert.ok(refinedNlr.residual <= 1e-9 && refinedNlr.coupled.blResidual <= 1e-9);
assert.ok(refinedNlr.coupled.displacementResidual <= 1e-8 && refinedNlr.coupled.edgeResidual <= 1e-8);
assert.equal(refinedNlr.boundaryLayer.complete, true);
// Numerical regression for the coupled root at this fixed resolution.
assert.ok(Math.abs(refinedNlr.cl - 2.4043444) < 2e-5);
assert.ok(Math.abs(refinedNlr.cm + .5532271) < 2e-5);
assert.ok(Math.abs(refinedNlr.viscousWakeDrag - .01181029) < 2e-6);

const rejected = spawnSync(binary, ['--topology', '0', '.1'], {
  input: encode([{ points: sharp.points.slice().reverse() }]), encoding: 'utf8',
});
assert.equal(rejected.status, 1);
assert.equal(JSON.parse(rejected.stdout).converged, undefined);
console.log('Native airfoil checks passed for coarse coupled cases, RAE 2822 Euler and coupled refinement, and NLR 7301 coupled refinement.');
