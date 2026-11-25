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

const cases = [
  { label: 'sharp NACA 2412', elements: [sharp], alpha: 0, mach: .1, surface: 8, bl: true },
  { label: 'finite-base NACA 2412', elements: [finite], alpha: 0, mach: .1, surface: 8, bl: true, wakes: true },
  { label: 'stacked NACA 0012', elements: stacked, alpha: 0, mach: .1, surface: 6, bl: true, symmetric: true },
  { label: 'coarse transonic RAE 2822', elements: getBenchmarkAirfoil('rae2822-mses').elements,
    alpha: 2.68, mach: .74, surface: 6, transonic: true },
];
for (const c of cases) {
  const args = ['--topology', ...(c.bl ? ['--bl',1e6,9] : []), c.alpha, c.mach, c.surface, 4].map(String);
  const output = spawnSync(binary, args, { input: encode(c.elements), encoding: 'utf8' });
  assert.equal(output.status, 0, `${c.label}: ${output.stderr}${output.stdout}`);
  const r = JSON.parse(output.stdout);
  assert.equal(r.invalidCells, 0);
  assert.equal(r.flowSeedStatus, 0);
  assert.ok(r.minimumCornerSine > 0 && r.nodes > 100 && r.unknowns > r.nodes);
  assert.equal(r.wakeWidths > 0, Boolean(c.wakes));
  assert.equal(r.converged, true);
  assert.ok(r.iterations > 1 && r.residual <= 1e-9);
  assert.ok(r.maxMach > r.mach && (c.transonic ? r.maxMach > 1.1 : r.maxMach < .2));
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
assert.equal(refinedEuler.intervals, 40);
assert.equal(refinedEuler.converged, true);
assert.equal(refinedEuler.invalidCells, 0);
assert.ok(refinedEuler.residual <= 1e-9 && refinedEuler.maxMach > 1.1);

const refinedViscous = spawnSync(binary, ['--topology', '--potential-grid', '--coarse-start', '8',
  '--coupled', '2700000', '4', '--trip', '.03', '.07', '2.68', '.74', '16', '4'], {
  input: raeInput, encoding: 'utf8',
});
const refinedCoupled = JSON.parse(refinedViscous.stdout);
assert.equal(refinedViscous.status, refinedCoupled.converged ? 0 : 1);
assert.equal(refinedCoupled.coarseStartApplied, true);
assert.equal(refinedCoupled.solveStatus, 0);
assert.equal(refinedCoupled.boundaryLayer.complete, true);
const coupling = refinedCoupled.coupled;
assert.ok(coupling.updates > 0 && coupling.blResidual <= 1e-9);
assert.equal(Object.values(coupling.rejections).reduce((sum, count) => sum + count, 0), coupling.backtracks);
if (refinedCoupled.converged) {
  assert.ok(coupling.displacementResidual <= 1e-8 && coupling.edgeResidual <= 1e-8);
} else {
  assert.equal(coupling.usedLocalInteraction, true);
  assert.equal(refinedCoupled.boundaryLayer.edgeSource, 'interaction-iterate');
  assert.ok(coupling.edgeResidual < .25);
}

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

const mesh = spawnSync(binary, ['--topology', '--mesh-only', '6', '.185', '8', '4'], {
  input: encode(getBenchmarkAirfoil('nlr7301').elements), encoding: 'utf8',
});
assert.equal(mesh.status, 0, mesh.stderr + mesh.stdout);
const report = JSON.parse(mesh.stdout);
assert.equal(report.invalidCells, 0);
assert.equal(report.flowSeedStatus, 0);
assert.ok(report.wakeWidths > 0 && report.minimumCornerSine > 0);
assert.equal(report.solveStatus, undefined);

const rejected = spawnSync(binary, ['--topology', '0', '.1'], {
  input: encode([{ points: sharp.points.slice().reverse() }]), encoding: 'utf8',
});
assert.equal(rejected.status, 1);
assert.equal(JSON.parse(rejected.stdout).converged, undefined);
console.log('Native airfoil checks passed for coarse coupled cases, the refined 16-interval RAE 2822 Euler solve, and the fine RAE 2822 and NLR 7301 start grids.');
