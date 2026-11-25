// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { naca4, naca4Standard, transform } from '../src/geometry/airfoil.js';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';

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
console.log('Native contour-to-Euler passed for four solves and the NLR 7301 start grid. Three cases passed fixed-edge and coupled BL checks.');
