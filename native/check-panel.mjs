// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { execFileSync, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { naca4, naca4Standard, transform } from '../src/geometry/airfoil.js';
import { solveInviscid } from '../src/inviscid/linear-vortex.js';

const binary = fileURLToPath(new URL('./alula-panel', import.meta.url));
const main = naca4('2412', 160);
const flap = transform(naca4('0012', 96), { chord: .30, angle: -25, x: .94, y: -.08 });
const angles = [-5, 0, 4, 10];
const sharp = points => ({ points });
const finite = points => ({
  points: [...points, points[0]],
  trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: points.length - 1 },
});
const encode = (elements, topology) => [elements.length, ...elements.flatMap(({ points, trailingEdge }) => {
  const header = !topology ? points.length : trailingEdge?.kind === 'finite-base' ?
    `${points.length} F ${trailingEdge.upperIndex} ${trailingEdge.lowerIndex}` : `${points.length} S`;
  return [header, ...points.map(({ x, y }) => `${x} ${y}`)];
})].join('\n') + '\n';

const nacaFinite = finite(naca4Standard('2412', 160));
const nlr = getBenchmarkAirfoil('nlr7301').elements;
const three = buildReliabilityCase({ preset: 'three', mode: 'streamtube-bl', changes: {} }).caseData.elements;
const rotateStart = (element, offset) => {
  const { points, trailingEdge } = element;
  const n = points.length - 1;
  return { points: [...points.slice(offset, n), ...points.slice(0, offset + 1)], trailingEdge: {
    kind: 'finite-base', upperIndex: (trailingEdge.upperIndex + n - offset) % n,
    lowerIndex: (trailingEdge.lowerIndex + n - offset) % n,
  } };
};
const cases = [
  ['sharp single', [sharp(main)]],
  ['sharp flap', [sharp(main), sharp(flap)]],
  ['sharp fine', [sharp(naca4('0012', 600))]],
  ['finite single', [nacaFinite]],
  ['mixed edges', [nacaFinite, sharp(flap)]],
  ['finite slat and flap', three],
  ['NLR 7301 with segmented bases', nlr],
  ['NLR 7301 with shifted contour starts', nlr.map(e => rotateStart(e, 47))],
];

for (const [name, elements] of cases) {
  const topology = elements.some(e => e.trailingEdge?.kind === 'finite-base');
  const input = encode(elements, topology);
  const output = execFileSync(binary, ['--cp', ...(topology ? ['--topology'] : []), ...angles.map(String)],
    { input, encoding: 'utf8' });
  const results = output.trim().split('\n').map(JSON.parse);
  assert.equal(results.length, angles.length);
  for (const [i, actual] of results.entries()) {
    const alpha = angles[i];
    const reference = solveInviscid({ elements, alpha });
    assert.equal(reference.status, 'solved');
    assert.equal(actual.alpha, alpha);
    assert.ok(actual.residualBound < 1e-9);
    for (const [key, value] of Object.entries({ cl: reference.cl, cm: reference.cm,
      pressureDrag: reference.diagnostics.pressureDrag, circulationLift: reference.diagnostics.circulationLift })) {
      assert.ok(Math.abs(actual[key] - value) < 1e-10, `${name} at ${alpha}°, ${key}: ${actual[key]} vs ${value}`);
    }
    const cp = reference.elements.flatMap(element => element.cp.map(p => p.cp));
    assert.equal(actual.cp.length, cp.length);
    for (const [j, value] of cp.entries())
      assert.ok(Math.abs(actual.cp[j] - value) < 1e-9, `${name} at ${alpha}°, Cp at panel ${j}: ${actual.cp[j]} vs ${value}`);
  }
}
assert.equal(spawnSync(binary, ['0'], { input: encode([sharp(main.slice().reverse())]) }).status, 1);
assert.equal(spawnSync(binary, ['nan'], { input: encode([sharp(main)]) }).status, 1);
console.log(`Native panel parity passed for ${cases.length * angles.length} single and multi-element solves.`);
