// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { execFileSync, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { getBenchmarkAirfoil } from '../src/geometry/benchmark-airfoils.js';
import { buildReliabilityCase } from '../scripts/validation/solver-reliability-cases.js';
import { naca4, naca4Standard, transform } from '../src/geometry/airfoil.js';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';
import { streamfunctionAt, streamfunctionBranchIncrement } from '../src/inviscid/streamfunction.js';
import { traceArcCurve } from '../src/inviscid/arc-curve.js';
import { tracePotentialCurve } from '../src/inviscid/potential-curve.js';
import { createSurfaceContourCurve } from '../src/geometry/contour-topology.js';

const binary = fileURLToPath(new URL('./alula-panel', import.meta.url));
const fieldBinary = fileURLToPath(new URL('./check-panel-field', import.meta.url));
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

let probeCount = 0;
let traceCount = 0;
for (const [name, elements] of cases.filter(([name]) => name !== 'sharp fine')) {
  const alpha = 4;
  const { field } = solveInviscid({ elements, alpha });
  const probes = [
    { x: -5, y: .27 }, { x: -.22, y: .01 }, { x: .08, y: .12 },
    { x: .52, y: -.16 }, { x: 1.28, y: -.065 }, { x: 6, y: 2 },
  ].map(point => ({ point }));
  probes.push({ point: field.panels[0].a, singularVelocity: true });
  const base = field.basePanels?.find(panel => panel.sourceStrength !== 0);
  if (base) {
    const d = base.cutDirection;
    const o = base.cutOrigin;
    const x = o.x + .6*d.x, y = o.y + .6*d.y;
    probes.push(
      { point: { x: x + .001*d.y, y: y - .001*d.x } },
      { point: { x: x - .001*d.y, y: y + .001*d.x }, crossesCut: true },
    );
  }
  const start = { x: -3, y: 0 };
  const traces = [
    { seed: { x: -.5, y: .36 }, initial: .1, potential: .1, end: 4, plane: 1.5 },
    { seed: { x: 1.5, y: .48 }, initial: .1, potential: .1, end: -4, plane: -.4 },
    { seed: { x: -.6, y: -.42 }, initial: .15, potential: -.5, end: 1.85, byArc: true },
    { seed: { x: -.6, y: .42 }, initial: .1, potential: -.3, end: 3, plane: 1.4, byArc: true },
    { seed: { x: -.6, y: .42 }, initial: 0, potential: 0, end: .1, plane: 1.4, missesPlane: true },
    { seed: { x: -.6, y: .42 }, initial: 0, potential: 0, end: 2, plane: 1.4, clip: .32 },
  ];
  const input = [
    `${elements.length} ${alpha} 0`,
    ...elements.flatMap(({ points, trailingEdge }) => [
      `${points.length} ${+(trailingEdge?.kind === 'finite-base')} ${trailingEdge?.upperIndex ?? 0} ${trailingEdge?.lowerIndex ?? 0}`,
      ...points.map(({ x, y }) => `${x} ${y}`),
    ]),
    probes.length,
    ...probes.map(({ point }, i) => {
      const a = probes[i - 1]?.point ?? start;
      return `${a.x} ${a.y} ${point.x} ${point.y}`;
    }),
    ...traces.map(({ byArc, seed, initial, potential, end, plane, clip }) =>
      [Number(!!byArc), seed.x, seed.y, initial, potential, end,
        Number(plane !== undefined), plane ?? 0, clip ?? 10].join(' ')),
  ].join('\n') + '\n';
  const rows = execFileSync(fieldBinary, { input, encoding: 'utf8' }).trim().split('\n');
  assert.equal(rows.length, probes.length + traces.length);
  for (const [i, row] of rows.slice(0, probes.length).entries()) {
    const { point, singularVelocity, crossesCut } = probes[i];
    const [hasVelocity, hasPsi, u, v, psi, jump] = row.split(' ').map(Number);
    assert.equal(hasVelocity, +!singularVelocity, `${name}, velocity availability at probe ${i}`);
    assert.equal(hasPsi, 1, `${name}, streamfunction availability at probe ${i}`);
    if (hasVelocity) {
      const expected = velocityAt(point, field);
      assert.ok(Math.hypot(u - expected.u, v - expected.v) < 2e-10,
        `${name}, velocity at probe ${i}: (${u}, ${v}) vs (${expected.u}, ${expected.v})`);
    }
    if (hasPsi) assert.ok(Math.abs(psi - streamfunctionAt(point, field)) < 2e-10,
      `${name}, streamfunction at probe ${i}`);
    const expectedJump = streamfunctionBranchIncrement([probes[i - 1]?.point ?? start, point], field);
    if (crossesCut) assert.ok(Math.abs(expectedJump) > 1e-8);
    assert.ok(Math.abs(jump - expectedJump) < 2e-12, `${name}, source branch at probe ${i}`);
  }
  probeCount += probes.length;
  for (const [i, row] of rows.slice(probes.length).entries()) {
    const { byArc, seed, initial, potential, end, plane, clip, missesPlane } = traces[i];
    const [status, accepted, rejected, x, y, parameter, phi] = row.split(' ').map(Number);
    if (clip !== undefined) {
      assert.equal(status, 3, `${name}, trace should reject a clipped segment`);
      assert.ok(accepted > 0 && rejected > 0 && x < clip);
      continue;
    }
    const controls = {
      seed, initialPotential: potential, velocity: point => velocityAt(point, field),
      tolerance: 1e-9, potentialTolerance: 1e-9, maxStep: .05, maxSpatialStep: .08,
      ...(plane !== undefined ? { plane: { normal: { x: 1, y: 0 }, offset: plane } } : {}),
    };
    const reference = byArc ? traceArcCurve({ ...controls, initialArc: initial, endArc: end }) :
      tracePotentialCurve({ ...controls, endPotential: end });
    const point = reference.points.at(-1);
    assert.equal(reference.converged, !missesPlane);
    assert.equal(status, missesPlane ? 6 : plane !== undefined ? 1 : 0, `${name}, trace status`);
    assert.ok(Math.hypot(x - point.x, y - point.y) < 2e-8, `${name}, trace position at case ${i}`);
    assert.ok(Math.abs(phi - point.potential) < 2e-8, `${name}, trace potential at case ${i}`);
    assert.ok(Math.abs(parameter - (byArc ? point.arc : point.potential)) < 2e-8);
    if (plane !== undefined && !missesPlane) assert.ok(Math.abs(x - plane) <= 1e-9);
    assert.ok(accepted > 0);
  }
  traceCount += traces.length;
}
console.log(`Native panel field parity passed for ${probeCount} velocity, streamfunction, and source-cut probes.`);
console.log(`Native panel streamline checks passed for ${traceCount} forward, backward, and clipped traces.`);

const guideCases = [
  ['sharp single', [sharp(main)], 4],
  ['sharp flap with reversed input order', [sharp(main), sharp(flap)], 4],
  ['mixed trailing edges', [nacaFinite, sharp(flap)], 4],
  ['three-element assembly', three, 4],
  ['RAE 2822', getBenchmarkAirfoil('rae2822-mses').elements, 2.68],
  ['NLR 7301 with flap', nlr, 6],
];
let bodiesChecked = 0;
for (const [name, elements, alpha] of guideCases) {
  const input = `${elements.length} ${alpha} 1\n` + elements.flatMap(({ points, trailingEdge }) => [
    `${points.length} ${+(trailingEdge?.kind === 'finite-base')} ${trailingEdge?.upperIndex ?? 0} ${trailingEdge?.lowerIndex ?? 0}`,
    ...points.map(({ x, y }) => `${x} ${y}`),
  ]).join('\n') + '\n';
  const rows = execFileSync(fieldBinary, { input, encoding: 'utf8' }).trim().split('\n')
    .map(row => row.split(' ').map(Number));
  assert.equal(rows[0][0], 0, `${name}, guides status ${rows[0]}`);
  const levels = rows[1];
  assert.equal(levels.length, elements.length + 2);
  assert.ok(levels.every((level, i) => i === 0 || level > levels[i - 1]));
  const panel = solveInviscid({ elements, alpha });
  const { field } = panel;
  assert.ok(Math.abs(levels[0] - streamfunctionAt({ x: -3, y: -3 }, field)) < 2e-10);
  assert.ok(Math.abs(levels.at(-1) - streamfunctionAt({ x: -3, y: 3 }, field)) < 2e-10);
  const seen = new Set();
  for (const [i, row] of rows.slice(2, 2 + elements.length).entries()) {
    const [element, stag, strain, wall, seedX, seedY, inletY, wakeSeedX, wakeSeedY, wakeX, wakeY,
      incomingMidY, wakeMidY] = row;
    assert.ok(!seen.has(element) && elements[element]);
    seen.add(element);
    const body = elements[element];
    const curve = createSurfaceContourCurve(body.points, body);
    const samples = panel.elements[element].cp.filter(p => !p.base);
    const j = samples.findIndex((p, k) => k && samples[k - 1].qt < 0 && p.qt >= 0);
    assert.ok(j > 0);
    const left = .5 * (curve.knots[j - 1] + curve.knots[j]);
    const right = .5 * (curve.knots[j] + curve.knots[j + 1]);
    const slope = (samples[j].qt - samples[j - 1].qt) / (right - left);
    const expected = left - samples[j - 1].qt / slope;
    const atStag = curve.evaluate(expected);
    const wallPsi = panel.diagnostics.surfaceStreamfunctions?.[element] ?? streamfunctionAt(atStag.point, field);
    assert.ok(Math.abs(stag - expected) < 2e-10, `${name}, stagnation`);
    assert.ok(Math.abs(strain - slope / Math.hypot(atStag.derivative.x, atStag.derivative.y)) < 2e-8, `${name}, strain`);
    assert.ok(Math.abs(wall - wallPsi) < 2e-10);
    assert.ok(Math.abs(streamfunctionAt({ x: seedX, y: seedY }, field) - wall) < 2e-10, `${name}, projected stagnation seed`);
    assert.ok(Math.abs(streamfunctionAt({ x: -3, y: inletY }, field) - levels[i + 1]) < 2e-8, `${name}, captured level drift`);
    for (const [seed, plane, y] of [
      [{ x: seedX, y: seedY }, -3, inletY],
      [{ x: seedX, y: seedY }, -.6, incomingMidY],
      [{ x: wakeSeedX, y: wakeSeedY }, 4, wakeY],
      [{ x: wakeSeedX, y: wakeSeedY }, 1.6, wakeMidY],
    ]) {
      const reference = tracePotentialCurve({
        seed, endPotential: Math.sign(plane - seed.x) * 28, velocity: p => velocityAt(p, field),
        tolerance: 2e-9, plane: { normal: { x: 1, y: 0 }, offset: plane },
      });
      assert.ok(reference.converged, `${name}, JavaScript trace: ${reference.reason}`);
      assert.ok(Math.abs(reference.points.at(-1).y - y) < 2e-8, `${name}, guide plane intersection`);
    }
    assert.ok(Math.abs(wakeX - 4) <= 2e-9);
  }
  assert.equal(seen.size, elements.length);
  bodiesChecked += seen.size;
  const passages = rows.slice(2 + elements.length);
  assert.equal(passages.length, elements.length + 1);
  for (const [i, [level, inletY, outletY]] of passages.entries()) {
    assert.ok(level > levels[i] && level < levels[i + 1]);
    const seed = { x: -3, y: inletY };
    assert.ok(Math.abs(streamfunctionAt(seed, field) - level) < 2e-10);
    const reference = tracePotentialCurve({ seed, endPotential: 28, velocity: p => velocityAt(p, field),
      tolerance: 2e-9, plane: { normal: { x: 1, y: 0 }, offset: 4 } });
    assert.ok(reference.converged, `${name}, interior passage trace: ${reference.reason}`);
    assert.ok(Math.abs(reference.points.at(-1).y - outletY) < 2e-8, `${name}, interior plane intersection`);
  }
}
console.log(`Native guide parity passed for ${bodiesChecked} bodies and ${bodiesChecked + guideCases.length} passages across ${guideCases.length} assemblies.`);
