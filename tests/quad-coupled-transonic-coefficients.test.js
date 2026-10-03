// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledTransonicCoefficients, quadCoupledTransonicCoefficientsFromResult,
  historicalCommonIsentropeWakeGas } from '../src/ui/quad-coupled-transonic-coefficients.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';

const close = (a, b, tolerance = 3e-12) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

// Two closed manufactured diamonds, with independent physical pressure and
// gas data. This fixture is not an Euler/BL root and invokes no solver.
function fixture() {
  const nx = 4, tubes = [1, 1, 1], mach = .4, gamma = 1.4, reynolds = 1e4;
  const pInf = 1 / (gamma * mach * mach), h0 = 1 / ((gamma - 1) * mach * mach) + .5;
  const bodies = [0, 1].map(body => ({ element: 1 - body, leadingIndex: 1, trailingIndex: 3 }));
  const lower = b => [0, 1, 2, 3, 4].map(i => ({ x: i - 1, y: 3 * b - (i === 2 ? .5 : 0) }));
  const upper = b => lower(b).map((p, i) => ({ ...p, y: i === 2 ? p.y + 1 : p.y }));
  const nodes = tubes.map((_, g) => Array.from({ length: nx + 1 }, (_, i) => [
    g === 0 ? { x: i - 1, y: -2 } : upper(g - 1)[i],
    g === 2 ? { x: i - 1, y: 6 } : lower(g)[i],
  ]));
  const gas = { rho: 1, q: 1, p: pInf, enthalpy: h0 - .5, machSquared: mach * mach };
  const cells = Array.from({ length: nx - 1 }, (_, k) => tubes.map((_, g) => [{
    interfacePressure: { lower: pInf - .03 * nodes[g][k + 1][0].y,
      upper: pInf - .03 * nodes[g][k + 1][1].y }, geometry: { normalAreas: [1, 1] },
  }]));
  const flow = { nodes, undisplacedNodes: structuredClone(nodes), cells,
    sections: Array.from({ length: nx }, () => tubes.map(() => [{ ...gas }])),
    allocation: { groups: [2, 3, 3].map(m => [{ massFlow: m }]) } };
  const stations = [], surfaces = [], wakes = [], packed = [];
  const append = station => {
    const id = stations.length; stations.push({ ...station, id });
    packed.push(.03, .1, .15, 1); return id; // theta=.001, delta*=.0015 in solver-length units.
  };
  bodies.forEach((body, b) => {
    for (const side of ['upper', 'lower']) surfaces.push({ body: b, side,
      ids: [2, 3].map(i => append({ body: b, side, i, kind: 'surface' })) });
    wakes.push({ body: b, ids: [3, 4].map(i => append({ body: b, i, kind: 'wake' })) });
  });
  const input = { flowModel: 'compressible', streamwiseMode: 'hybrid', mach, gamma, alpha: 0, bodies,
    outerLower: nodes[0].map(row => row[0]), weights: tubes.map(n => Array(n).fill(1)),
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  return { checkpoint: { version: 1, restart: { input, initialBL: packed,
    options: { reynolds, ncrit: 9, edgeMatching: 'section-velocity', blThermodynamics: 'historical-common-isentrope' } } },
    flow, bl: { scale: 1 / Math.sqrt(reynolds), stations, surfaces, wakes }, bodies,
    solverLength: 2, referenceChord: 2, momentReference: { x: .5, y: 0 }, alpha: 0, mach, gamma };
}

test('compact live frames preserve all coefficients and keep valid lift/moment when drag is unavailable', () => {
  const input = fixture(), before = structuredClone(input);
  const compact = structuredClone({ ...input, flow: observableFlow(input.flow), bl: observableBL(input.bl) });
  const expected = quadCoupledTransonicCoefficients(input), actual = quadCoupledTransonicCoefficients(compact);
  for (const key of ['cl', 'cd', 'cm', 'viscousDragCoefficient', 'eulerWaveDragCoefficient']) assert.equal(actual[key], expected[key]);
  assert.equal(compact.flow.sections[0], null);
  delete compact.flow.sections;
  const partial = quadCoupledTransonicCoefficients(compact);
  assert.equal(partial.cl, expected.cl); assert.equal(partial.cm, expected.cm); assert.equal(partial.cd, null);
  assert.deepEqual(input, before);
});

test('nominal physical-pressure coefficients and both wake losses use the correct length normalization', () => {
  const input = fixture(), before = structuredClone(input), r = quadCoupledTransonicCoefficients(input);
  close(r.cl, .06); close(r.cm, -.015); close(r.pressureIntegralDrag, 0);
  // Two wakes: thetaPhysical=.001*solverLength=.002. At Ue=1,
  // rho=1 and p=pInfinity, so Cd_v=2*(.002+.002)/referenceChord=.004.
  close(r.viscousDragCoefficient, .004); close(r.eulerWaveDragCoefficient, 0);
  assert.equal(r.cd, r.viscousDragCoefficient + r.eulerWaveDragCoefficient);
  assert.equal(r.decomposition.includesPressureIntegral, false);
  assert.equal(r.decomposition.includesAdditionalSkinFriction, false);
  assert.equal(r.eulerExit.massFlow, 8, 'Euler physical tube mass must not be multiplied by solverLength.');
  assert.equal(r.viscousExit.gasSource, 'historical-common-isentrope');
  assert.deepEqual(r.wakes.map(w => w.element), [0, 1]);
  assert.deepEqual(r.elements.map(e => e.element), [0, 1]);
  assert.equal(r.physicalAcceptance, false); assert.deepEqual(input, before);
});

test('common-isentrope wake gas agrees with the native BL station gas below and above sonic speed', () => {
  for (const mach of [.2, .4, .7]) for (const localMach of [.7, 1, 1.3]) {
    const gamma = 1.4, h0 = 1 / ((gamma - 1) * mach * mach) + .5;
    const h = h0 / (1 + .5 * (gamma - 1) * localMach * localMach), ue = Math.sqrt(2 * (h0 - h));
    const theta = .002, H = (1 + .113 * localMach * localMach) * 1.6 + .29 * localMach * localMach;
    const kernel = createIntegralKernel({ mach, gamma, reynolds: 1e6 });
    const native = kernel.station({ s: 1.3, theta, deltaStar: H * theta, ue, aux: .03 }, 'wake');
    const gas = historicalCommonIsentropeWakeGas(ue, { mach, gamma });
    close(gas.density, native.rho); close(gas.machSquared, native.machSquared);
    close(gas.pressure, (gamma - 1) / gamma * native.rho * h);
    assert.equal(gas.gasSource, 'historical-common-isentrope');
  }
});

test('physical wall pressure controls lift even with unchanged stored BL edge speeds', () => {
  const input = fixture(), first = quadCoupledTransonicCoefficients(input);
  // Double the transverse physical pressure gradient, preserving all BL data.
  for (let k = 0; k < input.flow.cells.length; k++) for (let g = 0; g < input.flow.cells[k].length; g++)
    for (const [side, bank] of [['lower', 0], ['upper', 1]])
      input.flow.cells[k][g][0].interfacePressure[side] -= .03 * input.flow.undisplacedNodes[g][k + 1][bank].y;
  const changed = quadCoupledTransonicCoefficients(input);
  close(changed.cl, 2 * first.cl); close(changed.cm, 2 * first.cm);
  assert.equal(changed.cd, first.cd);
  for (const wake of changed.wakeGas) assert.equal(wake.speed, 1);
});

test('a manufactured normal-shock entropy loss stays in Euler exits and is added once to viscous loss', () => {
  const input = fixture(), gamma = input.gamma, gm1 = gamma - 1, M1 = 1.4;
  const h0 = 1 / (gm1 * input.mach ** 2) + .5, h1 = h0 / (1 + .5 * gm1 * M1 ** 2);
  const q1 = Math.sqrt(gm1 * h1 * M1 ** 2);
  const upstream = historicalCommonIsentropeWakeGas(q1, input);
  const densityRatio = (gamma + 1) * M1 ** 2 / (gm1 * M1 ** 2 + 2);
  const pressureRatio = 1 + 2 * gamma / (gamma + 1) * (M1 ** 2 - 1);
  const q2 = q1 / densityRatio, rho2 = upstream.density * densityRatio, p2 = upstream.pressure * pressureRatio;
  const downstream = { q: q2, rho: rho2, p: p2, enthalpy: h0 - .5 * q2 * q2,
    machSquared: rho2 * q2 * q2 / (gamma * p2) };
  input.flow.sections.at(-1)[1][0] = downstream;
  input.flow.allocation.groups[1][0].massFlow = 3 * rho2 * q2;
  const r = quadCoupledTransonicCoefficients(input), eta = p2 / (upstream.pressure)
    * (1 + .5 * gm1 * downstream.machSquared) ** (gamma / gm1)
    / (1 + .5 * gm1 * M1 * M1) ** (gamma / gm1);
  const entropy = -Math.log(eta), pInf = 1 / (gamma * input.mach ** 2);
  const hRecovered = gamma / gm1 * p2 / rho2 * (pInf / p2) ** (gm1 / gamma);
  const qRecovered = Math.sqrt(2 * (h0 - hRecovered));
  const expectedWave = input.flow.allocation.groups[1][0].massFlow * (1 - qRecovered) / (.5 * input.referenceChord);
  assert.ok(entropy > 0); close(r.eulerExit.sections[1].entropyOverR, entropy);
  assert.ok(r.eulerWaveDragCoefficient > 0); close(r.eulerWaveDragCoefficient, expectedWave);
  close(r.viscousDragCoefficient, .004); close(r.cd, expectedWave + .004);
  // Stored common-isentrope BL Ue remains 1 despite physical Euler loss.
  for (const w of r.wakeGas) { close(w.density, 1); close(w.pressure, pInf); }
});

test('geometric unit scaling preserves coefficients and scales each physical Euler mass only once', () => {
  const input = fixture(), before = quadCoupledTransonicCoefficients(input), k = 3;
  const scaled = structuredClone(input), point = p => ({ x: k * p.x, y: k * p.y });
  scaled.solverLength *= k; scaled.referenceChord *= k; scaled.momentReference = point(scaled.momentReference);
  for (const key of ['nodes', 'undisplacedNodes']) scaled.flow[key] = scaled.flow[key].map(g => g.map(row => row.map(point)));
  scaled.flow.allocation.groups.forEach(g => g.forEach(t => { t.massFlow *= k; }));
  const result = quadCoupledTransonicCoefficients(scaled);
  for (const key of ['cl', 'cm', 'cd', 'pressureIntegralDrag', 'viscousDragCoefficient', 'eulerWaveDragCoefficient']) close(result[key], before[key]);
  close(result.eulerExit.massFlow, k * before.eulerExit.massFlow);
  close(result.wakeGas[0].theta, k * before.wakeGas[0].theta);
  const chord = quadCoupledTransonicCoefficients({ ...input, referenceChord: 2 * input.referenceChord });
  close(chord.cl, before.cl / 2); close(chord.cd, before.cd / 2); close(chord.cm, before.cm / 4);
});

test('small negative physical Euler entropy contributions remain signed in the total', () => {
  const input = fixture(), state = input.flow.sections.at(-1)[0][0], factor = 1.00001;
  state.rho *= factor; state.p *= factor; input.flow.allocation.groups[0][0].massFlow *= factor;
  const r = quadCoupledTransonicCoefficients(input);
  assert.ok(r.eulerExit.sections[0].entropyOverR < 0);
  assert.ok(r.eulerWaveDragCoefficient < 0);
  assert.equal(r.cd, r.viscousDragCoefficient + r.eulerWaveDragCoefficient);
  assert.ok(r.cd < r.viscousDragCoefficient);
});

test('one invalid or missing wake makes total drag unavailable, preserving physical lift and valid Euler diagnostics', () => {
  for (const mutate of [x => { x.checkpoint.restart.initialBL[4 * x.bl.wakes[1].ids.at(-1) + 1] = 0; },
    x => { x.bl.wakes.pop(); }, x => { x.bl.wakes[1].body = x.bl.wakes[0].body; }]) {
    const input = fixture(); mutate(input); const r = quadCoupledTransonicCoefficients(input);
    assert.equal(r.cd, null); assert.equal(r.viscousDragCoefficient, null); assert.equal(r.decomposition.complete, false);
    assert.ok(Number.isFinite(r.cl)); assert.ok(Number.isFinite(r.eulerWaveDragCoefficient));
    assert.equal(r.wakes.length, 0); assert.ok(r.errors.viscousExit); assert.ok(r.warnings.some(w => /Viscous exit/.test(w)));
  }
  const input = fixture(); input.flow.sections.at(-1)[0][0].rho = 0;
  const invalidEuler = quadCoupledTransonicCoefficients(input);
  assert.equal(invalidEuler.cd, null); assert.equal(invalidEuler.eulerWaveDragCoefficient, null);
  close(invalidEuler.viscousDragCoefficient, .004); assert.ok(invalidEuler.errors.eulerExit);
});

test('raw-result and iterate helpers agree, use actual Mach and reject mismatched physical metadata', () => {
  const input = fixture(), raw = { conditions: { mach: input.mach, blThermodynamics: 'historical-common-isentrope' },
    mach: .6, checkpoint: input.checkpoint, solverInput: input.checkpoint.restart.input, flow: input.flow,
    kernelReynolds: input.checkpoint.restart.options.reynolds, boundaryLayer: input.bl,
    solverLength: input.solverLength, referenceChord: input.referenceChord };
  assert.deepEqual(quadCoupledTransonicCoefficientsFromResult(raw, { momentReference: input.momentReference }),
    quadCoupledTransonicCoefficients(input));
  assert.throws(() => quadCoupledTransonicCoefficients({ ...input, mach: .6 }), /actual checkpoint/);
  assert.throws(() => quadCoupledTransonicCoefficients({ ...input, bl: { ...input.bl, scale: 1 } }), /normalization/);
  const bad = structuredClone(input); delete bad.checkpoint.restart.options.blThermodynamics;
  assert.throws(() => quadCoupledTransonicCoefficients(bad), /explicit hybrid/);
  delete raw.checkpoint; assert.throws(() => quadCoupledTransonicCoefficientsFromResult(raw), /checkpoint/);
});
