// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createUpwindStreamtubeChannel } from '../src/euler/tests/streamtube-upwind-channel.js';
import { initializeUpwindStreamtubeChannel } from '../src/euler/tests/streamtube-channel-startup.js';

const close = (a, b, label = '', tolerance = 2e-12) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${label}: ${a} != ${b}`);
function inclined(gamma = 1.4) {
  const x = [0, .3, .8, 1.4, 2], slope = .25, height = .24;
  return { x, lower: x.map(x => .03 + slope * x), upper: x.map(x => .03 + height + slope * x),
    massFlows: [.06, .14], stagnationDensity: [1.2, 1.05], stagnationEnthalpy: [4, 5],
    gamma, referenceDensity: .9, referencePressure: 1.2, outletPressure: 1.1,
    inletSlopes: [slope], upwind: { mucon: 1, mcrit: .6, boundary: { kind: 'unfiltered-first-two' } } };
}
function contraction() {
  const p = inclined(), lower = [0, .015, .04, .03, .02], heights = [.3, .26, .18, .24, .3];
  return { ...p, lower, upper: lower.map((v, i) => v + heights[i]), inletSlopes: [0] };
}
// Direct section-face projection, independent of streamtubeCellGeometry.
function normalAreas(nodes, tube) {
  return nodes.slice(1).map((row, i) => {
    const previous = nodes[i];
    const dx = .5 * (row[tube].x + row[tube + 1].x - previous[tube].x - previous[tube + 1].x);
    const dy = .5 * (row[tube].y + row[tube + 1].y - previous[tube].y - previous[tube + 1].y);
    const gx = .5 * (row[tube + 1].x - row[tube].x + previous[tube + 1].x - previous[tube].x);
    const gy = .5 * (row[tube + 1].y - row[tube].y + previous[tube + 1].y - previous[tube].y);
    return (dx * gy - dy * gx) / Math.hypot(dx, dy);
  });
}

test('startup uses independent closed-form sonic capacity and uniform reservoir density', () => {
  for (const gamma of [1.4, 5 / 3]) {
    const p = inclined(gamma), system = createUpwindStreamtubeChannel(p);
    const original = { initial: system.initial.slice(), conditions: structuredClone(system.conditions), input: structuredClone(p) };
    const result = initializeUpwindStreamtubeChannel(system), { diagnostics: d, flow, initial } = result;
    for (let j = 0; j < system.nt; j++) {
      // gamma=7/5: rho*/rho0=(5/6)^(5/2), q*²=h0/3.
      // gamma=5/3: rho*/rho0=(3/4)^(3/2), q*²=h0/2.
      const densityRatio = gamma === 1.4 ? (5 / 6) ** (5 / 2) : (3 / 4) ** (3 / 2);
      const sonicDensity = p.stagnationDensity[j] * densityRatio;
      const sonicSpeed = Math.sqrt(p.stagnationEnthalpy[j] / (gamma === 1.4 ? 3 : 2));
      const area = .24 * (j === 0 ? .3 : .7) / Math.sqrt(1 + .25 ** 2);
      const mass = area * sonicDensity * sonicSpeed;
      close(d.tubes[j].sonicDensity, sonicDensity);
      close(d.tubes[j].sonicSpeed, sonicSpeed);
      close(d.tubes[j].minimumNormalArea, area);
      close(d.capacityMassFlows[j], mass);
      close(flow.massFlows[j], mass);
      close(initial[system.massIndex(j)], Math.log(mass / p.massFlows[j]));
      for (let i = 0; i < system.nx; i++) {
        close(initial[i * system.nt + j], Math.log(p.stagnationDensity[j] / p.referenceDensity));
        close(flow.sections[i][j].rho, p.stagnationDensity[j]);
        close(flow.sections[i][j].q, densityRatio * sonicSpeed);
      }
    }
    assert.equal(d.initialGuessOnly, true); assert.equal(d.targetEquationsUnchanged, true);
    assert.equal(d.massFlowsAreUnknowns, true);
    assert.ok(d.initialMaxMach > 0 && d.initialMaxMach < 1);
    assert.ok(d.residual > 0, 'Reservoir-density startup is an unsolved initial guess.');
    close(d.residual, Math.max(...flow.residual.map(Math.abs)));
    assert.deepEqual(system.residual(initial), flow.residual);
    assert.deepEqual(system.initial, original.initial);
    assert.deepEqual(system.conditions, original.conditions);
    assert.deepEqual(p, original.input);
  }
});

test('startup chooses the smallest signed physical section capacity on a curved contraction', () => {
  const p = contraction(), system = createUpwindStreamtubeChannel(p), result = initializeUpwindStreamtubeChannel(system);
  const { nodes } = system.decode(system.initial);
  for (let j = 0; j < system.nt; j++) {
    const areas = normalAreas(nodes, j), minimum = Math.min(...areas), tube = result.diagnostics.tubes[j];
    assert.ok(areas.every(v => v > 0));
    close(tube.minimumNormalArea, minimum);
    assert.equal(tube.minimumAreaSection, areas.indexOf(minimum));
    close(result.flow.massFlows[j], minimum * p.stagnationDensity[j] * (5 / 6) ** (5 / 2) * Math.sqrt(p.stagnationEnthalpy[j] / 3));
  }
  assert.ok(result.flow.residual.every(Number.isFinite));
  assert.ok(result.flow.sections[0].every(s => s.machSquared < 1));
  assert.ok(result.flow.sections.at(-1).every(s => s.machSquared < 1));
});

test('arbitrary reference-mass scaling including an inadmissible 100x guess does not determine the physical seed', () => {
  const p = contraction(), baseline = initializeUpwindStreamtubeChannel(createUpwindStreamtubeChannel(p));
  for (const factor of [.01, 100]) {
    const system = createUpwindStreamtubeChannel({ ...p, massFlows: p.massFlows.map(m => factor * m) });
    if (factor === 100) assert.throws(() => system.evaluate(system.initial));
    const result = initializeUpwindStreamtubeChannel(system);
    result.flow.massFlows.forEach((m, j) => close(m, baseline.flow.massFlows[j]));
    result.flow.sections.forEach((row, i) => row.forEach((section, j) => {
      for (const key of ['rho', 'q', 'p', 'enthalpy', 'machSquared']) close(section[key], baseline.flow.sections[i][j][key], key);
    }));
    for (let j = 0; j < system.nt; j++) close(result.initial[system.massIndex(j)],
      baseline.initial[system.massIndex(j)] - Math.log(factor));
  }
});

test('startup needs no back-pressure or shock-location guess and preserves geometric unit transformations', () => {
  const p = contraction(), baseline = initializeUpwindStreamtubeChannel(createUpwindStreamtubeChannel(p));
  for (const outletPressure of [.7, 1.4]) {
    const result = initializeUpwindStreamtubeChannel(createUpwindStreamtubeChannel({ ...p, outletPressure }));
    assert.deepEqual(result.initial, baseline.initial);
    assert.deepEqual(result.flow.massFlows, baseline.flow.massFlows);
    assert.notDeepEqual(result.flow.residual, baseline.flow.residual);
  }
  for (const scale of [.2, 7]) {
    const transformed = { ...p, x: p.x.map(x => 3 + scale * x), lower: p.lower.map(y => -2 + scale * y),
      upper: p.upper.map(y => -2 + scale * y), massFlows: p.massFlows.map(m => scale * m) };
    const result = initializeUpwindStreamtubeChannel(createUpwindStreamtubeChannel(transformed));
    result.initial.forEach((v, i) => close(v, baseline.initial[i]));
    result.flow.massFlows.forEach((m, j) => close(m / scale, baseline.flow.massFlows[j]));
    result.flow.sections.forEach((row, i) => row.forEach((s, j) => {
      for (const key of ['rho', 'q', 'p', 'enthalpy', 'machSquared']) close(s[key], baseline.flow.sections[i][j][key], key);
    }));
    close(result.diagnostics.minimumNormalArea / scale, baseline.diagnostics.minimumNormalArea);
    close(result.diagnostics.initialMaxMach, baseline.diagnostics.initialMaxMach);
  }
});
