// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import { linearizeStreamtubeCell } from '../../src/euler/streamtube-linearization.js';
import { streamtubeEdgeVelocityTangent } from '../../src/euler/streamtube-edge-velocity.js';

const oracle = cases => {
  const r = spawnSync('python3', ['scripts/validation/decimal-streamtube-cell.py'], {
    input: JSON.stringify({ cases }), encoding: 'utf8', maxBuffer: 4 * 1024 * 1024 });
  assert.equal(r.status, 0, r.stderr); return JSON.parse(r.stdout);
};
const error = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));
const close = (a, b, tolerance = 1e-12) => assert.ok(error(Number(a), Number(b)) < tolerance, `${a} != ${b}`);
const uniform = { lower: [{ x: 0, y: 0 }, { x: .5, y: 0 }, { x: 1, y: 0 }],
  upper: [{ x: 0, y: .25 }, { x: .5, y: .25 }, { x: 1, y: .25 }],
  densities: [1, 1], massFlow: .25, stagnationEnthalpy: 8, gamma: 1.5, pressureCorrectionFactor: .125 };
const features = (local, v) => ({ isentropicResidual: v.isentropicResidual, streamwiseResidual: v.streamwiseResidual,
  lowerPressure: v.interfacePressure.lower, upperPressure: v.interfacePressure.upper,
  pressureCorrection: v.pressureCorrection, pressureCurvature: v.geometry.pressureCurvature,
  area: v.geometry.area, q0: v.states[0].q, q1: v.states[1].q,
  edgeVelocity: streamtubeEdgeVelocityTangent(local.value, v) });

test('80-digit oracle recovers exact uniform-flow values and independently known mass derivatives', () => {
  const r = oracle([{ id: 'uniform', parameters: uniform, directions: [{ direction: 0, tangent: { massFlow: .25 } }] }]);
  assert.equal(r.precision, 80); assert.equal(r.inputConversion, 'exact binary64');
  const c = r.cases[0];
  for (const key of ['isentropicResidual', 'streamwiseResidual', 'pressureCorrection', 'pressureCurvature']) close(c.value[key], 0);
  close(c.value.area, .125); close(c.value.lowerPressure, 2.5); close(c.value.upperPressure, 2.5);
  close(c.value.edgeVelocity, 1);
  for (const estimate of c.directions[0].derivatives) {
    close(estimate.value.lowerPressure, -1 / 3); close(estimate.value.upperPressure, -1 / 3);
    close(estimate.value.edgeVelocity, 1); close(estimate.value.q0, 1); close(estimate.value.q1, 1);
    close(estimate.value.isentropicResidual, 0); close(estimate.value.streamwiseResidual, 0);
  }
});

test('independent tiny-step reference resolves all local coordinates and flow derivatives on both smooth Mach branches', () => {
  const p = { lower: [{ x: 0, y: 0 }, { x: .45, y: .025 }, { x: 1.1, y: .08 }],
    upper: [{ x: .1, y: .25 }, { x: .58, y: .3 }, { x: 1.18, y: .4 }],
    densities: [1.05, 1.02], massFlow: .27, stagnationEnthalpy: 8, gamma: 1.4, pressureCorrectionFactor: .1 };
  const tangents = [];
  for (const side of ['lower', 'upper']) for (let i = 0; i < 3; i++) for (const key of ['x', 'y']) {
    const row = Array.from({ length: 3 }, () => ({ x: 0, y: 0 })); row[i][key] = 1; tangents.push({ [side]: row });
  }
  tangents.push({ densities: [1, 0] }, { densities: [0, 1] }, { massFlow: 1 }, { stagnationEnthalpy: 1 }, { pressureCorrectionFactor: 1 });
  const cases = [.27, .65].map(massFlow => ({ id: String(massFlow), parameters: { ...p, massFlow },
    directions: tangents.map((tangent, direction) => ({ direction, tangent })) }));
  const result = oracle(cases);
  cases.forEach((c, i) => {
    const local = linearizeStreamtubeCell(c.parameters);
    c.directions.forEach((d, j) => {
      const v = local.apply(d.tangent), expected = result.cases[i].directions[j].derivatives;
      const values = features(local, v);
      for (const [key, value] of Object.entries(values)) {
        close(expected[0].value[key], expected[1].value[key], 1e-15);
        for (const e of expected) close(value, e.value[key]);
      }
      // A meaningful negative control: a missing density/mass term cannot
      // pass merely because both finite-difference steps agree.
      if (d.tangent.massFlow) assert.ok(error(0, Number(expected[0].value.edgeVelocity)) > .1);
    });
  });
});

test('retained leading-edge, trailing-edge and wake thin-cell tangents agree with independent high precision', () => {
  const { cases } = JSON.parse(fs.readFileSync('tests/fixtures/coupled-thin-cells.json'));
  assert.equal(cases.length, 10);
  const result = oracle(cases);
  cases.forEach((c, i) => {
    const local = linearizeStreamtubeCell(c.parameters);
    c.directions.forEach((d, j) => {
      const values = features(local, local.apply(d.tangent)), estimates = result.cases[i].directions[j].derivatives;
      for (const [key, value] of Object.entries(values)) {
        close(estimates[0].value[key], estimates[1].value[key], 1e-15);
        for (const estimate of estimates) close(value, estimate.value[key], 1e-7);
      }
    });
  });
});
