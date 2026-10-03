// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { streamtubeCellGeometry } from '../src/euler/streamtube-cell.js';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';
import { streamtubeEdgeVelocity, streamtubeEdgeVelocityTangent } from '../src/euler/streamtube-edge-velocity.js';

test('representable thin-channel sections retain their exact widths after a large translation', () => {
  const gap = 2 ** -42;
  const geometries = [0, 1024].map(offset => {
    const lower = [0, .25, .5].map(x => ({ x, y: offset }));
    const upper = lower.map((p, i) => ({ x: p.x, y: offset + (i + 1) * gap }));
    upper.forEach((p, i) => assert.equal(p.y - lower[i].y, (i + 1) * gap, 'The input shape must be retained exactly.'));
    const geometry = streamtubeCellGeometry(lower, upper);
    assert.deepEqual(geometry.sections, [{ x: 0, y: 1.5 * gap }, { x: 0, y: 2.5 * gap }]);
    return geometry;
  });
  assert.deepEqual(geometries[0], geometries[1]);
});

const close = (a, b, limit) => assert.ok(Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b)) < limit, `${a} != ${b}`);
const features = (v, edge) => ({ isentropicResidual: v.isentropicResidual, streamwiseResidual: v.streamwiseResidual,
  lowerPressure: v.interfacePressure.lower, upperPressure: v.interfacePressure.upper,
  pressureCorrection: v.pressureCorrection, pressureCurvature: v.geometry.pressureCurvature,
  area: v.geometry.area, q0: v.states[0].q, q1: v.states[1].q, edgeVelocity: edge });

test('retained thin Euler and wake cells preserve high-precision values and derivatives without a Python runtime', () => {
  const fixture = JSON.parse(fs.readFileSync('tests/fixtures/coupled-thin-cells-reference.json'));
  assert.equal(fixture.referenceSource.precision, 80); assert.equal(fixture.cases.length, 10);
  for (const c of fixture.cases) {
    const local = linearizeStreamtubeCell(c.parameters);
    for (const [key, v] of Object.entries(features(local.value, streamtubeEdgeVelocity(local.value).ue))) close(v, Number(c.reference.value[key]), 1e-12);
    for (const d of c.directions) {
      const applied = local.apply(d.tangent), expected = c.reference.directions.find(v => v.direction === d.direction);
      for (const [key, v] of Object.entries(features(applied, streamtubeEdgeVelocityTangent(local.value, applied))))
        for (const e of expected.derivatives) close(v, Number(e.value[key]), 1e-10);
    }
  }
});
