import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';
import { streamtubeEdgeVelocity, streamtubeEdgeVelocityTangent } from '../src/euler/streamtube-edge-velocity.js';

test('thin surface-cell velocity derivatives in a full Newton direction match 80-digit arithmetic', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/mrchdu-wall-edge-decimal.json', import.meta.url)));
  for (const { row, input, expected } of fixture.cases) {
    const point = ([x, y]) => ({ x, y }); let residual = input.ue, derivative = input.ueTangent;
    for (const c of input.cells) {
      const cell = linearizeStreamtubeCell({ lower: c.nodes[0].map(point), upper: c.nodes[1].map(point),
        densities: c.densities, massFlow: c.mass, gamma: input.gamma, stagnationEnthalpy: input.h0 });
      const tangent = cell.apply({ lower: c.nodeTangents[0].map(point), upper: c.nodeTangents[1].map(point),
        densities: c.densities.map((rho, i) => rho * c.logDensityTangents[i]), massFlow: c.massTangent });
      residual -= streamtubeEdgeVelocity(cell.value).ue / input.cells.length;
      derivative -= streamtubeEdgeVelocityTangent(cell.value, tangent) / input.cells.length;
    }
    const residualError = Math.abs(residual - Number(expected.residual));
    const derivativeError = Math.abs(derivative - Number(expected.derivative)) / Math.max(1, Math.abs(Number(expected.derivative)));
    assert.ok(residualError < 5e-10); assert.ok(derivativeError < 1e-8);
    t.diagnostic(JSON.stringify({ row, residualError, derivativeError }));
  }
});
