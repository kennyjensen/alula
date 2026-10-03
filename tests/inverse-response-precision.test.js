import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';
import { streamtubeEdgeVelocity, streamtubeEdgeVelocityTangent } from '../src/euler/streamtube-edge-velocity.js';

test('retained surface and wake stencils match independent 80-digit values and directional derivatives', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/inverse-response-edge-decimal.json', import.meta.url)));
  const sha = value => createHash('sha256').update(value).digest('hex');
  assert.equal(sha(fs.readFileSync(new URL('../' + fixture.oracle.path, import.meta.url))), fixture.oracle.sha256);
  let valueError = 0, derivativeError = 0;
  for (const c of fixture.cases) {
    const cell = linearizeStreamtubeCell(c.parameters), v = streamtubeEdgeVelocity(cell.value).ue;
    valueError = Math.max(valueError, Math.abs(v - Number(c.expected.value.edgeVelocity)));
    for (const [i, direction] of c.directions.entries()) {
      const derivative = streamtubeEdgeVelocityTangent(cell.value, cell.apply(direction.tangent));
      for (const e of c.expected.directions[i].derivatives) derivativeError = Math.max(derivativeError,
        Math.abs(derivative - Number(e.value.edgeVelocity)) / Math.max(1, Math.abs(derivative), Math.abs(Number(e.value.edgeVelocity))));
    }
  }
  assert.equal(fixture.cases.length, 3); assert.ok(valueError < 5e-10 && derivativeError < 1e-8);
  t.diagnostic(JSON.stringify({ cells: fixture.cases.length, valueError, derivativeError }));
});
