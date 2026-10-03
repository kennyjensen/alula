import test from 'node:test';
import assert from 'node:assert/strict';
import { createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('traced boundary-arc stations retain walls and shared cuts, report physical potential, and obey their arc growth bound', () => {
  const input = intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 12, tubes: 5, tubeGrowth: 3, surfaceSpacing: 'cosine', contourPanels: 40 });
  const original = structuredClone(input), controls = { recordPotentialCoordinates: true };
  const base = createPanelStreamtubeGrid(input, controls);
  const candidate = createPanelStreamtubeGrid(input, { ...controls, interiorStationPlacement: 'boundary-arc' });
  assert.deepEqual(input, original);
  assert.deepEqual(candidate.input, base.input);
  assert.ok(candidate.diagnostics.arcTraces > 0);
  assert.ok(candidate.diagnostics.maxStreamfunctionDrift < 2e-8, String(candidate.diagnostics.maxStreamfunctionDrift));
  const { nodes, potentialCoordinates, arcCoordinates, system } = candidate;
  nodes.forEach((rows, g) => {
    const nt = rows[0].length - 1;
    rows.forEach((row, i) => {
      assert.deepEqual(row[0], base.nodes[g][i][0]);
      assert.deepEqual(row[nt], base.nodes[g][i][nt]);
      for (let j = 0; j <= nt; j++) {
        if (i) { assert.ok(arcCoordinates[g][i][j] > arcCoordinates[g][i - 1][j]); assert.ok(potentialCoordinates[g][i][j].x > potentialCoordinates[g][i - 1][j].x); }
        if (g && j === 0 && !system.layout.active(g - 1, i)) assert.deepEqual(row[j], nodes[g - 1][i].at(-1));
      }
    });
    const d = candidate.diagnostics.boundaryArcTransfer[g];
    assert.ok(d.maximumInteriorArcRatio <= d.boundaryMaximumAdjacentRatio * (1 + 1e-10));
  });
  assert.throws(() => createPanelStreamtubeGrid(input, { interiorStationPlacement: 'boundary-arc', interiorInitialization: 'linear' }), /requires traced/);
});
