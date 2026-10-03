import test from 'node:test';
import assert from 'node:assert/strict';
import { harmonicStagnationPoint, runHarmonicStagnationStations } from '../scripts/validation/harmonic-stagnation-stations.js';

test('harmonic signed-arc station map preserves exact mass, bounded spacing and near-wall correspondence', () => {
  for (const psi of [0, 1e-14, .0001, .1, 1]) {
    const points = Array.from({ length: 41 }, (_, i) => harmonicStagnationPoint(i / 20 - 1, psi));
    points.forEach((p, i) => {
      assert.ok(Math.abs(2 * p.x * p.y - psi) < 5e-16);
      assert.ok(Math.abs(p.x - p.y - (i / 20 - 1)) < 5e-16);
      if (i) { const d = Math.hypot(p.x - points[i - 1].x, p.y - points[i - 1].y); assert.ok(d >= .05 / Math.sqrt(2) - 1e-14 && d <= .05 + 1e-14); }
    });
    if (psi) { const p = harmonicStagnationPoint(.5, psi); assert.ok(Math.abs(p.x - .5 - p.y) < 3e-16); }
  }
});

test('Giles SLOR refines physical stagnation mass while retaining a regular station coordinate', t => {
  const evidence = [];
  for (const controls of [{ offset: .25 }, { stretch: 0 }, { stretch: 3 }]) {
    const cases = [8, 16].map(nt => runHarmonicStagnationStations({ nt, ...controls }));
    for (const r of cases) {
      assert.equal(r.converged, true, r.reason); assert.equal(r.quality.valid, true);
      assert.equal(r.boundaryMovement, 0); assert.ok(r.stationError < 1e-12);
      assert.ok(r.maximumAdjacentRatio < Math.sqrt(2));
      if (!r.offset) assert.ok(Math.abs(r.nearWallShiftOverGap - 1) < 1e-10);
    }
    assert.ok(cases[1].massError < .65 * cases[0].massError, JSON.stringify(cases));
    evidence.push(...cases);
  }
  t.diagnostic(JSON.stringify(evidence));
});
