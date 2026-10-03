import test from 'node:test';
import assert from 'node:assert/strict';
import { auditHarmonicGrid } from '../src/geometry/harmonic-grid-audit.js';
import { potentialGridQuality } from '../src/geometry/potential-plane-grid.js';

test('exact cylinder streamline nodes and a fixed-polygon harmonic reference are different boundary problems', t => {
  const sqrt = (x, y) => ({ x: Math.sqrt(Math.max(0, (Math.hypot(x, y) + x) / 2)), y: Math.sqrt(Math.max(0, (Math.hypot(x, y) - x) / 2)) });
  const cylinder = w => {
    const a = sqrt(w.x - 2, w.y), b = sqrt(w.x + 2, w.y);
    return { x: .5 * (w.x + a.x * b.x - a.y * b.y), y: .5 * (w.y + a.x * b.y + a.y * b.x) };
  };
  const nx = 12, nt = 6, span = 1.4, width = span / nt;
  const psi = p => p.y * (1 - 1 / (p.x * p.x + p.y * p.y));
  const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => cylinder({ x: -3 + 6 * i / nx, y: width * j })));
  assert.equal(potentialGridQuality(nodes).valid, true);
  let exactNodalError = 0, boundaryChordError = 0;
  nodes.forEach(row => row.forEach((p, j) => { exactNodalError = Math.max(exactNodalError, Math.abs(psi(p) - width * j) / width); }));
  for (let i = 0; i < nx; i++) {
    const p = nodes[i][0], q = nodes[i + 1][0];
    boundaryChordError = Math.max(boundaryChordError, Math.abs(psi({ x: .5 * (p.x + q.x), y: .5 * (p.y + q.y) })) / width);
  }
  const reference = auditHarmonicGrid({ nodes, massFlows: Array(nt).fill(width) });
  const etaDiscrepancy = reference.levels.at(-1).maximum.tubeIntervals;
  assert.ok(exactNodalError < 1e-12);
  assert.ok(boundaryChordError > .1, 'Straight chords do not retain the analytic curved-wall boundary data.');
  assert.ok(etaDiscrepancy > .05, 'A polygon-domain error cannot be called an error in these exact nodal streamlines.');
  assert.ok(reference.referenceChange.tubeIntervals < .01);
  // The FE solve is valid for its polygon domain. This counterexample must
  // not be used to certify an inaccurate foil grid or to erase its errors.
  t.diagnostic(JSON.stringify({ exactNodalError, boundaryChordError, etaDiscrepancy,
    lastEtaReferenceChange: reference.referenceChange.tubeIntervals }));
});
