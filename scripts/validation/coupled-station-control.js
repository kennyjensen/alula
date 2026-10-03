// SPDX-License-Identifier: GPL-2.0-or-later
// Small complete multielement control shared by Node and browser verification.
import { intrinsicBodyFixture } from '../../tests/fixtures/intrinsic-body.js';
import { solveCoupledStreamtubeIses } from '../../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';
import { redistributeCoupledSurfaceStations } from '../../src/euler/tests/streamtube-coupled-redistribution.js';

export function coupledStationControl() {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    wakeGeometry: 'independent-banks', streamwiseMode: 'isentropic' };
  const root = solveCoupledStreamtubeIses(input, { reynolds: 1e5, edgeMatching: 'section-velocity', tolerance: 1e-10, stepAcceptance: 'admissible' });
  if (!root.converged) throw new Error(root.reason);
  const ne = root.x.length - 4 * root.boundaryLayer.stations.length;
  const source = createCoupledStreamtubeBody(root.solverInput, { ...root.coupledOptions, initialEuler: root.flow, initialBL: root.x.slice(ne) });
  const surfaceFractions = source.euler.fractions.map((body, b) => Object.fromEntries(['upper', 'lower'].map(side =>
    [side, body[side].map(f => f === 0 || f === 1 ? f : f + (side === 'upper' ? 1 : -1) * (b + 1) * .01 * Math.sin(Math.PI * f) ** 3)])));
  const mapped = redistributeCoupledSurfaceStations(root.solverInput, source, { surfaceFractions });
  // Newton rebases its geometry chart; retain the physical seed before it runs.
  const seed = mapped.system.evaluate(mapped.system.initial);
  const result = solveCoupledStreamtubeBody(mapped.system, { tolerance: 1e-10, maxIterations: 12 });
  return { input: root.solverInput, source, mapped, seed, result };
}
