// SPDX-License-Identifier: GPL-2.0-or-later
import { intrinsicBodyFixture } from '../../tests/fixtures/intrinsic-body.js';
import { solveCoupledStreamtubeIses } from '../../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';
import { truncateCoupledDownstreamDomain } from '../../src/euler/tests/streamtube-coupled-truncation.js';

export function coupledDownstreamControl() {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    wakeGeometry: 'independent-banks', wakeOutlet: 'banks', streamwiseMode: 'isentropic' };
  const root = solveCoupledStreamtubeIses(input, { reynolds: 1e5, edgeMatching: 'section-velocity', tolerance: 1e-10, stepAcceptance: 'admissible' });
  if (!root.converged) throw new Error(root.reason);
  const ne = root.x.length - 4 * root.boundaryLayer.stations.length;
  const source = createCoupledStreamtubeBody(root.solverInput, { ...root.coupledOptions, initialEuler: root.flow, initialBL: root.x.slice(ne) });
  const mapped = truncateCoupledDownstreamDomain(root.solverInput, source, { endIndex: source.euler.layout.nx - 1 });
  const seed = mapped.system.evaluate(mapped.system.initial);
  const result = solveCoupledStreamtubeBody(mapped.system, { tolerance: 1e-10, maxIterations: 12 });
  return { input: root.solverInput, source, mapped, seed, result };
}
