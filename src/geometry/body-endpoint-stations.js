// SPDX-License-Identifier: GPL-2.0-or-later
// Joint wall/cut endpoint lengths, retaining counts, actual hits and natural
// intervals at internal block joins. Reconstructed spacing policy, not MSET
// source or a physical continuity condition at stagnation.
import { distributeSurfaceByDensity } from './surface-density-stations.js';

export function reconcileBodyEndpointStations({ upper, lower, upstream, wake }) {
  const components = { upper, lower, upstream, wake }, prepared = {};
  for (const [name, c] of Object.entries(components)) {
    if (!c || !Array.isArray(c.positions) || c.positions.length < 4 || !Array.isArray(c.anchors) || c.anchors.length < 2
      || c.positions.some((s, i) => !Number.isFinite(s) || i && !(s > c.positions[i - 1])))
      throw new Error('Body endpoint reconciliation needs ordered physical positions and actual hit anchors.');
    const firstIndex = c.anchors[0].index;
    if (c.anchors.some(a => c.positions[a.index - firstIndex] !== a.position)
      || c.anchors.at(-1).index - firstIndex !== c.positions.length - 1)
      throw new Error('Body endpoint anchors must coincide exactly with supplied stations.');
    prepared[name] = { ...c, endpointSpacings: c.anchors.slice(1).map((a, k) => {
      const i = c.anchors[k].index - firstIndex, j = a.index - firstIndex;
      return { firstSpacing: c.positions[i + 1] - c.positions[i], lastSpacing: c.positions[j] - c.positions[j - 1] };
    }) };
  }
  const joins = [
    { edge: 'LE', occurrences: [['upper', 0, 'firstSpacing'], ['lower', 0, 'firstSpacing'], ['upstream', -1, 'lastSpacing']] },
    { edge: 'TE', occurrences: [['upper', -1, 'lastSpacing'], ['lower', -1, 'lastSpacing'], ['wake', 0, 'firstSpacing']] },
  ].map(({ edge, occurrences }) => {
    const natural = occurrences.map(([name, k, key]) => prepared[name].endpointSpacings.at(k)[key]);
    const target = Math.exp(natural.reduce((sum, h) => sum + Math.log(h), 0) / natural.length);
    occurrences.forEach(([name, k, key]) => { prepared[name].endpointSpacings.at(k)[key] = target; });
    return { edge, natural, target };
  });
  const fitted = Object.fromEntries(Object.entries(prepared).map(([name, c]) => {
    try { return [name, distributeSurfaceByDensity({ requested: c.requested ?? c.positions, anchors: c.anchors, endpointSpacings: c.endpointSpacings })]; }
    catch (error) { throw new Error(`${name} endpoint fit: ${error.message}`, { cause: error }); }
  }));
  return { ...fitted, joins, policy: 'minimize sum of squared log changes in the three incident endpoint lengths',
    internalJoinIntervalsRetained: true, facingReconciled: false, joinDerivativeMatched: false };
}
