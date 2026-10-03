// SPDX-License-Identifier: GPL-2.0-or-later
// Resample one requested surface distribution between its actual physical
// passage hits. Counts remain shared; the density shape is local to this
// surface. This alone does not reconcile competing facing distributions or
// guarantee compatible derivatives where consecutive hit blocks meet.
import { createMonotoneCubicMap } from '../numerics/monotone-cubic.js';
import { createVinokurStations } from '../numerics/vinokur-stations.js';

export function distributeSurfaceByDensity({ requested, anchors, endpointSpacings }) {
  if (!Array.isArray(requested) || requested.length < 3
    || requested.some((s, i) => !Number.isFinite(s) || i && !(s > requested[i - 1]))
    || !Array.isArray(anchors) || anchors.length < 2
    || anchors.some((a, i) => !a || !Number.isInteger(a.index) || a.index < 0 || !Number.isFinite(a.position)
      || i && (!(a.index > anchors[i - 1].index) || !(a.position > anchors[i - 1].position)))
    || anchors[0].position !== requested[0] || anchors.at(-1).position !== requested.at(-1))
    throw new Error('Ordered requested surface positions and endpoint-preserving hit indices are required.');
  if (endpointSpacings !== undefined && (!Array.isArray(endpointSpacings) || endpointSpacings.length !== anchors.length - 1
    || endpointSpacings.some((e, k) => !e || ![e.firstSpacing, e.lastSpacing].every(h => Number.isFinite(h) && h > 0)
      || !(e.firstSpacing < anchors[k + 1].position - anchors[k].position - e.lastSpacing))))
    throw new Error('One feasible pair of endpoint spacings is required for every surface block.');
  const labels = requested.map((_, i) => i), densityMap = createMonotoneCubicMap(labels, requested, { derivatives: 'pchip' });
  const inverse = s => {
    const exact = requested.indexOf(s);
    if (exact >= 0) return exact;
    let lo = 0, hi = labels.at(-1);
    for (let k = 0; k < 64; k++) {
      const mid = lo + .5 * (hi - lo);
      if (mid === lo || mid === hi) break;
      if (densityMap.value(mid) < s) lo = mid; else hi = mid;
    }
    return lo + .5 * (hi - lo);
  };
  const hits = anchors.map(a => ({ ...a, densityCoordinate: inverse(a.position) }));
  const positions = [anchors[0].position], coordinates = [hits[0].densityCoordinate], blocks = [];
  for (let k = 1; k < hits.length; k++) {
    const a = hits[k - 1], b = hits[k], count = b.index - a.index;
    const dq = (b.densityCoordinate - a.densityCoordinate) / count;
    if (!(dq > 0 && Number.isFinite(dq))) throw new Error('Unresolved surface density per interval.');
    const ends = endpointSpacings?.[k - 1];
    const stretch = ends && createVinokurStations({ length: b.densityCoordinate - a.densityCoordinate, intervals: count,
      firstSpacing: inverse(a.position + ends.firstSpacing) - a.densityCoordinate,
      lastSpacing: b.densityCoordinate - inverse(b.position - ends.lastSpacing) });
    const start = positions.length - 1;
    for (let i = 1; i <= count; i++) {
      const q = i === count ? b.densityCoordinate : a.densityCoordinate + (stretch ? stretch.positions[i] : i * dq);
      coordinates.push(q); positions.push(i === count ? b.position : densityMap.value(q));
    }
    blocks.push({ fromIndex: a.index, toIndex: b.index, intervals: count,
      requestedDensityIntervals: b.densityCoordinate - a.densityCoordinate,
      densityCoordinatePerInterval: dq,
      ...(stretch ? { endpointFit: { ...ends, branch: stretch.branch, parameter: stretch.parameter,
        logAsymmetry: stretch.logAsymmetry, metric: 'requested-density coordinate',
        uniformDensityIncrements: false } } : {}),
      firstSpacing: positions[start + 1] - positions[start], lastSpacing: positions.at(-1) - positions.at(-2) });
  }
  if (positions.some((s, i) => !Number.isFinite(s) || i && !(s > positions[i - 1])))
    throw new Error('Resampled surface stations are not resolved and increasing.');
  return { positions, coordinates, firstIndex: anchors[0].index, anchors: hits, blocks,
    interpolation: 'PCHIP of requested physical positions against requested node index',
    placement: endpointSpacings ? 'two-ended Vinokur distribution in requested-density coordinate within each actual hit block'
      : 'uniform requested-density increments within each actual hit block',
    countPolicy: 'retain supplied hit indices; report achieved/requested density scale',
    facingReconciled: false, joinDerivativeMatched: false };
}
