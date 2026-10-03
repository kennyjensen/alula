// SPDX-License-Identifier: GPL-2.0-or-later
// Reconcile local spacing requests in normalized physical-arc progress.
// This explicit reconstruction reports a count deficit instead of implying
// that a fixed count can satisfy every facing request.
import { createSpacingEnvelope } from './spacing-envelope.js';
import { createVinokurStations } from '../numerics/vinokur-stations.js';

export function distributeFacingDensity({ profiles, ranges, intervals, maximumIntervals = 10000, endpointSpacings }) {
  if (!Array.isArray(profiles) || !profiles.length || !Number.isInteger(maximumIntervals) || maximumIntervals < 2
    || intervals !== undefined && (!Number.isInteger(intervals) || intervals < 2 || intervals > maximumIntervals)
    || profiles.some(p => !Array.isArray(p) || p.length < 3 || p.some((s, i) => !Number.isFinite(s) || i && !(s > p[i - 1]))))
    throw new Error('Facing density requires ordered physical station requests and valid interval counts.');
  ranges ??= profiles.map(p => [p[0], p.at(-1)]);
  if (!Array.isArray(ranges) || ranges.length !== profiles.length || ranges.some((r, k) => !Array.isArray(r) || r.length !== 2
    || !r.every(Number.isFinite) || !(r[0] < r[1]) || r[0] < profiles[k][0] || r[1] > profiles[k].at(-1)))
    throw new Error('Facing ranges must lie inside their original requests.');
  const lengths = ranges.map(([a, b]) => b - a);
  const normalized = profiles.map((p, k) => p.map(s => s === ranges[k][0] ? 0 : s === ranges[k][1] ? 1 : (s - ranges[k][0]) / lengths[k]));
  const envelope = createSpacingEnvelope({ profiles: normalized, start: 0, end: 1 });
  const requestedIntervals = envelope.total, roundoff = 64 * Number.EPSILON * Math.max(1, requestedIntervals);
  const requiredIntervals = Math.max(2, Math.ceil(requestedIntervals - roundoff));
  const count = intervals ?? requiredIntervals;
  if (count > maximumIntervals) throw new Error('Facing density exceeds the supplied station budget.');
  let stretch;
  if (endpointSpacings) {
    const { firstSpacing, lastSpacing } = endpointSpacings;
    if (![firstSpacing, lastSpacing].every(v => Number.isFinite(v) && v > 0) || !(firstSpacing + lastSpacing < 1))
      throw new Error('Facing endpoint spacings must be feasible normalized physical intervals.');
    stretch = createVinokurStations({ length: requestedIntervals, intervals: count,
      firstSpacing: envelope.metric(firstSpacing), lastSpacing: requestedIntervals - envelope.metric(1 - lastSpacing) });
  }
  const coordinates = Array.from({ length: count + 1 }, (_, i) => stretch ? stretch.positions[i] : requestedIntervals * i / count);
  const progress = coordinates.map((q, i) => !i ? 0 : i === count ? 1 : envelope.inverse(q));
  if (progress.some((u, i) => !Number.isFinite(u) || i && !(u > progress[i - 1]))) throw new Error('Unresolved facing-density stations.');
  const maximumDensityInterval = Math.max(...coordinates.slice(1).map((q, i) => q - coordinates[i]));
  return { progress, intervals: count, requestedIntervals, requiredIntervals,
    densityScale: count / requestedIntervals, countCapacitySatisfied: count + roundoff >= requestedIntervals,
    resolutionSatisfied: maximumDensityInterval <= 1 + roundoff, maximumDensityInterval,
    ...(stretch ? { endpointFit: { ...endpointSpacings, branch: stretch.branch, parameter: stretch.parameter,
      metric: 'common facing-density integral' } } : {}),
    positions: ranges.map(([a, b], k) => progress.map((u, i) => !i ? a : i === count ? b : a + lengths[k] * u)),
    metric: 'lower envelope of positive interval requests in normalized physical arc',
    exactModernMsetLaw: false, endpointJoinMatched: false };
}
