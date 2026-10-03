// SPDX-License-Identifier: GPL-2.0-or-later
// Common potential-based cross-lines. Surface spacing demands are combined
// by their lower envelope, not by unioning nearly coincident node positions.
// This supplies compatible counts at every LE/TE block boundary.
import { createSpacingEnvelope } from './spacing-envelope.js';
const blend = t => t * t * (3 - 2 * t);

export function createBodyPotentialMap({ start, end, leading, upperTrailing, lowerTrailing, inlet, outletIncrement }) {
  const trailing = .5 * (upperTrailing + lowerTrailing);
  if (![start, end, leading, upperTrailing, lowerTrailing, inlet, outletIncrement].every(Number.isFinite)
    || !(start < leading && leading < trailing && trailing < end)
    || inlet >= leading || Math.min(upperTrailing, lowerTrailing) <= leading || outletIncrement <= 0)
    throw new Error('Invalid potential block anchors.');
  const offsets = { upper: upperTrailing - trailing, lower: lowerTrailing - trailing };
  const inletOffset = inlet - start, outletOffset = trailing + outletIncrement - end;
  // The smooth cubic transitions have |B'| <= 1.5. Certify strict
  // monotonicity before using inverse maps; never sort a folded mapping.
  const minimumDerivative = Math.min(1 - 1.5 * Math.max(0, inletOffset) / (leading - start),
    ...Object.values(offsets).map(d => 1 + 1.5 * Math.min(0, d) / (trailing - leading)),
    1 + 1.5 * Math.min(0, outletOffset) / (end - trailing));
  if (!(minimumDerivative > 0)) throw new Error('Potential cross-line map is not monotone; revise the block topology.');
  const value = (side, chi) => {
    if (!Object.hasOwn(offsets, side) || !Number.isFinite(chi) || chi < start || chi > end) throw new Error('Potential cross-line coordinate is outside its block map.');
    if (chi <= leading) return chi + inletOffset * blend((leading - chi) / (leading - start));
    if (chi <= trailing) return chi + offsets[side] * blend((chi - leading) / (trailing - leading));
    return chi + offsets[side] + outletOffset * blend((chi - trailing) / (end - trailing));
  };
  const inverse = (side, phi) => {
    if (!Number.isFinite(phi) || phi < value(side, start) || phi > value(side, end)) throw new Error('Requested potential lies outside its block map.');
    let lo = start, hi = end;
    for (let i = 0; i < 60; i++) { const mid = .5 * (lo + hi); if (value(side, mid) < phi) lo = mid; else hi = mid; }
    return .5 * (lo + hi);
  };
  return { start, end, leading, trailing, offsets, minimumDerivative, value, inverse };
}

export function createPotentialCrosslines({ profiles, start, end, growth = .25, maxIntervals = 2000, inletIntervals, outletIntervals, blockIntervals }) {
  if (!Array.isArray(profiles) || !profiles.length || ![start, end, growth].every(Number.isFinite)
    || !(start < end) || growth <= 0 || !Number.isInteger(maxIntervals) || maxIntervals < 4
    || [inletIntervals, outletIntervals].some(n => n !== undefined && (!Number.isInteger(n) || n < 2))
    || profiles.some(row => !Array.isArray(row) || row.length < 3 || !row.every(Number.isFinite)
      || row[0] <= start || row.at(-1) >= end || row.some((v, i) => i && v <= row[i - 1])))
    throw new Error('Invalid potential spacing profiles.');
  const { total, anchors, metric, inverse } = createSpacingEnvelope({ profiles, start, end, growth });
  if (blockIntervals !== undefined && (!Array.isArray(blockIntervals) || blockIntervals.length !== anchors.length - 1
    || blockIntervals.some(n => !Number.isInteger(n) || n < 2))) throw new Error('One valid interval count is required per potential block.');
  const x = [start], blocks = [];
  for (let i = 1; i < anchors.length; i++) {
    const a = anchors[i - 1], b = anchors[i], ma = metric(a), mb = metric(b), length = mb - ma;
    const requested = blockIntervals?.[i - 1] ?? (i === 1 ? inletIntervals : i === anchors.length - 1 ? outletIntervals : undefined);
    const count = requested ?? Math.max(2, Math.ceil(length - 64 * Number.EPSILON * Math.max(1, length)));
    if (x.length - 1 + count > maxIntervals) throw new Error('Potential cross-line budget exceeded.');
    blocks.push({ start: a, end: b, intervals: count, metricLength: length });
    for (let k = 1; k <= count; k++) x.push(k === count ? b : inverse(ma + length * k / count));
  }
  if (x.some((v, i) => !Number.isFinite(v) || (i && v <= x[i - 1]))) throw new Error('Nonmonotone potential cross-lines.');
  return { x, anchors, blocks, metricLength: total };
}
