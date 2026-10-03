// SPDX-License-Identifier: GPL-2.0-or-later
// Reparameterize prescribed outer streamlines through fixed block targets.
// Interior blocks use standard SLATEC/PCHIP interpolation. Optional Vinokur
// end segments match the interior join slope and have zero far-end curvature.
// This is boundary preparation, not modern MSET's interior X-spacing force.
import { createMonotoneCubicMap } from '../numerics/monotone-cubic.js';
import { createVinokurStretch } from '../numerics/vinokur-stretch.js';

function withVinokurEnds(base) {
  const { knots: x, values: y } = base, last = x.length - 1;
  // A single anchored interval is already affine, the zero-parameter case.
  if (last === 1) return { ...base, endSegments: [] };
  const segments = [0, last - 1].map((i, side) => {
    const width = x[i + 1] - x[i], length = y[i + 1] - y[i];
    const stretch = createVinokurStretch({ length, logicalWidth: width,
      initialDerivative: base.slopes[side === 0 ? 1 : last - 1] });
    const evaluate = station => {
      const t = side === 0 ? (x[i + 1] - station) / width : (station - x[i]) / width;
      const p = stretch.evaluate(t);
      return { value: station === x[i] ? y[i] : station === x[i + 1] ? y[i + 1]
        : side === 0 ? y[i + 1] - p.value : y[i] + p.value,
      derivative: p.derivative, secondDerivative: side === 0 ? -p.secondDerivative : p.secondDerivative };
    };
    return { evaluate, record: { interval: i, branch: stretch.branch, parameter: stretch.parameter,
      logicalWidth: width, length, orientation: side === 0 ? 'reflected inlet' : 'outlet',
      left: evaluate(x[i]), right: evaluate(x[i + 1]) } };
  });
  const evaluate = station => {
    if (!Number.isFinite(station) || station < x[0] || station > x[last])
      throw new Error('Outer spacing argument lies outside its knots.');
    if (station < x[1]) return segments[0].evaluate(station);
    if (station >= x[last - 1]) return segments[1].evaluate(station);
    return base.evaluate(station);
  };
  const slopes = base.slopes.slice();
  slopes[0] = segments[0].record.left.derivative;
  slopes[last] = segments[1].record.right.derivative;
  return { ...base, slopes, evaluate, value: x => evaluate(x).value,
    endSegments: segments.map(s => s.record) };
}

export function createOuterStreamtubeSpacing({ ranks, anchors, targets, spread = .85,
  derivatives = 'pchip', endSpacing = 'cubic' }) {
  const increasing = values => Array.isArray(values) && values.length >= 2
    && Array.from(values).every((v, i) => Number.isFinite(v) && (!i || v > values[i - 1]));
  if (!increasing(ranks) || !increasing(anchors) || !Array.isArray(targets)
    || !increasing(targets.map(p => p.rank)) || !increasing(targets.map(p => p.potential))
    || !Number.isFinite(spread) || spread < 0 || spread > 1
    || !['pchip', 'harmonic'].includes(derivatives) || !['cubic', 'vinokur'].includes(endSpacing))
    throw new Error('Invalid outer streamtube spacing data.');
  ranks = ranks.slice(); anchors = anchors.slice(); targets = targets.map(p => ({ ...p }));
  if (anchors[0] !== ranks[0] || anchors.at(-1) !== ranks.at(-1)
    || targets[0].rank !== ranks[0] || targets.at(-1).rank !== ranks.at(-1)
    || targets.some(p => !anchors.includes(p.rank)))
    throw new Error('Outer spacing targets must lie on its anchored domain.');
  const indices = anchors.map(rank => ranks.indexOf(rank));
  if (indices.some(i => i < 0)) throw new Error('Every outer block anchor needs a station.');
  const rankMap = createMonotoneCubicMap(targets.map(p => p.rank), targets.map(p => p.potential));
  const anchorPotentials = anchors.map(rankMap.value);
  const cubic = createMonotoneCubicMap(indices, anchorPotentials, { derivatives });
  const indexMap = endSpacing === 'vinokur' ? withVinokurEnds(cubic) : cubic;
  const potentials = ranks.map((rank, i) => indices.includes(i) ? rankMap.value(rank)
    : (1 - spread) * rankMap.value(rank) + spread * indexMap.value(i));
  if (!increasing(potentials)) throw new Error('Unresolved outer streamtube station spacing.');
  return { potentials, rankMap, indexMap, diagnostics: {
    indexInterpolation: derivatives, endSpacing, spread, indices, anchorPotentials,
    endSegments: indexMap.endSegments ?? [],
    source: 'prescribed block targets; interpolated targets at shielded stations',
    exactModernMsetX: false,
  } };
}
