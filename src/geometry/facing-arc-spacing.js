// SPDX-License-Identifier: GPL-2.0-or-later
// Reconcile facing wall parameterizations before equidistributing demands.
// At each shared physical hit, fit the ratio of positive arc derivatives to
// the incident blocks' arc-length ratios in least squared logarithms. Then
// choose the nearest pair of original derivatives in logarithms, subject to
// the sufficient positive-Hermite cone m <= 2*min(adjacent secants).
// This is a stated reconstruction of MSET's facing-spacing stage, not its
// unavailable source formula. Shared cut schedules are not changed here.
export function reconcileFacingArcSlopes(first, second) {
  const prepare = map => {
    const { knots, values, slopes } = map;
    if (![knots, values, slopes].every(Array.isArray) || knots.length < 2
      || values.length !== knots.length || slopes.length !== knots.length
      || !Array.from(knots).every((v, i) => Number.isFinite(v) && (!i || v > knots[i - 1]))
      || !Array.from(values).every((v, i) => Number.isFinite(v) && (!i || v > values[i - 1]))
      || !Array.from(slopes).every(v => Number.isFinite(v) && v > 0))
      throw new Error('Invalid facing arc map.');
    const secants = knots.slice(1).map((v, i) => (values[i + 1] - values[i]) / (v - knots[i]));
    if (!secants.every(v => v > 0 && Number.isFinite(v))) throw new Error('Unresolved facing arc secant.');
    return { knots: knots.slice(), slopes: slopes.slice(), secants,
      upper: knots.map((_, i) => 2 * Math.min(i ? secants[i - 1] : Infinity, i < secants.length ? secants[i] : Infinity)) };
  };
  const a = prepare(first), b = prepare(second), start = Math.max(a.knots[0], b.knots[0]), end = Math.min(a.knots.at(-1), b.knots.at(-1));
  const result = { first: a.slopes, second: b.slopes, knots: [], start, end,
    method: 'least logarithmic change of positive facing arc slopes', exactModernMsetLaw: false };
  if (!(start < end)) return result;
  const anchors = [...new Set([...a.knots, ...b.knots].filter(v => v >= start && v <= end))].sort((u, v) => u - v);
  if (anchors.some(v => !a.knots.includes(v) || !b.knots.includes(v)))
    throw new Error('Facing arc blocks require shared physical hit ranks.');
  for (const [k, rank] of anchors.entries()) {
    const ia = a.knots.indexOf(rank), ib = b.knots.indexOf(rank), ratios = [];
    if (k) ratios.push(Math.log(a.secants[ia - 1]) - Math.log(b.secants[ib - 1]));
    if (k < anchors.length - 1) ratios.push(Math.log(a.secants[ia]) - Math.log(b.secants[ib]));
    const logRatio = ratios.reduce((sum, v) => sum + v, 0) / ratios.length;
    const logProduct = Math.log(first.slopes[ia]) + Math.log(second.slopes[ib]);
    const logA = .5 * (logProduct + logRatio), logB = .5 * (logProduct - logRatio);
    const shift = Math.min(0, Math.log(a.upper[ia]) - logA, Math.log(b.upper[ib]) - logB);
    a.slopes[ia] = Math.min(a.upper[ia], Math.exp(logA + shift));
    b.slopes[ib] = Math.min(b.upper[ib], Math.exp(logB + shift));
    if (!(a.slopes[ia] > 0 && b.slopes[ib] > 0)) throw new Error('Unresolved positive facing arc slopes.');
    result.knots.push({ rank, firstIndex: ia, secondIndex: ib,
      before: [first.slopes[ia], second.slopes[ib]], after: [a.slopes[ia], b.slopes[ib]],
      incidentLogRatios: ratios, targetLogRatio: logRatio, commonLogAdjustment: shift });
  }
  return result;
}
