// SPDX-License-Identifier: GPL-2.0-or-later
// Positive geometric intervals with a prescribed near-body spacing and an
// exact far boundary. Fit the common ratio instead of clipping a final step,
// which would introduce a small last cell and an abrupt spacing reversal.
export function gradedIntervals({ extent, firstSpacing, intervals }) {
  if (![extent, firstSpacing].every(Number.isFinite) || extent <= 0 || firstSpacing <= 0
    || !Number.isInteger(intervals) || intervals < 2 || intervals > 512)
    throw new Error('Invalid graded interval controls.');
  // If the requested first interval exceeds the average, a uniform grid is
  // already finer everywhere. Never shrink the intervals toward infinity.
  const first = Math.min(firstSpacing, extent / intervals);
  const sum = logarithm => {
    let value = 0;
    for (let k = 0; k < intervals; k++) value += first * Math.exp(logarithm * k);
    return value;
  };
  let lower = 0, upper = Math.log(extent / first) / (intervals - 1);
  for (let k = 0; k < 64; k++) {
    const mid = .5 * (lower + upper);
    if (sum(mid) < extent) lower = mid; else upper = mid;
  }
  const logarithm = .5 * (lower + upper), x = [0];
  for (let k = 0; k < intervals; k++) x.push(k === intervals - 1 ? extent : x.at(-1) + first * Math.exp(logarithm * k));
  if (x.some((v, i) => !Number.isFinite(v) || (i && v <= x[i - 1])))
    throw new Error('Unresolved graded intervals.');
  return { x, firstSpacing: first, growth: Math.exp(logarithm), intervals };
}
