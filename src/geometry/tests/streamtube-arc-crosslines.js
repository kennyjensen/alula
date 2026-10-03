// SPDX-License-Identifier: GPL-2.0-or-later
// Potential defines the LE/TE block anchors. Surface spacing inside each
// block is reconciled in physical arc, so stagnation's vanishing velocity
// cannot turn an O(ds) surface demand into an O(ds^2) demand on another foil.
// Ordinary positive C1 boundary interpolation removes artificial derivative
// jumps at foreign LE/TE anchors. This is not a recovered MSET spacing law.
import { createMonotoneCubicMap } from '../../numerics/monotone-cubic.js';

export function createArcCrosslineMap({ anchors, arcAt }) {
  if (!Array.isArray(anchors) || anchors.length < 2 || typeof arcAt !== 'function'
    || anchors.some((v, i) => !Number.isFinite(v) || i && v <= anchors[i - 1]))
    throw new Error('Invalid surface arc block anchors.');
  const coordinate = Object.freeze(anchors.slice()), arc = Object.freeze(coordinate.map(arcAt));
  if (arc.some((v, i) => !Number.isFinite(v) || i && v <= arc[i - 1]))
    throw new Error('Surface arc must increase through its block anchors.');
  const map = createMonotoneCubicMap(coordinate, arc);
  const inverse = value => {
    if (!Number.isFinite(value) || value < arc[0] || value > arc.at(-1)) throw new Error('Surface arc cross-line lies outside its blocks.');
    let lo = 0, hi = arc.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (arc[mid] <= value) lo = mid; else hi = mid; }
    if (value === arc[lo]) return coordinate[lo];
    if (value === arc[hi]) return coordinate[hi];
    // Invert the actual forward cubic. Interpolating the swapped data would
    // define a second cubic and would not be its mathematical inverse.
    let left = coordinate[lo], right = coordinate[hi];
    for (let iteration = 0; iteration < 64; iteration++) {
      const mid = left + .5 * (right - left);
      if (mid === left || mid === right) break;
      const actual = map.value(mid);
      if (actual === value) return mid;
      if (actual < value) left = mid; else right = mid;
    }
    return Math.abs(map.value(left) - value) <= Math.abs(map.value(right) - value) ? left : right;
  };
  return { coordinate, arc, value: map.value, evaluate: map.evaluate, inverse,
    interpolation: 'positive C1 cubic through fixed block anchors' };
}
