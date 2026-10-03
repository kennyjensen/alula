// SPDX-License-Identifier: GPL-2.0-or-later
// A nested station map: retain every old station and subdivide only the
// requested parent intervals. Coordinates refer to the old station index.
export function nestedIntervalMap(intervals, { factor, subdivisions } = {}) {
  if (!Number.isInteger(intervals) || intervals < 1
    || factor !== undefined && subdivisions !== undefined)
    throw new Error('Invalid nested refinement interval controls.');
  const uniform = factor === undefined ? 2 : factor;
  if (subdivisions === undefined && (!Number.isInteger(uniform) || uniform < 1 || uniform > 4))
    throw new Error('Invalid nested refinement factor.');
  const counts = subdivisions === undefined ? Array(intervals).fill(uniform) : subdivisions;
  if (!Array.isArray(counts) || counts.length !== intervals
    || counts.some(n => !Number.isInteger(n) || n < 1 || n > 4))
    throw new Error('Invalid nested refinement subdivisions.');
  const retained = [0], parents = [0], fractions = [0];
  for (let j = 0; j < intervals; j++) {
    for (let k = 1; k <= counts[j]; k++) {
      parents.push(k === counts[j] ? j + 1 : j);
      fractions.push(k === counts[j] ? 0 : k / counts[j]);
    }
    retained.push(retained.at(-1) + counts[j]);
  }
  const uniformFactor = counts.every(n => n === counts[0]) ? counts[0] : null;
  // Preserve the previous i/factor arithmetic for uniform refinements,
  // including branch-local BL coordinates with a nonzero LE/TE origin.
  const coordinate = (i, origin = 0) => uniformFactor === null
    ? parents[i] - origin + fractions[i] : (i - retained[origin]) / uniformFactor;
  return { counts: [...counts], retained, uniformFactor, coordinate,
    coordinates: parents.map((_, i) => coordinate(i)) };
}
