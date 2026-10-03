// SPDX-License-Identifier: GPL-2.0-or-later
// Integer prefix differences: I[to]-I[from] >= minimum. Fixed individual
// block counts add reverse edges. Longest paths give the least feasible
// prefix indices, hence minimum total count for the supplied lower bounds.
export function allocateOrderedIntervalCounts({ minimumCounts, requirements = [], fixedCounts = [], maximumTotal = 2000 }) {
  const n = minimumCounts?.length;
  if (!Array.isArray(minimumCounts) || !n || Array.from(minimumCounts).some(v => !Number.isSafeInteger(v) || v < 1)
    || !Array.isArray(requirements) || requirements.some(r => !r || !Number.isInteger(r.from) || !Number.isInteger(r.to)
      || r.from < 0 || r.from >= r.to || r.to > n || !Number.isSafeInteger(r.minimum) || r.minimum < 1)
    || !Array.isArray(fixedCounts) || fixedCounts.some(r => !r || !Number.isInteger(r.block) || r.block < 0 || r.block >= n
      || !Number.isSafeInteger(r.intervals) || r.intervals < 1)
    || new Set(fixedCounts.map(r => r.block)).size !== fixedCounts.length
    || !Number.isSafeInteger(maximumTotal) || maximumTotal < 1) throw new Error('Invalid ordered interval constraints.');
  const edges = minimumCounts.map((minimum, from) => ({ from, to: from + 1, minimum }));
  for (const { block, intervals } of fixedCounts) edges.push({ from: block, to: block + 1, minimum: intervals }, { from: block + 1, to: block, minimum: -intervals });
  edges.push(...requirements);
  const indices = [0, ...Array(n).fill(-Infinity)];
  let changed = false;
  for (let pass = 0; pass <= n; pass++) {
    changed = false;
    for (const e of edges) {
      const target = indices[e.from] + e.minimum;
      if (target > indices[e.to]) { indices[e.to] = target; changed = true; }
    }
    if (!changed) break;
  }
  if (changed || indices[0] !== 0) throw new Error('Fixed interval counts conflict with the required passage counts.');
  if (indices[n] > maximumTotal) throw new Error(`Required interval total ${indices[n]} exceeds budget ${maximumTotal}.`);
  const counts = indices.slice(1).map((v, i) => v - indices[i]);
  if (!indices.every(Number.isSafeInteger) || edges.some(e => indices[e.to] - indices[e.from] < e.minimum))
    throw new Error('Unresolved integer interval allocation.');
  return { counts, indices, total: indices[n], minimumTotalForBounds: true,
    tieBreak: 'least feasible prefix indices; additional underdetermined intervals go downstream' };
}
