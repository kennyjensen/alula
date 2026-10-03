// SPDX-License-Identifier: GPL-2.0-or-later
// A stated reconstruction of facing-spacing reconciliation, not MSET source.
// Fit incident block secant ratios in logarithms; merge shared cut slopes;
// solve the resulting forest exactly and choose the least logarithmic change
// of all original occurrences under the positive Hermite cone. Shared copies
// must use the same physical metric units (additive origins may differ).
export function reconcileFacingOutlineSlopes({ maps, pairs, shared = [] } = {}) {
  if (!Array.isArray(maps) || !maps.length || !Array.isArray(pairs) || !Array.isArray(shared))
    throw new Error('Invalid facing outline maps, pairs or shared cuts.');
  const prepared = Array.from(maps, map => {
    const { knots, values, slopes } = map ?? {};
    if (![knots, values, slopes].every(Array.isArray) || knots.length < 2
      || values.length !== knots.length || slopes.length !== knots.length
      || Array.from(knots).some((v, i) => !Number.isFinite(v) || i && !(v > knots[i - 1]))
      || Array.from(values).some((v, i) => !Number.isFinite(v) || i && !(v > values[i - 1]))
      || Array.from(slopes).some(v => !Number.isFinite(v) || !(v > 0)))
      throw new Error('Invalid facing outline map.');
    const secants = knots.slice(1).map((v, i) => (values[i + 1] - values[i]) / (v - knots[i]));
    const upper = knots.map((_, i) => 2 * Math.min(i ? secants[i - 1] : Infinity,
      i < secants.length ? secants[i] : Infinity));
    if (![...secants, ...upper].every(v => Number.isFinite(v) && v > 0))
      throw new Error('Unresolved finite positive facing outline cone.');
    return { knots: knots.slice(), slopes: slopes.slice(), secants, upper };
  });
  const validMap = index => Number.isInteger(index) && index >= 0 && index < prepared.length;
  const ranks = new Map(), pairKeys = new Set();
  const atRank = rank => {
    if (!ranks.has(rank)) ranks.set(rank, { rank, pairs: [], shared: [] });
    return ranks.get(rank);
  };
  for (const pair of pairs) {
    if (!Array.isArray(pair) || pair.length !== 2 || !Array.from(pair).every(validMap) || pair[0] === pair[1])
      throw new Error('Invalid facing outline pair.');
    const [first, second] = pair, key = [...pair].sort((a, b) => a - b).join(':');
    if (pairKeys.has(key)) throw new Error('Duplicate facing outline pair.');
    pairKeys.add(key);
    const a = prepared[first], b = prepared[second];
    const start = Math.max(a.knots[0], b.knots[0]), end = Math.min(a.knots.at(-1), b.knots.at(-1));
    if (!(start < end)) continue;
    const anchors = [...new Set([...a.knots, ...b.knots].filter(v => v >= start && v <= end))].sort((u, v) => u - v);
    if (anchors.some(v => !a.knots.includes(v) || !b.knots.includes(v)))
      throw new Error('Facing outline blocks require shared physical hit ranks.');
    anchors.forEach((rank, k) => {
      const ia = a.knots.indexOf(rank), ib = b.knots.indexOf(rank), ratios = [], logs = [];
      for (const offset of [...(k ? [-1] : []), ...(k < anchors.length - 1 ? [0] : [])]) {
        const la = Math.log(a.secants[ia + offset]), lb = Math.log(b.secants[ib + offset]);
        ratios.push(la - lb); logs.push(la, lb);
      }
      atRank(rank).pairs.push({ maps: [...pair], indices: [ia, ib], incidentLogRatios: ratios,
        targetLogRatio: ratios.reduce((sum, v) => sum + v, 0) / ratios.length,
        roundoffTolerance: 64 * Number.EPSILON * Math.max(1, ...logs.map(Math.abs)) });
    });
  }
  for (const item of shared) {
    if (!item || !Number.isFinite(item.rank) || !Array.isArray(item.maps) || item.maps.length < 2
      || !Array.from(item.maps).every(validMap) || new Set(item.maps).size !== item.maps.length
      || item.maps.some(m => !prepared[m].knots.includes(item.rank)))
      throw new Error('Shared cuts require a present rank and distinct valid maps.');
    atRank(item.rank).shared.push([...item.maps]);
  }
  const slopes = prepared.map(m => m.slopes.slice()), knots = [];
  for (const record of [...ranks.values()].sort((a, b) => a.rank - b.rank)) {
    const active = [...new Set([...record.pairs.flatMap(p => p.maps), ...record.shared.flat()])].sort((a, b) => a - b);
    const parents = new Map(active.map(m => [m, m]));
    const root = m => { while (parents.get(m) !== m) m = parents.get(m); return m; };
    for (const group of record.shared) for (const m of group.slice(1)) {
      const a = root(group[0]), b = root(m); parents.set(Math.max(a, b), Math.min(a, b));
    }
    const variables = new Map();
    for (const m of active) {
      const id = root(m), index = prepared[m].knots.indexOf(record.rank);
      if (!variables.has(id)) variables.set(id, { id, occurrences: [], edges: [] });
      variables.get(id).occurrences.push({ map: m, index, before: prepared[m].slopes[index], upper: prepared[m].upper[index] });
    }
    record.pairs.forEach((pair, edge) => {
      const a = root(pair.maps[0]), b = root(pair.maps[1]);
      if (a === b) {
        if (Math.abs(pair.targetLogRatio) > pair.roundoffTolerance)
          throw new Error(`Inconsistent shared-cut self-edge demand at rank ${record.rank}.`);
        pair.redundantSharedConstraint = true;
      } else {
        variables.get(a).edges.push({ to: b, delta: -pair.targetLogRatio, edge });
        variables.get(b).edges.push({ to: a, delta: pair.targetLogRatio, edge });
      }
    });
    const visited = new Set(), components = [];
    for (const first of variables.keys()) {
      if (visited.has(first)) continue;
      const relative = new Map([[first, 0]]), stack = [{ id: first, incoming: -1 }], component = [];
      visited.add(first);
      while (stack.length) {
        const { id, incoming } = stack.pop(), variable = variables.get(id);
        component.push(variable);
        for (const edge of variable.edges) {
          if (edge.edge === incoming) continue;
          if (visited.has(edge.to)) throw new Error(`Facing outline cycles are unsupported at rank ${record.rank}.`);
          const offset = relative.get(id) + edge.delta;
          if (!Number.isFinite(offset)) throw new Error('Unresolved facing outline log ratios.');
          relative.set(edge.to, offset); visited.add(edge.to); stack.push({ id: edge.to, incoming: edge.edge });
        }
      }
      let count = 0, sum = 0, cap = Infinity;
      for (const variable of component) for (const occurrence of variable.occurrences) {
        sum += Math.log(occurrence.before) - relative.get(variable.id); count++;
        cap = Math.min(cap, Math.log(occurrence.upper) - relative.get(variable.id));
      }
      const unconstrainedLogShift = sum / count, shift = Math.min(unconstrainedLogShift, cap);
      if (![unconstrainedLogShift, shift].every(Number.isFinite)) throw new Error('Unresolved facing outline common shift.');
      let squaredLogChange = 0;
      for (const variable of component) {
        const bound = Math.min(...variable.occurrences.map(o => o.upper));
        const slope = Math.min(bound, Math.exp(relative.get(variable.id) + shift));
        if (!(slope > 0) || !Number.isFinite(slope)) throw new Error('Unresolved positive facing outline slope.');
        for (const occurrence of variable.occurrences) {
          slopes[occurrence.map][occurrence.index] = slope; occurrence.after = slope;
          squaredLogChange += (Math.log(slope) - Math.log(occurrence.before)) ** 2;
        }
      }
      components.push({ maps: component.flatMap(v => v.occurrences.map(o => o.map)).sort((a, b) => a - b),
        variables: component.map(v => v.occurrences.map(o => o.map)), occurrences: count,
        unconstrainedLogShift, commonLogShift: shift, commonLogAdjustment: shift - unconstrainedLogShift, squaredLogChange });
    }
    for (const pair of record.pairs) {
      pair.achievedLogRatio = Math.log(slopes[pair.maps[0]][pair.indices[0]]) - Math.log(slopes[pair.maps[1]][pair.indices[1]]);
      pair.incidentLogErrors = pair.incidentLogRatios.map(r => pair.achievedLogRatio - r);
      pair.rmsIncidentLogMismatch = Math.sqrt(pair.incidentLogErrors.reduce((s, v) => s + v * v, 0) / pair.incidentLogErrors.length);
    }
    knots.push({ ...record, occurrences: [...variables.values()].flatMap(v => v.occurrences), components });
  }
  for (const [m, row] of slopes.entries()) if (row.some((v, i) => !Number.isFinite(v) || !(v > 0) || v > prepared[m].upper[i]))
    throw new Error('Facing outline result lies outside the finite positive Hermite cone.');
  return { slopes, knots, method: 'least logarithmic change of shared facing outline slopes on a forest', exactModernMsetLaw: false };
}
