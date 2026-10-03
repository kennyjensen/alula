// SPDX-License-Identifier: GPL-2.0-or-later
// Physical lower-to-upper indexing of the cut intrinsic body grid. The dummy
// farfield strip in Drela's logical wrap is represented by two outer pressure
// boundaries. Outside a body, the two cut nodes share one position unknown;
// pressure continuity remains an explicit equation. No body cell is created.
export function createStreamtubeBodyLayout({ segments, tubes, bodies, primaryBody = 0, densityUnknowns = true, independentWakeBanks = false }) {
  const nx = segments, elements = bodies?.length;
  if (!Number.isInteger(nx) || nx < 4 || !Array.isArray(bodies) || !elements
    || !Array.isArray(tubes) || tubes.length !== elements + 1 || !tubes.every(n => Number.isInteger(n) && n > 0)
    || !Number.isInteger(primaryBody) || primaryBody < 0 || primaryBody >= elements || typeof densityUnknowns !== 'boolean' || typeof independentWakeBanks !== 'boolean')
    throw new Error('Invalid intrinsic body topology dimensions.');
  bodies = bodies.map(b => ({ ...b, roundedLeadingEdge: b.roundedLeadingEdge !== false })); tubes = tubes.slice();
  if (bodies.some(b => !Number.isInteger(b.leadingIndex) || !Number.isInteger(b.trailingIndex)
    || b.leadingIndex < 1 || b.trailingIndex >= nx || b.leadingIndex >= b.trailingIndex))
    throw new Error('Keep each body leading/trailing cross-line strictly inside the domain.');
  const active = (body, i) => i >= bodies[body].leadingIndex && i <= bodies[body].trailingIndex;
  const nt = tubes.reduce((s, n) => s + n, 0), groupOffsets = [0];
  for (const count of tubes) groupOffsets.push(groupOffsets.at(-1) + count);
  const densityIndex = (i, group, tube) => {
    if (!densityUnknowns) throw new Error('Incompressible density is eliminated, not an unknown.');
    return i * nt + groupOffsets[group] + tube;
  };
  const densityCount = densityUnknowns ? nx * nt : 0, positions = [], sharedCuts = new Map();
  let next = densityCount;
  const nodes = tubes.map((count, group) => Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: count + 1 }, (_, j) => {
    const body = j === 0 && group > 0 ? group - 1 : j === count && group < elements ? group : null;
    const side = j === 0 ? 'upper' : 'lower';
    if (body !== null && active(body, i)) return { kind: 'wall', body, side, i, column: null };
    if (body !== null) {
      const independent = independentWakeBanks && i > bodies[body].trailingIndex;
      const key = `${body}:${i}${independent ? `:${side}` : ''}`;
      if (!sharedCuts.has(key)) {
        const node = { kind: 'cut', body, i, ...(independent ? { side } : {}), column: next++ }; sharedCuts.set(key, node); positions.push(node);
      }
      return sharedCuts.get(key);
    }
    const node = { kind: (group === 0 && j === 0) || (group === elements && j === count) ? 'farfield' : 'interior',
      group, j, i, column: next++ };
    positions.push(node); return node;
  })));
  const positionCount = next - densityCount, globalOffset = next;
  const globals = { circulation: next++, stagnation: bodies.map(b => b.roundedLeadingEdge ? next++ : null),
    capture: bodies.map((_, b) => b !== primaryBody ? next++ : null), source: next++, doubletX: next++, doubletY: next++ };
  const n = next, rows = [], push = row => { row.index = rows.length; rows.push(row); };
  if (densityUnknowns) for (let i = 1; i < nx; i++) for (let group = 0; group < tubes.length; group++)
    for (let tube = 0; tube < tubes[group]; tube++) push({ kind: 'streamwise', i, group, tube });
  if (densityUnknowns) for (let group = 0; group < tubes.length; group++) for (let tube = 0; tube < tubes[group]; tube++)
    push({ kind: 'inletDensity', group, tube });
  for (let i = 1; i < nx; i++) {
    for (let group = 0; group < tubes.length; group++) for (let j = 1; j < tubes[group]; j++)
      push({ kind: 'internalPressure', i, group, j });
    push({ kind: 'farfieldPressure', i, side: 'lower' }); push({ kind: 'farfieldPressure', i, side: 'upper' });
    for (let body = 0; body < elements; body++) if (!active(body, i)) push({ kind: 'cutPressure', i, body });
  }
  // One centerline outlet tangency plus a gap equation, not two bank
  // tangencies plus a gap (which would overdetermine the wake endpoint).
  for (const i of [0, nx]) for (const node of positions.filter(node => node.i === i && !(independentWakeBanks && i === nx && node.kind === 'cut' && node.side === 'upper')))
    push({ kind: 'endTangency', i, node });
  if (independentWakeBanks) for (let body = 0; body < elements; body++) for (let i = bodies[body].trailingIndex + 1; i <= nx; i++)
    push({ kind: 'wakeGap', body, i });
  for (let body = 0; body < elements; body++) push({ kind: 'trailingKutta', i: bodies[body].trailingIndex, body });
  for (let body = 0; body < elements; body++) if (bodies[body].roundedLeadingEdge)
    push({ kind: 'leadingKutta', i: bodies[body].leadingIndex, body });
  for (const mode of ['source', 'doubletX', 'doubletY']) push({ kind: 'farfieldMatch', mode });
  if (rows.length !== n) throw new Error(`Intrinsic body count is not square: ${n} unknowns, ${rows.length} rows.`);
  const rowCounts = {};
  for (const row of rows) rowCounts[row.kind] = (rowCounts[row.kind] ?? 0) + 1;
  return { nx, nt, elements, tubes, bodies, primaryBody, n, densityUnknowns, independentWakeBanks, densityCount, positionCount, globalOffset,
    globalCount: n - globalOffset, globals, nodes, positions, rows, rowCounts, groupOffsets, densityIndex, active };
}

// Pair each physical equation with its unknown before AMD constructs A+A^T.
// The residual's family-wise row numbering is unrelated to column numbering;
// using it directly creates artificial long-range edges and excessive LU fill.
// This is an equation permutation only, not a change of variables or equations.
export function streamtubeEquationOrder(layout) {
  const order = new Int32Array(layout.n).fill(-1);
  for (const row of layout.rows) {
    let column;
    switch (row.kind) {
      case 'streamwise': column = layout.densityIndex(row.i, row.group, row.tube); break;
      case 'inletDensity': column = layout.densityIndex(0, row.group, row.tube); break;
      case 'internalPressure': column = layout.nodes[row.group][row.i][row.j].column; break;
      case 'farfieldPressure': {
        const group = row.side === 'lower' ? 0 : layout.elements;
        column = layout.nodes[group][row.i][row.side === 'lower' ? 0 : layout.tubes[group]].column;
        break;
      }
      case 'cutPressure': column = layout.nodes[row.body][row.i][layout.tubes[row.body]].column; break;
      case 'endTangency': column = row.node.column; break;
      case 'wakeGap': column = layout.nodes[row.body + 1][row.i][0].column; break;
      case 'leadingKutta': column = layout.globals.stagnation[row.body]; break;
      case 'trailingKutta': column = row.body === layout.primaryBody ? layout.globals.circulation : layout.globals.capture[row.body]; break;
      case 'farfieldMatch': column = layout.globals[row.mode]; break;
    }
    if (!Number.isInteger(column) || column < 0 || column >= layout.n || order[column] !== -1)
      throw new Error(`Cannot pair streamtube equation ${row.kind}.`);
    order[column] = row.index;
  }
  if (order.some(row => row < 0)) throw new Error('Incomplete streamtube equation ordering.');
  return order;
}

// Cumulative captured mass levels delimit E+1 fluid groups. Total mass and
// the primary dividing level are prescribed; E-1 other levels are unknown.
// Fixed positive weights distribute each group's mass among its streamtubes.
export function allocateStreamtubeMasses(levels, weights, primaryBody = 0) {
  if (!Array.isArray(levels) || !levels.every(Number.isFinite) || levels.length < 3
    || levels.some((v, i) => i && v <= levels[i - 1]) || !Array.isArray(weights) || weights.length !== levels.length - 1
    || !Number.isInteger(primaryBody) || primaryBody < 0 || primaryBody >= levels.length - 2)
    throw new Error('Invalid captured streamtube mass levels.');
  const groups = weights.map((row, group) => {
    if (!Array.isArray(row) || !row.length || !row.every(w => Number.isFinite(w) && w > 0)) throw new Error('Invalid streamtube mass weights.');
    const sum = row.reduce((s, w) => s + w, 0), mass = levels[group + 1] - levels[group];
    if (!Number.isFinite(sum) || !Number.isFinite(mass)) throw new Error('Nonfinite streamtube mass allocation.');
    return row.map(w => {
      const fraction = w / sum, derivatives = new Map();
      if (group > 0 && group - 1 !== primaryBody) derivatives.set(group - 1, -fraction);
      if (group < weights.length - 1 && group !== primaryBody) derivatives.set(group, fraction);
      return { massFlow: fraction * mass, derivatives };
    });
  });
  return { groups, totalMass: levels.at(-1) - levels[0], primaryLevel: levels[primaryBody + 1] };
}
