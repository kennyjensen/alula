// SPDX-License-Identifier: GPL-2.0-or-later
// Optional station ordering from the ACTUAL numerical sparsity pattern. No
// governing system, matrix coefficients, scaling or numerical pivot changes.
import { createStreamtubeBodyLayout } from './streamtube-body-layout.js';

const fail = message => { throw new Error(`Station ordering: ${message}`); };
const positiveInteger = n => Number.isInteger(n) && n > 0 && n < 0x80000000;

function enumerateStations(layout) {
  const stations = [];
  for (let body = 0; body < layout.elements; body++) {
    const b = layout.bodies[body];
    for (const side of ['upper', 'lower']) for (let i = b.leadingIndex + 1; i <= b.trailingIndex; i++)
      stations.push({ kind: 'surface', body, side, i });
    for (let i = b.trailingIndex; i <= layout.nx; i++) stations.push({ kind: 'wake', body, i });
  }
  return stations;
}

function liveStationMetadata(layout, stations) {
  if (!layout || !positiveInteger(layout.nx) || !positiveInteger(layout.n) || layout.densityUnknowns !== true
    || !Array.isArray(layout.bodies) || layout.elements !== layout.bodies.length || !layout.elements
    || !Array.isArray(layout.tubes) || layout.tubes.length !== layout.elements + 1
    || !Array.from(layout.tubes).every(positiveInteger) || !Array.isArray(layout.rows) || layout.rows.length !== layout.n
    || !Array.isArray(layout.positions) || typeof layout.densityIndex !== 'function'
    || !Array.isArray(stations)) fail('complete live compressible layout and BL station list required.');
  if (Array.from(layout.rows).some((row, index) => !row || row.index !== index)
    || Array.from(layout.bodies).some(b => !b || !positiveInteger(b.leadingIndex)
      || !positiveInteger(b.trailingIndex) || b.leadingIndex >= b.trailingIndex || b.trailingIndex >= layout.nx))
    fail('invalid live layout row or body indices.');
  const canonical = enumerateStations(layout);
  if (stations.length !== canonical.length || canonical.some((expected, id) => {
    const actual = stations[id];
    return !actual || actual.id !== id || actual.kind !== expected.kind || actual.body !== expected.body
      || actual.side !== expected.side || actual.i !== expected.i
      || (actual.k !== undefined && actual.k !== actual.i - (expected.kind === 'surface'
        ? layout.bodies[expected.body].leadingIndex : layout.bodies[expected.body].trailingIndex));
  })) fail('live BL stations must retain complete ordered ids and four-slot surface/wake topology.');
  // Regime, phase, arc and physical state are deliberately not metadata for
  // this ordering. Keep the same canonical records as checkpoint enumeration.
  return buildMetadata(layout, canonical);
}

// Pure index enumeration, mirroring the four BL entries per station in
// streamtube-boundary-layers.js and the edge slot in streamtube-coupled.js.
export function streamtubeStationMetadata(checkpoint) {
  const r = checkpoint?.restart, input = r?.input, nodes = r?.initialEuler?.nodes;
  if (checkpoint?.version !== 1 || !input || !Array.isArray(nodes) || nodes.length < 2
    || !Array.isArray(nodes[0]) || !Array.isArray(nodes[0][0])) fail('complete checkpoint topology required.');
  const nx = nodes[0].length - 1, tubes = nodes.map(g => g?.[0]?.length - 1);
  if (!positiveInteger(nx) || !tubes.every(positiveInteger)
    || nodes.some((g, k) => !Array.isArray(g) || g.length !== nx + 1
      || Array.from(g).some(row => !Array.isArray(row) || row.length !== tubes[k] + 1))) fail('inconsistent node topology.');
  const layout = createStreamtubeBodyLayout({ segments: nx, tubes, bodies: input.bodies,
    primaryBody: input.primaryBody, densityUnknowns: true,
    independentWakeBanks: input.wakeGeometry === 'independent-banks' });
  if (r.initialEuler.x?.length !== layout.n) fail('Euler packed length does not match layout.');
  const stations = enumerateStations(layout);
  if (r.initialBL?.length !== 4 * stations.length) fail('BL packed length does not match station enumeration.');
  return buildMetadata(layout, stations);
}

function buildMetadata(layout, stations) {
  const { nx, tubes } = layout;
  const n = layout.n + 4 * stations.length, rows = [], columns = Array(n), preferred = new Int32Array(n).fill(-1);
  for (let i = 0; i < nx; i++) for (let g = 0; g < tubes.length; g++) for (let j = 0; j < tubes[g]; j++)
    columns[layout.densityIndex(i, g, j)] = { kind: 'density', i, group: g, tube: j };
  for (const p of layout.positions) columns[p.column] = { ...p };
  for (const [name, values] of Object.entries(layout.globals))
    (Array.isArray(values) ? values : [values]).forEach((column, body) => {
      if (column !== null) columns[column] = { kind: 'global', name, ...(Array.isArray(values) ? { body } : {}), i: nx + 1 };
    });
  const nodeColumn = (g, i, j) => layout.nodes[g][i][j].column;
  for (const row of layout.rows) {
    const i = row.i ?? (row.kind === 'inletDensity' ? 0 : nx + 1);
    rows.push({ ...row, i });
    let column = null;
    if (row.kind === 'streamwise') column = layout.densityIndex(i, row.group, row.tube);
    else if (row.kind === 'inletDensity') column = layout.densityIndex(0, row.group, row.tube);
    else if (row.kind === 'internalPressure') column = nodeColumn(row.group, i, row.j);
    else if (row.kind === 'farfieldPressure') column = row.side === 'lower' ? nodeColumn(0, i, 0) : nodeColumn(tubes.length - 1, i, tubes.at(-1));
    else if (row.kind === 'cutPressure') column = nodeColumn(row.body, i, tubes[row.body]);
    else if (row.kind === 'endTangency') column = row.node.column;
    else if (row.kind === 'wakeGap') column = nodeColumn(row.body + 1, i, 0);
    else if (row.kind === 'leadingKutta') column = layout.globals.stagnation[row.body];
    else if (row.kind === 'farfieldMatch') column = layout.globals[row.mode];
    // In particular, trailing Kutta has no assumed direct circulation entry.
    if (column !== null && column !== undefined) preferred[row.index] = column;
  }
  stations.forEach((station, id) => {
    for (let k = 0; k < 4; k++) {
      const index = layout.n + 4 * id + k;
      columns[index] = { ...station, station: id, kind: `bl-${['aux', 'theta', 'delta', 'ue'][k]}` };
      rows.push({ ...station, index, station: id, kind: `bl-${['aux', 'momentum', 'energy', 'edge'][k]}` });
      preferred[index] = index;
    }
  });
  if (rows.length !== n || Array.from(columns).some(q => !q)) fail('incomplete row/column enumeration.');
  return { n, rows, columns, preferred, layout: { nx, tubes, bodies: layout.bodies.map(({ leadingIndex, trailingIndex }) => ({ leadingIndex, trailingIndex })),
    eulerUnknowns: layout.n, densityUnknowns: layout.densityCount, positionUnknowns: layout.positionCount,
    globalUnknowns: layout.globalCount, boundaryLayerStations: stations.length } };
}

export function createStreamtubeStationOrdering({ matrix, checkpoint, layout, stations } = {}) {
  const fromCheckpoint = checkpoint !== undefined, fromLive = layout !== undefined || stations !== undefined;
  if (fromCheckpoint === fromLive || fromLive && (layout === undefined || stations === undefined))
    fail('supply either checkpoint or both live layout and stations, exclusively.');
  const metadata = fromCheckpoint ? streamtubeStationMetadata(checkpoint) : liveStationMetadata(layout, stations);
  const { n, rows, columns, preferred } = metadata;
  const { rowPtr, colIndex, values } = matrix ?? {};
  if (matrix?.n !== n || rowPtr?.length !== n + 1 || rowPtr[0] !== 0
    || rowPtr[n] !== values?.length || colIndex?.length !== values?.length) fail('invalid square CSR dimensions.');
  const neighbors = Array.from({ length: n }, () => []), columnDegrees = new Int32Array(n);
  const rowMax = new Float64Array(n);
  let exactZeros = 0, numericalNonzeros = 0;
  for (let r = 0; r < n; r++) {
    if (!Number.isInteger(rowPtr[r]) || !Number.isInteger(rowPtr[r + 1]) || rowPtr[r + 1] < rowPtr[r]
      || rowPtr[r + 1] > values.length) fail('invalid CSR row pointer.');
    const seen = new Set();
    for (let p = rowPtr[r]; p < rowPtr[r + 1]; p++) {
      const c = colIndex[p], v = values[p];
      if (!Number.isInteger(c) || c < 0 || c >= n || !Number.isFinite(v) || seen.has(c))
        fail('CSR columns must be unique within each row, in range, and finite.');
      seen.add(c); if (v === 0) { exactZeros++; continue; }
      neighbors[r].push(c); columnDegrees[c]++; numericalNonzeros++; rowMax[r] = Math.max(rowMax[r], Math.abs(v));
    }
  }
  // Preference is topological, not a magnitude-based pivot or equilibration:
  // own physical variable if present, nonglobal/local station first, then
  // lower column degree and original column index for deterministic ties.
  for (let r = 0; r < n; r++) neighbors[r].sort((a, b) =>
    Number(b === preferred[r]) - Number(a === preferred[r])
    || Number(columns[a].kind === 'global') - Number(columns[b].kind === 'global')
    || Math.abs(rows[r].i - columns[a].i) - Math.abs(rows[r].i - columns[b].i)
    || columnDegrees[a] - columnDegrees[b] || a - b);
  const rowToColumn = new Int32Array(n).fill(-1), columnToRow = new Int32Array(n).fill(-1);
  let naturalMatches = 0, missingNaturalEntries = 0, greedyMatches = 0, augmentations = 0, maximumAugmentingRows = 0;
  const assign = (r, c) => { rowToColumn[r] = c; columnToRow[c] = r; };
  for (let r = 0; r < n; r++) {
    const c = preferred[r];
    if (c < 0) continue;
    if (!neighbors[r].includes(c)) { missingNaturalEntries++; continue; }
    if (columnToRow[c] === -1) { assign(r, c); naturalMatches++; }
  }
  for (let r = 0; r < n; r++) if (rowToColumn[r] === -1) {
    const c = neighbors[r].find(c => columnToRow[c] === -1);
    if (c !== undefined) { assign(r, c); greedyMatches++; }
  }
  // Breadth-first alternating paths, without recursion or a finite search
  // depth. A failed search is an actual unmatched component of this pattern;
  // no zero perturbation or discarded equation repairs it.
  const queue = new Int32Array(n), seenRows = new Int32Array(n), parentRow = new Int32Array(n), viaColumn = new Int32Array(n);
  let stamp = 0;
  function augment(root) {
    let head = 0, tail = 0; stamp++; queue[tail++] = root; seenRows[root] = stamp;
    while (head < tail) {
      const r = queue[head++];
      for (const c of neighbors[r]) {
        const owner = columnToRow[c];
        if (owner === -1) {
          let row = r, column = c, length = 0;
          while (true) {
            assign(row, column); length++;
            if (row === root) break;
            column = viaColumn[row]; row = parentRow[row];
          }
          augmentations++; maximumAugmentingRows = Math.max(maximumAugmentingRows, length); return true;
        }
        if (seenRows[owner] !== stamp) {
          seenRows[owner] = stamp; parentRow[owner] = r; viaColumn[owner] = c; queue[tail++] = owner;
        }
      }
    }
    return false;
  }
  for (let r = 0; r < n; r++) if (rowToColumn[r] === -1) augment(r);
  const unmatchedRows = [], unmatchedColumns = [];
  for (let k = 0; k < n; k++) {
    if (rowToColumn[k] === -1) unmatchedRows.push({ index: k, ...rows[k] });
    if (columnToRow[k] === -1) unmatchedColumns.push({ index: k, ...columns[k] });
  }
  if (unmatchedRows.length) {
    const error = new Error('Station ordering: actual nonzero pattern has no perfect row/column matching.');
    error.code = 'STREAMTUBE_STATION_MATCHING_FAILED'; error.diagnostics = { unmatchedRows, unmatchedColumns, naturalMatches, greedyMatches, augmentations };
    throw error;
  }
  const pairs = Array.from(rowToColumn, (c, r) => ({ row: r, column: c,
    station: Math.max(rows[r].i, columns[c].i),
    border: columns[c].kind === 'global' || rows[r].kind === 'farfieldMatch' }));
  pairs.sort((a, b) => Number(a.border) - Number(b.border) || a.station - b.station || a.row - b.row);
  const Puser = Int32Array.from(pairs, q => q.row), Quser = Int32Array.from(pairs, q => q.column);
  const borderPairs = [], displacedNaturalPairs = [];
  let minimumDiagonalMagnitude = Infinity, minimumDiagonalToRowMaximum = Infinity, maximumLocalStationMismatch = 0;
  for (const q of pairs) {
    const { row: r, column: c } = q;
    if (columnToRow[c] !== r || !neighbors[r].includes(c)) fail('internal matching verification failed.');
    let value;
    for (let p = rowPtr[r]; p < rowPtr[r + 1]; p++) if (colIndex[p] === c) { value = values[p]; break; }
    if (!Number.isFinite(value) || value === 0) fail('matched diagonal is not an actual finite nonzero.');
    minimumDiagonalMagnitude = Math.min(minimumDiagonalMagnitude, Math.abs(value));
    minimumDiagonalToRowMaximum = Math.min(minimumDiagonalToRowMaximum, Math.abs(value) / rowMax[r]);
    const record = { row: r, column: c, rowMetadata: rows[r], columnMetadata: columns[c], value,
      stationMismatch: Math.abs(rows[r].i - columns[c].i) };
    if (q.border) borderPairs.push(record); else maximumLocalStationMismatch = Math.max(maximumLocalStationMismatch, record.stationMismatch);
    if (preferred[r] !== c) displacedNaturalPairs.push(record);
  }
  return { Puser, Quser, layout: metadata.layout, diagnostics: { n, numericalNonzeros, exactZeros,
    naturalMatches, missingNaturalEntries, greedyMatches, augmentations, maximumAugmentingRows,
    allSelectedDiagonalsFiniteNonzero: true, bijective: true, minimumDiagonalMagnitude, minimumDiagonalToRowMaximum,
    maximumLocalStationMismatch, borderPairCount: borderPairs.length, borderPairs, displacedNaturalPairs,
    convention: 'Puser[k],Quser[k] are original row/column indices at new position k; Aordered[k,l]=A[Puser[k],Quser[l]].',
    stationOrder: 'Matched pair max(row station,column station), then original row; border pairs last.',
    border: 'Union of matched pairs containing a global column or dense farfieldMatch row. Kutta rows have no prescribed global diagonal.',
    cumulativeWakeBL: 'All actual nonlocal entries retained; only matched pairs are ordered.',
    numericalPivotingPerformed: false, factorizationPerformed: false, speedupClaim: false } };
}
