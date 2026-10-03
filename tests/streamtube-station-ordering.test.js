import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodyLayout } from '../src/euler/streamtube-body-layout.js';
import { createStreamtubeStationOrdering, streamtubeStationMetadata } from '../scripts/validation/streamtube-station-ordering.js';

function fixture({ multi = false, independent = true } = {}) {
  const nx = multi ? 9 : 6, tubes = multi ? [2, 1, 2] : [2, 2];
  const bodies = multi ? [{ leadingIndex: 2, trailingIndex: 5 }, { leadingIndex: 3, trailingIndex: 6, roundedLeadingEdge: false }]
    : [{ leadingIndex: 2, trailingIndex: 4 }];
  const layout = createStreamtubeBodyLayout({ segments: nx, tubes, bodies, independentWakeBanks: independent });
  const stations = bodies.reduce((sum, b) => sum + 2 * (b.trailingIndex - b.leadingIndex) + nx - b.trailingIndex + 1, 0);
  const checkpoint = { version: 1, restart: { input: { bodies, primaryBody: 0, flowModel: 'compressible',
    wakeGeometry: independent ? 'independent-banks' : 'centerline' }, initialEuler: { x: Array(layout.n).fill(0),
    nodes: tubes.map(t => Array.from({ length: nx + 1 }, () => Array.from({ length: t + 1 }, () => ({ x: 0, y: 0 })))) },
    initialBL: Array(4 * stations).fill(0) } };
  const metadata = streamtubeStationMetadata(checkpoint), used = new Set(), matches = [];
  for (const c of metadata.preferred) { if (c >= 0) assert(!used.has(c)); if (c >= 0) used.add(c); matches.push(c); }
  const remaining = Array.from({ length: metadata.n }, (_, c) => c).filter(c => !used.has(c));
  matches.forEach((c, r) => { if (c < 0) matches[r] = remaining.shift(); });
  assert.equal(remaining.length, 0);
  const graph = matches.map((c, r) => [[c, r % 2 ? -2 : 3]]);
  return { checkpoint, metadata, matches, graph };
}
function csr(graph) {
  const rowPtr = [0], colIndex = [], values = [];
  graph.forEach(row => { row.forEach(([c, v]) => { colIndex.push(c); values.push(v); }); rowPtr.push(values.length); });
  return { n: graph.length, rowPtr: Int32Array.from(rowPtr), colIndex: Int32Array.from(colIndex), values: Float64Array.from(values) };
}
function verify(matrix, result) {
  const { n } = matrix, { Puser: p, Quser: q } = result;
  assert.equal(p.length, n); assert.equal(q.length, n);
  assert.deepEqual([...p].sort((a, b) => a - b), Array.from({ length: n }, (_, k) => k));
  assert.deepEqual([...q].sort((a, b) => a - b), Array.from({ length: n }, (_, k) => k));
  for (let k = 0; k < n; k++) {
    const entries = [];
    for (let j = matrix.rowPtr[p[k]]; j < matrix.rowPtr[p[k] + 1]; j++) if (matrix.colIndex[j] === q[k]) entries.push(matrix.values[j]);
    assert.equal(entries.length, 1); assert(entries[0] !== 0 && Number.isFinite(entries[0]));
  }
}

test('two staggered bodies: actual matched pairs keep every global and dense constraint at the end', () => {
  const f = fixture({ multi: true }), { metadata } = f;
  // Manufacture the structural fact that a Kutta row has no direct Gamma
  // entry: pair it with a physical column and let a local row reach Gamma.
  const kutta = metadata.rows.findIndex(r => r.kind === 'trailingKutta');
  const local = metadata.rows.findIndex(r => r.kind === 'farfieldPressure');
  [f.graph[kutta], f.graph[local]] = [f.graph[local], f.graph[kutta]];
  // Dense global constraint deliberately pairs with a local unknown too.
  const dense = metadata.rows.findIndex(r => r.kind === 'farfieldMatch');
  const local2 = metadata.rows.findIndex(r => r.kind === 'internalPressure');
  [f.graph[dense], f.graph[local2]] = [f.graph[local2], f.graph[dense]];
  const matrix = csr(f.graph), before = structuredClone({ checkpoint: f.checkpoint, matrix });
  const result = createStreamtubeStationOrdering({ matrix, checkpoint: f.checkpoint }); verify(matrix, result);
  assert.deepEqual({ checkpoint: f.checkpoint, matrix }, before);
  const border = result.diagnostics.borderPairs;
  assert.equal(border.length, metadata.layout.globalUnknowns + 1);
  assert.equal(border.filter(p => p.columnMetadata.kind === 'global').length, metadata.layout.globalUnknowns);
  assert.equal(border.filter(p => p.rowMetadata.kind === 'farfieldMatch').length, 3);
  assert.deepEqual([...result.Puser.slice(-border.length)], border.map(p => p.row));
  assert(!border.some(p => p.row === kutta), 'Kutta is not assigned an invented global diagonal.');
  const stationPairs = [...result.Puser.slice(0, -border.length)].map((r, k) => Math.max(metadata.rows[r].i, metadata.columns[result.Quser[k]].i));
  assert(stationPairs.every((s, k) => k === 0 || s >= stationPairs[k - 1]));
  assert(metadata.rows.some(r => r.kind === 'wakeGap'));
});

test('a zero natural diagonal is repaired by an actual alternating path, deterministically', () => {
  const f = fixture(), a = 0, b = 1, ca = f.matches[a], cb = f.matches[b];
  f.graph[a] = [[ca, 4], [cb, -1]]; f.graph[b] = [[cb, 0], [ca, 2]];
  const matrix = csr(f.graph), run = () => createStreamtubeStationOrdering({ matrix, checkpoint: f.checkpoint });
  const first = run(), second = run(); verify(matrix, first);
  assert.deepEqual(first, second); assert.equal(first.diagnostics.exactZeros, 1);
  assert.equal(first.diagnostics.augmentations, 1); assert.equal(first.diagnostics.maximumAugmentingRows, 2);
  const inverseP = [...first.Puser].indexOf(a); assert.equal(first.Quser[inverseP], cb);
  assert.equal(first.Quser[[...first.Puser].indexOf(b)], ca);
});

test('full matching retains long-range cumulative wake couplings without forcing their diagonal', () => {
  const f = fixture({ multi: true, independent: false });
  const wake = f.metadata.rows.findIndex(r => r.kind === 'bl-energy' && r.station !== undefined && r.i === 9);
  const upstream = f.metadata.columns.findIndex(c => c.kind === 'bl-delta' && c.i === 3);
  assert(wake >= 0 && upstream >= 0); f.graph[wake].push([upstream, 1e-40]);
  const matrix = csr(f.graph), before = structuredClone(matrix);
  const result = createStreamtubeStationOrdering({ matrix, checkpoint: f.checkpoint }); verify(matrix, result);
  assert.deepEqual(matrix, before); assert.equal(result.diagnostics.numericalNonzeros, matrix.n + 1);
  assert.equal(result.Quser[[...result.Puser].indexOf(wake)], f.matches[wake]);
  assert.equal(result.diagnostics.factorizationPerformed, false);
});

test('structural deficiency fails explicitly without perturbation or partial permutations', () => {
  const f = fixture(); f.graph[1] = [[f.matches[0], 2]];
  const matrix = csr(f.graph), before = structuredClone(matrix);
  assert.throws(() => createStreamtubeStationOrdering({ matrix, checkpoint: f.checkpoint }), error => {
    assert.equal(error.code, 'STREAMTUBE_STATION_MATCHING_FAILED');
    assert.equal(error.diagnostics.unmatchedRows.length, 1); assert.equal(error.diagnostics.unmatchedColumns.length, 1); return true;
  });
  assert.deepEqual(matrix, before);
});

test('malformed CSR, duplicate entries and inconsistent checkpoint dimensions are rejected', () => {
  const f = fixture(), matrix = csr(f.graph), run = (a, cp = f.checkpoint) => createStreamtubeStationOrdering({ matrix: a, checkpoint: cp });
  for (const mutation of [a => { a.values[0] = NaN; }, a => { a.colIndex[0] = a.n; },
    a => { a.rowPtr[1] = -1; }, a => { a.n--; }]) {
    const bad = structuredClone(matrix); mutation(bad); assert.throws(() => run(bad));
  }
  f.graph[0].push(f.graph[0][0]); assert.throws(() => run(csr(f.graph)), /unique/);
  const cp = structuredClone(f.checkpoint); cp.restart.initialBL.pop(); assert.throws(() => run(matrix, cp), /BL packed length/);
});
