import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodyLayout } from '../src/euler/streamtube-body-layout.js';
import { createStreamtubeStationOrdering, streamtubeStationMetadata } from '../src/euler/streamtube-station-ordering.js';
import { createStreamtubeStationOrdering as prototypeOrdering } from '../scripts/validation/streamtube-station-ordering.js';

// Only topology and sparse arrays: no Euler/BL constructors or evaluations.
function fixture({ independent = true, multi = true } = {}) {
  const nx = 9, tubes = multi ? [2, 1, 3] : [2, 2];
  const bodies = (multi ? [[2, 5], [3, 6]] : [[2, 5]]).map(([leadingIndex, trailingIndex], body) => ({
    leadingIndex, trailingIndex, roundedLeadingEdge: body === 0,
    trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 4 } }));
  const layout = createStreamtubeBodyLayout({ segments: nx, tubes, bodies, independentWakeBanks: independent });
  const stations = [];
  for (let body = 0; body < bodies.length; body++) {
    const b = bodies[body];
    for (const side of ['upper', 'lower']) for (let i = b.leadingIndex + 1; i <= b.trailingIndex; i++) {
      const k = i - b.leadingIndex;
      stations.push({ kind: 'surface', body, side, i, k, id: stations.length,
        regime: k === 1 ? 'similarity' : k === 2 ? 'transition' : 'turbulent' });
    }
    for (let i = b.trailingIndex; i <= nx; i++) stations.push({ kind: 'wake', body, i,
      k: i - b.trailingIndex, id: stations.length, regime: i === b.trailingIndex ? 'trailing-edge' : 'wake' });
  }
  const checkpoint = { version: 1, restart: { input: { bodies, primaryBody: 0, flowModel: 'compressible',
    wakeGeometry: independent ? 'independent-banks' : 'centerline' }, initialEuler: { x: Array(layout.n).fill(0),
    nodes: tubes.map(t => Array.from({ length: nx + 1 }, () => Array.from({ length: t + 1 }, () => ({ x: 0, y: 0 })))) },
    initialBL: Array(4 * stations.length).fill(0) } };
  const metadata = streamtubeStationMetadata(checkpoint), used = new Set([...metadata.preferred].filter(c => c >= 0));
  const available = Array.from({ length: metadata.n }, (_, i) => i).filter(i => !used.has(i));
  const matched = [...metadata.preferred].map(c => c >= 0 ? c : available.shift());
  const graph = matched.map(c => [[c, 2]]);
  // A real augmenting-path repair, not merely a diagonal input.
  graph[0].push([matched[1], -1]); graph[1] = [[matched[1], 0], [matched[0], 3]];
  // Cumulative wake-arc coupling remains in the matrix, independent of the
  // station enumeration and physical finite-base gap/phase fields.
  const wakeRow = layout.n + 4 * (stations.length - 1) + 2;
  graph[wakeRow].push([layout.n + 2, 1e-40]);
  const rowPtr = [0], colIndex = [], values = [];
  graph.forEach(row => { row.forEach(([c, v]) => { colIndex.push(c); values.push(v); }); rowPtr.push(values.length); });
  const matrix = { n: metadata.n, rowPtr: Int32Array.from(rowPtr), colIndex: Int32Array.from(colIndex), values: Float64Array.from(values) };
  return { checkpoint, layout, stations, matrix };
}

test('runtime checkpoint matching is exactly the frozen prototype, and live metadata returns the same P/Q', () => {
  for (const options of [{}, { independent: false }, { multi: false }]) {
    const f = fixture(options), source = JSON.stringify({ checkpoint: f.checkpoint, layout: f.layout, stations: f.stations }), matrixBefore = structuredClone(f.matrix);
    const old = prototypeOrdering({ matrix: f.matrix, checkpoint: f.checkpoint });
    const saved = createStreamtubeStationOrdering({ matrix: f.matrix, checkpoint: f.checkpoint });
    const live = createStreamtubeStationOrdering({ matrix: f.matrix, layout: f.layout, stations: f.stations });
    assert.deepEqual(saved, old); assert.deepEqual(live, saved);
    assert.equal(live.diagnostics.augmentations, 1); assert.equal(live.diagnostics.allSelectedDiagonalsFiniteNonzero, true);
    assert.equal(JSON.stringify({ checkpoint: f.checkpoint, layout: f.layout, stations: f.stations }), source);
    assert.deepEqual(f.matrix, matrixBefore);
  }
});

test('live metadata never reconstructs a layout or reads a synthetic checkpoint', async () => {
  const source = fs.readFileSync(new URL('../src/euler/streamtube-station-ordering.js', import.meta.url), 'utf8');
  const imported = "import { createStreamtubeBodyLayout } from './streamtube-body-layout.js';";
  assert.equal(source.split(imported).length, 2);
  const isolated = await import(`data:text/javascript;base64,${Buffer.from(source.replace(imported,
    "const createStreamtubeBodyLayout = () => { throw new Error('unexpected layout construction'); };")).toString('base64')}`);
  const f = fixture();
  assert.deepEqual(isolated.createStreamtubeStationOrdering({ matrix: f.matrix, layout: f.layout, stations: f.stations }),
    createStreamtubeStationOrdering({ matrix: f.matrix, checkpoint: f.checkpoint }));
  assert.throws(() => isolated.createStreamtubeStationOrdering({ matrix: f.matrix, checkpoint: f.checkpoint }), /unexpected layout construction/);
});

test('live finite-base wake station order is bound while auxiliary phase metadata is irrelevant', () => {
  const f = fixture(), call = stations => createStreamtubeStationOrdering({ matrix: f.matrix, layout: f.layout, stations });
  const baseline = call(f.stations);
  const changedPhase = f.stations.map(s => ({ ...s, regime: 'diagnostic-only', theta: .01, deltaStar: .02, aux: 123, s: s.id }));
  assert.deepEqual(call(changedPhase), baseline);
  for (const mutate of [s => s.pop(), s => { s[0].id = 1; }, s => { [s[1], s[2]] = [s[2], s[1]]; },
    s => { delete s[2]; }, s => { s.find(a => a.kind === 'wake').i++; },
    s => { s.find(a => a.kind === 'wake').side = 'upper'; }, s => { s[0].k++; }]) {
    const broken = structuredClone(f.stations); mutate(broken); assert.throws(() => call(broken), /ordered ids/);
  }
});

test('checkpoint and live forms are mutually exclusive and matrix dimensions keep the four-slot contract', () => {
  const f = fixture();
  for (const args of [{}, { layout: f.layout }, { stations: f.stations },
    { checkpoint: f.checkpoint, layout: f.layout, stations: f.stations }, { checkpoint: null, layout: f.layout, stations: f.stations }])
    assert.throws(() => createStreamtubeStationOrdering({ matrix: f.matrix, ...args }), /exclusively/);
  assert.throws(() => createStreamtubeStationOrdering({ matrix: { ...f.matrix, n: f.matrix.n - 1 }, layout: f.layout, stations: f.stations }), /square CSR/);
  assert.throws(() => createStreamtubeStationOrdering({ matrix: f.matrix, layout: { ...f.layout, densityUnknowns: false }, stations: f.stations }), /compressible/);
});
