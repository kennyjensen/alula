// SPDX-License-Identifier: GPL-2.0-or-later
// Sparse linear systems and topology only. No flow/Jacobian constructors.
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodyLayout } from '../src/euler/streamtube-body-layout.js';
import { createStreamtubeStationOrdering, streamtubeStationMetadata } from '../src/euler/streamtube-station-ordering.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { solveCoupledLinearSystem } from '../src/euler/streamtube-coupled-linear-solve.js';

function fixture() {
  const nx = 7, tubes = [2, 1, 2], bodies = [{ leadingIndex: 2, trailingIndex: 4 },
    { leadingIndex: 3, trailingIndex: 5, roundedLeadingEdge: false }];
  const layout = createStreamtubeBodyLayout({ segments: nx, tubes, bodies, independentWakeBanks: true });
  const stations = [];
  bodies.forEach((b, body) => {
    for (const side of ['upper', 'lower']) for (let i = b.leadingIndex + 1; i <= b.trailingIndex; i++)
      stations.push({ id: stations.length, kind: 'surface', body, side, i });
    for (let i = b.trailingIndex; i <= nx; i++) stations.push({ id: stations.length, kind: 'wake', body, i });
  });
  const checkpoint = { version: 1, restart: { input: { bodies, wakeGeometry: 'independent-banks' },
    initialEuler: { x: Array(layout.n).fill(0), nodes: tubes.map(t => Array.from({ length: nx + 1 },
      () => Array.from({ length: t + 1 }, () => ({ x: 0, y: 0 })))) }, initialBL: Array(stations.length * 4).fill(0) } };
  const metadata = streamtubeStationMetadata(checkpoint), used = new Set([...metadata.preferred].filter(c => c >= 0));
  const remaining = Array.from({ length: metadata.n }, (_, c) => c).filter(c => !used.has(c));
  const matched = [...metadata.preferred].map(c => c >= 0 ? c : remaining.shift());
  const graph = matched.map((c, i) => [[c, 2 + (i % 5) / 7]]);
  // Nonsymmetric block forces a numerical pivot; the matrix is not diagonal.
  graph[0] = [[matched[0], 1e-12], [matched[1], 2]];
  graph[1] = [[matched[0], 3], [matched[1], 1]];
  graph.at(-2).push([layout.n + 2, .025]); // Nonlocal wake-to-surface coupling.
  const rowPtr = [0], colIndex = [], values = [];
  graph.forEach(row => {
    row.sort(([a], [b]) => a - b).forEach(([c, value]) => { colIndex.push(c); values.push(value); });
    rowPtr.push(values.length);
  });
  const matrix = { n: metadata.n, rowPtr: Int32Array.from(rowPtr), colIndex: Int32Array.from(colIndex), values: Float64Array.from(values) };
  const exact = Float64Array.from({ length: matrix.n }, (_, i) => Math.cos(i / 13) + i / 100);
  const rhs = Float64Array.from({ length: matrix.n }, (_, row) => {
    let sum = 0;
    for (let k = rowPtr[row]; k < rowPtr[row + 1]; k++) sum += values[k] * exact[colIndex[k]];
    return sum;
  });
  return { matrix, rhs, exact, layout, stations };
}
const given = (ordering, pivotTolerance = .001) => ({ ordering: 'given', rowPermutation: ordering.Puser,
  columnPermutation: ordering.Quser, btf: false, pivotTolerance, pivotFallback: false });
const controls = f => ({ layout: f.layout, stations: f.stations });
function verify(f, solved) {
  assert(solved.relativeResidual <= 1e-10);
  assert(Math.max(...solved.x.map((v, i) => Math.abs(v - f.exact[i]))) < 1e-12);
}
function timing(policy) {
  for (const value of Object.values(policy.timings)) assert(Number.isFinite(value) && value >= 0);
  const t = policy.timings;
  assert(t.totalMilliseconds >= t.matchingMilliseconds + t.stationSolveMilliseconds + t.autoSolveMilliseconds);
}

let isolatedId = 0;
async function injected(overrides = {}) {
  const key = `__coupledLinearPolicy${++isolatedId}`;
  globalThis[key] = { solveSparseDirect, createStreamtubeStationOrdering, ...overrides };
  const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-linear-solve.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  try {
    return (await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`)).solveCoupledLinearSystem;
  } finally { delete globalThis[key]; }
}

test('legacy auto retains the complete existing KLU result and does not require station metadata', () => {
  const f = fixture(), before = structuredClone({ matrix: f.matrix, rhs: f.rhs });
  for (const extra of [{}, { preferredOrdering: 'colamd', pivotTolerance: 1 }]) {
    const expected = solveSparseDirect(f.matrix, f.rhs, { preferredOrdering: 'amd', pivotTolerance: .001, ...extra });
    const result = solveCoupledLinearSystem(f.matrix, f.rhs, extra);
    assert.deepEqual(result, { linear: expected, ordering: null }); verify(f, result.linear);
  }
  assert.deepEqual({ matrix: f.matrix, rhs: f.rhs }, before);
});

test('explicit station retains matching, options and full linear payload without a fallback policy', () => {
  const f = fixture(), order = createStreamtubeStationOrdering({ matrix: f.matrix, ...controls(f) });
  for (const pivotTolerance of [.001, 1]) {
    const expected = solveSparseDirect(f.matrix, f.rhs, given(order, pivotTolerance));
    const actual = solveCoupledLinearSystem(f.matrix, f.rhs, { ...controls(f), mode: 'station', pivotTolerance });
    assert.deepEqual(actual, { linear: expected, ordering: order }); verify(f, actual.linear);
  }
});

test('station-auto accepts the original full system with numerical pivoting and an unchanged accuracy gate', () => {
  const f = fixture(), before = structuredClone({ matrix: f.matrix, rhs: f.rhs, stations: f.stations });
  const layoutBefore = JSON.stringify(f.layout);
  const legacy = solveCoupledLinearSystem(f.matrix, f.rhs, { ...controls(f), mode: 'station' });
  const result = solveCoupledLinearSystem(f.matrix, f.rhs, { ...controls(f), mode: 'station-auto' });
  assert.deepEqual(result.linear, legacy.linear); assert.deepEqual(result.ordering, legacy.ordering);
  assert.equal(result.linear.attempts.length, 1); assert.equal(result.linear.btf, false);
  assert.equal(result.stationPolicy.attempted, true); assert.equal(result.stationPolicy.accepted, true);
  assert.equal(result.stationPolicy.fallback, false); assert.equal(result.stationPolicy.recommendation, 'station-auto');
  assert.deepEqual(result.stationPolicy.selected, { ordering: 'given', pivotTolerance: .001 });
  assert.equal(result.stationPolicy.timings.autoSolveMilliseconds, 0); timing(result.stationPolicy); verify(f, result.linear);
  assert.deepEqual({ matrix: f.matrix, rhs: f.rhs, stations: f.stations }, before);
  assert.equal(JSON.stringify(f.layout), layoutBefore);
});

test('typed matching rejection falls back once, preserves the original matrix and recommends skipping further attempts', async () => {
  const f = fixture(), seen = [], rejection = Object.assign(new Error('Controlled incomplete matching'), {
    code: 'STREAMTUBE_STATION_MATCHING_FAILED', diagnostics: { unmatchedRows: [{ index: 3 }] } });
  let matches = 0;
  const solve = await injected({ createStreamtubeStationOrdering() { matches++; throw rejection; },
    solveSparseDirect(a, b, options) { assert.equal(a, f.matrix); assert.equal(b, f.rhs); seen.push(options); return solveSparseDirect(a, b, options); } });
  const result = solve(f.matrix, f.rhs, { ...controls(f), mode: 'station-auto', preferredOrdering: 'colamd', pivotTolerance: 1 });
  assert.equal(matches, 1); assert.deepEqual(seen, [{ preferredOrdering: 'colamd', pivotTolerance: 1 }]);
  assert.equal(result.ordering, null); verify(f, result.linear); timing(result.stationPolicy);
  assert.equal(result.stationPolicy.accepted, false); assert.equal(result.stationPolicy.fallback, true);
  assert.equal(result.stationPolicy.recommendation, 'auto');
  assert.equal(result.stationPolicy.stationFailure.stage, 'matching');
  rejection.diagnostics.unmatchedRows[0].index = 99;
  assert.equal(result.stationPolicy.stationFailure.diagnostics.unmatchedRows[0].index, 3);
  const next = solve(f.matrix, f.rhs, { mode: 'station-auto', stationEnabled: false, stationSkipReason: 'retained-policy' });
  assert.equal(matches, 1); assert.equal(seen.length, 2);
  assert.equal(next.stationPolicy.attempted, false); assert.equal(next.stationPolicy.fallback, false);
  assert.equal(next.stationPolicy.skippedReason, 'retained-policy'); assert.equal(next.stationPolicy.recommendation, 'auto');
  assert.equal(next.stationPolicy.timings.matchingMilliseconds, 0); assert.equal(next.stationPolicy.timings.stationSolveMilliseconds, 0);
});

test('typed station accuracy/factor failures fall back to the existing automatic options and actual pivot', async () => {
  const f = fixture(), before = structuredClone({ matrix: f.matrix, rhs: f.rhs });
  for (const rejection of [Object.assign(new Error('Controlled factor singularity'), { code: 'KLU_GIVEN_FACTORIZATION', status: 1 }),
    Object.assign(new Error('Controlled cleaned-up factor allocation failure'), { code: 'KLU_GIVEN_FACTORIZATION', status: -2 }),
    Object.assign(new Error('Controlled residual rejection'), { code: 'KLU_RESIDUAL_LIMIT', relativeResidual: 2e-10,
      diagnostics: { tolerance: 1e-10, residualEvaluation: 'compensated-original-system' }, attempts: [{ ordering: 'given', relativeResidual: 2e-10 }] })]) {
    const seen = [], solve = await injected({ solveSparseDirect(a, b, options) {
      assert.equal(a, f.matrix); assert.equal(b, f.rhs); seen.push(options);
      if (options.ordering === 'given') throw rejection;
      return solveSparseDirect(a, b, options);
    } });
    const result = solve(f.matrix, f.rhs, { ...controls(f), mode: 'station-auto', pivotTolerance: 1 });
    assert.equal(seen.length, 2); assert.equal(seen[0].pivotTolerance, 1);
    assert.deepEqual(seen[1], { preferredOrdering: 'amd', pivotTolerance: 1 });
    assert.equal(result.ordering, null); assert.equal(result.stationPolicy.stationFailure.code, rejection.code);
    assert.deepEqual(result.stationPolicy.selected, { ordering: result.linear.ordering, pivotTolerance: result.linear.pivotTolerance });
    assert.equal(result.stationPolicy.stationFailure.stage, 'station-solve'); verify(f, result.linear); timing(result.stationPolicy);
  }
  assert.deepEqual({ matrix: f.matrix, rhs: f.rhs }, before);
});

test('matching bugs, resource errors, invalid controls, and cancellations are never converted to fallback solves', async () => {
  const f = fixture();
  const failures = [new Error('Unexpected implementation failure'), new RangeError('Invalid matrix'),
    new Error('Unclassified allocation failure'),
    ...[-3, -4].map(status => Object.assign(new Error('Native resource or contract failure'), { code: 'KLU_GIVEN_FACTORIZATION', status })),
    Object.assign(new Error('Aborted'), { name: 'AbortError', code: 'KLU_RESIDUAL_LIMIT' }),
    Object.assign(new Error('Stopped'), { code: 'ABORT_ERR' }),
    Object.assign(new Error('Unclassified solve failure'), { code: 'KLU_GIVEN_SOLVE' })];
  for (const rejection of failures) {
    let factors = 0;
    const solve = await injected({ solveSparseDirect() { factors++; throw rejection; } });
    assert.throws(() => solve(f.matrix, f.rhs, { ...controls(f), mode: 'station-auto' }), error => error === rejection);
    assert.equal(factors, 1); assert.equal(rejection.stationPolicy, undefined);
  }
  for (const rejection of [new Error('Station ordering: malformed metadata'),
    Object.assign(new Error('Cancelled matching'), { name: 'AbortError', code: 'STREAMTUBE_STATION_MATCHING_FAILED' })]) {
    let calls = 0;
    const solve = await injected({ createStreamtubeStationOrdering() { throw rejection; }, solveSparseDirect() { calls++; } });
    assert.throws(() => solve(f.matrix, f.rhs, { ...controls(f), mode: 'station-auto' }), error => error === rejection);
    assert.equal(calls, 0);
  }
});

test('terminal automatic numerical failure preserves both attempts, original code and original error identity', async () => {
  const f = fixture(), first = Object.assign(new Error('Station certificate rejected'), { code: 'KLU_RESIDUAL_LIMIT',
    relativeResidual: 2e-10, attempts: [{ ordering: 'given', pivotTolerance: .001, relativeResidual: 2e-10 }] });
  const last = Object.assign(new Error('Automatic certificate rejected'), { code: 'KLU_RESIDUAL_LIMIT', relativeResidual: 3e-10,
    diagnostics: { tolerance: 1e-10, residualEvaluation: 'compensated-original-system' },
    attempts: [{ ordering: 'colamd', pivotTolerance: 1, relativeResidual: 3e-10 }] });
  const solve = await injected({ solveSparseDirect(a, b, options) { throw options.ordering === 'given' ? first : last; } });
  assert.throws(() => solve(f.matrix, f.rhs, { ...controls(f), mode: 'station-auto' }), error => {
    assert.equal(error, last); assert.equal(error.code, 'KLU_RESIDUAL_LIMIT'); assert.equal(error.relativeResidual, 3e-10);
    assert.equal(error.diagnostics.tolerance, 1e-10); assert.equal(error.stationPolicy.stationFailure.relativeResidual, 2e-10);
    assert.equal(error.stationPolicy.autoFailure.relativeResidual, 3e-10);
    assert.deepEqual(error.stationPolicy.autoFailure.attempts, last.attempts);
    assert.equal(error.stationPolicy.recommendation, 'auto'); timing(error.stationPolicy);
    return true;
  });
  assert.equal(first.stationPolicy, undefined);
});

test('cancellation in the fallback path propagates unchanged without another attempt', async () => {
  const f = fixture(), first = Object.assign(new Error('Controlled matching failure'), { code: 'STREAMTUBE_STATION_MATCHING_FAILED' });
  const abort = Object.assign(new Error('Cancelled during automatic solve'), { name: 'AbortError', code: 'ABORT_ERR' });
  let calls = 0;
  const solve = await injected({ createStreamtubeStationOrdering() { throw first; },
    solveSparseDirect() { calls++; throw abort; } });
  assert.throws(() => solve(f.matrix, f.rhs, { ...controls(f), mode: 'station-auto' }), error => error === abort);
  assert.equal(calls, 1); assert.equal(abort.stationPolicy, undefined); assert.equal(abort.diagnostics, undefined);
});

test('real deficient pattern stays a failure after bounded automatic KLU attempts', () => {
  const f = fixture(), bad = structuredClone(f.matrix);
  // Two single-entry rows now share one column; another column is unmatched.
  assert.equal(bad.rowPtr[3] - bad.rowPtr[2], 1); assert.equal(bad.rowPtr[4] - bad.rowPtr[3], 1);
  bad.colIndex[bad.rowPtr[3]] = bad.colIndex[bad.rowPtr[2]];
  const before = structuredClone(bad);
  assert.throws(() => solveCoupledLinearSystem(bad, f.rhs, { ...controls(f), mode: 'station-auto' }), error => {
    assert.equal(error.stationPolicy.stationFailure.code, 'STREAMTUBE_STATION_MATCHING_FAILED');
    assert(error.attempts.every(attempt => attempt.status === 1)); assert(error.attempts.length <= 4);
    assert.equal(error.stationPolicy.recommendation, 'auto'); return true;
  });
  assert.deepEqual(bad, before);
});

test('malformed actual topology and CSR remain terminal; explicit station never falls back', async () => {
  const f = fixture();
  assert.throws(() => solveCoupledLinearSystem(f.matrix, f.rhs, { ...controls(f), stations: f.stations.slice(1), mode: 'station-auto' }), /ordered ids/);
  const matrix = structuredClone(f.matrix); matrix.values[0] = NaN;
  assert.throws(() => solveCoupledLinearSystem(matrix, f.rhs, { ...controls(f), mode: 'station-auto' }), /finite/);
  for (const options of [{ mode: 'other' }, { mode: 'station-auto', stationEnabled: 1 },
    { mode: 'station-auto', pivotTolerance: 0 }, { mode: 'station-auto', preferredOrdering: 'given' },
    { mode: 'station-auto', stationEnabled: false, stationSkipReason: '' }])
    assert.throws(() => solveCoupledLinearSystem(f.matrix, f.rhs, options), /policy|controls/);
  let calls = 0;
  const rejection = Object.assign(new Error('Controlled explicit failure'), { code: 'KLU_RESIDUAL_LIMIT' });
  const solve = await injected({ solveSparseDirect() { calls++; throw rejection; } });
  assert.throws(() => solve(f.matrix, f.rhs, { ...controls(f), mode: 'station' }), error => error === rejection);
  assert.equal(calls, 1);
});
