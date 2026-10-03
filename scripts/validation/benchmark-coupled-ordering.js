// SPDX-License-Identifier: GPL-2.0-or-later
// node scripts/validation/benchmark-coupled-ordering.js checkpoint.json
// Alternating solves of ONE frozen full Jacobian, with unchanged accuracy.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import { fileURLToPath } from 'node:url';
import { createCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';
import { createStreamtubeStationOrdering } from '../../src/euler/streamtube-station-ordering.js';
import { solveSparseDirect, solveSparseDirectAligned } from '../../src/numerics/klu.js';
import { sparseProduct } from '../../src/numerics/sparse.js';

const source = process.argv[2] ?? fileURLToPath(new URL('../../docs/solver-reliability/nlr64x24-subsonic/performance/late-coupled-checkpoint.json.gz', import.meta.url));
const bytes = fs.readFileSync(source);
const checkpoint = JSON.parse(source.endsWith('.gz') ? gunzipSync(bytes) : bytes), r = checkpoint.restart;
let started = performance.now();
const system = createCoupledStreamtubeBody(r.input, { ...r.options, initialEuler: r.initialEuler, initialBL: r.initialBL });
const constructionMilliseconds = performance.now() - started;
started = performance.now();
const matrix = system.jacobian(system.initial), rhs = system.evaluate(system.initial).residual.map(v => -v);
const jacobianAndResidualMilliseconds = performance.now() - started;
started = performance.now();
const matching = createStreamtubeStationOrdering({ matrix, layout: system.euler.layout, stations: system.bl.stations });
const rows = new Int32Array(matrix.n);
for (let k = 0; k < matrix.n; k++) rows[matching.Quser[k]] = matching.Puser[k];
const matchingMilliseconds = performance.now() - started;
const solve = (method, b) => method === 'station'
  ? solveSparseDirect(matrix, b, { ordering: 'given', rowPermutation: matching.Puser, columnPermutation: matching.Quser, btf: false })
  : solveSparseDirectAligned(matrix, b, rows);
const results = []; let reference;
for (const method of ['station', 'aligned', 'aligned', 'station']) {
  started = performance.now(); const result = solve(method, rhs), milliseconds = performance.now() - started;
  reference ??= result.x;
  let difference = 0, scale = 0;
  for (let i = 0; i < matrix.n; i++) {
    difference = Math.max(difference, Math.abs(result.x[i] - reference[i]));
    scale = Math.max(scale, Math.abs(reference[i]));
  }
  assert(result.relativeResidual <= 1e-10);
  assert(difference / (scale || 1) < 1e-6);
  results.push({ method, milliseconds, factorNonzeros: result.factorNonzeros,
    relativeResidual: result.relativeResidual, relativeDirectionDifference: difference / (scale || 1) });
}
const exact = Float64Array.from({ length: matrix.n }, (_, i) => Math.sin(.37 * i));
const manufactured = solve('aligned', sparseProduct(matrix, exact));
let manufacturedForwardError = 0;
for (let i = 0; i < matrix.n; i++) manufacturedForwardError = Math.max(manufacturedForwardError, Math.abs(manufactured.x[i] - exact[i]));
assert(manufactured.relativeResidual <= 1e-10);
assert(manufacturedForwardError < 1e-5);
const hash = createHash('sha256');
for (const data of [matrix.rowPtr, matrix.colIndex, matrix.values, Float64Array.from(rhs)])
  hash.update(Buffer.from(data.buffer, data.byteOffset, data.byteLength));
console.log(JSON.stringify({ source, matrixAndRhsHash: hash.digest('hex'), unknowns: matrix.n, nonzeros: matrix.values.length,
  constructionMilliseconds, jacobianAndResidualMilliseconds, matchingMilliseconds, manufacturedForwardError, results }, null, 2));
