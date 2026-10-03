// Compare identical frozen Euler systems; no operating-point or solver changes.
// Usage: node scripts/validation/benchmark-euler-ordering.js [checkpoint.json]
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../../src/euler/streamtube-body.js';
import { streamtubeEquationOrder } from '../../src/euler/streamtube-body-layout.js';
import { solveSparseDirect, solveSparseDirectAligned } from '../../src/numerics/klu.js';
import { sparseProduct } from '../../src/numerics/sparse.js';

const source = process.argv[2] ?? 'docs/solver-reliability/rae-inviscid-audit/fixtures/rae64x9-harmonic-cold-euler-root.json';
const checkpoint = JSON.parse(fs.readFileSync(source, 'utf8'));
const system = createStreamtubeBodySystem(checkpoint.input);
const state = system.adoptGeometry(Float64Array.from(checkpoint.initialEuler.x), checkpoint.initialEuler.nodes);
const matrix = system.jacobian(state, { sparse: true });
const rhs = system.evaluate(state).residual.map(v => -v);
const order = streamtubeEquationOrder(system.layout), results = [];
const run = (method, b) => method === 'aligned'
  ? solveSparseDirectAligned(matrix, b, order) : solveSparseDirect(matrix, b);
// Warm both routes before reporting, then alternate to reduce timing bias.
for (const method of ['original', 'aligned']) run(method, rhs);
let original;
for (const method of ['original', 'aligned', 'original', 'aligned']) {
  const start = performance.now(), result = run(method, rhs);
  const milliseconds = performance.now() - start;
  assert.ok(result.relativeResidual <= 1e-10);
  if (method === 'original') original = result.x;
  let directionDifference = 0, directionScale = 0;
  for (let i = 0; i < matrix.n; i++) {
    directionDifference = Math.max(directionDifference, Math.abs(result.x[i] - original[i]));
    directionScale = Math.max(directionScale, Math.abs(original[i]));
  }
  results.push({ method, milliseconds, factorNonzeros: result.factorNonzeros,
    relativeResidual: result.relativeResidual, relativeDirectionDifference: directionDifference / (directionScale || 1),
    equationOrdering: result.equationOrdering ?? 'original' });
}
// A nonzero manufactured solution also checks correctness at an almost-zero root RHS.
const exact = Float64Array.from({ length: matrix.n }, (_, i) => Math.sin(.37 * i));
const manufactured = run('aligned', sparseProduct(matrix, exact));
let forwardError = 0;
for (let i = 0; i < matrix.n; i++) forwardError = Math.max(forwardError, Math.abs(manufactured.x[i] - exact[i]));
assert.ok(forwardError < 1e-6, `Manufactured forward error ${forwardError}`);
console.log(JSON.stringify({ source, unknowns: matrix.n, matrixNonzeros: matrix.values.length,
  manufacturedForwardError: forwardError, results }, null, 2));
