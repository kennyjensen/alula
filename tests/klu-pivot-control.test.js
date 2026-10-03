import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { sparseProduct } from '../src/numerics/sparse.js';

test('maximum partial pivoting resolves a retained coupled AMD failure without changing its matrix or residual bound', t => {
  const f = JSON.parse(fs.readFileSync('tests/fixtures/klu-refined-coupled.json'));
  const a = { ...f.matrix, rowPtr: Int32Array.from(f.matrix.rowPtr), colIndex: Int32Array.from(f.matrix.colIndex), values: Float64Array.from(f.matrix.values) };
  const b = Float64Array.from(f.rhs), matrixBefore = a.values.slice(), rhsBefore = b.slice();
  assert.throws(() => solveSparseDirect(a, b, { ordering: 'amd', pivotFallback: false }), e => e.code === 'KLU_RESIDUAL_LIMIT');
  const r = solveSparseDirect(a, b, { ordering: 'amd' }), ax = sparseProduct(a, r.x);
  assert.ok(r.relativeResidual < 1e-10); assert.equal(r.ordering, 'amd'); assert.equal(r.attempts.length, 2);
  assert.equal(r.attempts[0].pivotTolerance, .001); assert.ok(r.attempts[0].relativeResidual > 1e-10);
  assert.equal(r.pivotTolerance, 1);
  assert.ok(Math.hypot(...b.map((v, i) => v - ax[i])) / Math.hypot(...b) < 1e-10);
  const expected = Float64Array.from({ length: a.n }, (_, i) => Math.sin(.37 * i));
  const known = solveSparseDirect(a, sparseProduct(a, expected), { ordering: 'amd', pivotTolerance: 1 });
  const error = Math.max(...known.x.map((v, i) => Math.abs(v - expected[i]))); assert.ok(error < 2e-7);
  assert.deepEqual(a.values, matrixBefore); assert.deepEqual(b, rhsBefore);
  t.diagnostic(JSON.stringify({ unknowns: a.n, attempts: r.attempts, linearResidual: r.relativeResidual, manufacturedForwardError: error }));
});
