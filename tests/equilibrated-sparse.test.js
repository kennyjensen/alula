import test from 'node:test';
import assert from 'node:assert/strict';
import { solveEquilibratedSparse } from '../src/numerics/tests/equilibrated-sparse.js';
import { sparseMatrix, sparseProduct, sparseDense } from '../src/numerics/sparse.js';
import { solveLinear } from '../src/numerics/linear.js';

test('column equilibration preserves an independently solved nonsymmetric system across 24 orders of variable units', () => {
  const n = 37, base = sparseMatrix(Array.from({ length: n }, (_, i) => [i, (i + 1) % n, (i + 7) % n]));
  for (let i = 0; i < n; i++) for (let p = base.rowPtr[i]; p < base.rowPtr[i + 1]; p++)
    base.values[p] = base.colIndex[p] === i ? 3 : Math.sin(.7 * i + base.colIndex[p]);
  const rhs = Float64Array.from({ length: n }, (_, i) => Math.cos(i)), expected = solveLinear(sparseDense(base), rhs);
  const units = Float64Array.from({ length: n }, (_, i) => 10 ** (i % 25 - 12));
  const matrix = { ...base, values: base.values.map((v, p) => v * units[base.colIndex[p]]) };
  const original = structuredClone(matrix), beforeRhs = rhs.slice();
  const r = solveEquilibratedSparse(matrix, rhs);
  assert.ok(r.relativeResidual < 1e-12 && r.scaledRelativeResidual < 1e-12);
  for (let i = 0; i < n; i++) assert.ok(Math.abs(r.x[i] * units[i] - expected[i]) < 2e-14);
  const product = sparseProduct(matrix, r.x);
  assert.ok(Math.hypot(...product.map((v, i) => v - rhs[i])) / Math.hypot(...rhs) < 1e-12);
  assert.deepEqual(matrix, original); assert.deepEqual(rhs, beforeRhs);
});

test('equilibrated WASM solve rejects singular/nonfinite inputs and invalid controls, then recovers', () => {
  const a = sparseMatrix([[0, 1], [0, 1]]), b = Float64Array.of(1, 0);
  a.values.set([1, 0, 2, 0]); assert.throws(() => solveEquilibratedSparse(a, b), /zero matrix column/);
  a.values.set([1, 1, 2, 2]); assert.throws(() => solveEquilibratedSparse(a, b), /factorization failed/);
  a.values.set([0, 2, 3, 0]);
  for (const options of [{ tolerance: NaN }, { tolerance: 0 }, { maxRefinements: -1 }, { ordering: 'bad' }])
    assert.throws(() => solveEquilibratedSparse(a, b, options), /controls/);
  assert.throws(() => solveEquilibratedSparse(a, Float64Array.of(NaN, 0)), /right-hand side/);
  a.values[0] = Infinity; assert.throws(() => solveEquilibratedSparse(a, b), /matrix entry/);
  a.values.set([0, 2, 3, 0]);
  assert.deepEqual(Array.from(solveEquilibratedSparse(a, Float64Array.of(4, 9)).x), [3, 2]);
  assert.deepEqual(Array.from(solveEquilibratedSparse(a, new Float64Array(2)).x), [0, 0]);
});
