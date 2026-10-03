// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { solveSparseDirectAlignedMany } from '../src/numerics/klu.js';
import { sparseMatrix, sparseProduct } from '../src/numerics/sparse.js';

test('shared aligned factors preserve each original equation and all inputs', () => {
  const n = 23, a = sparseMatrix(Array.from({ length: n }, (_, i) => [i, (i + 1) % n, (i + 7) % n]));
  for (let i = 0; i < n; i++) for (let k = a.rowPtr[i]; k < a.rowPtr[i + 1]; k++)
    a.values[k] = (a.colIndex[k] === i ? 2 : Math.sin(i + a.colIndex[k])) * 10 ** (i % 9 - 4);
  const exact = [Float64Array.from({ length: n }, (_, i) => Math.sin(i)),
    Float64Array.from({ length: n }, (_, i) => 1e-7 * Math.cos(i)), new Float64Array(n)];
  const rhs = exact.map(x => sparseProduct(a, x)), before = structuredClone({ a, rhs });
  const order = Array.from({ length: n }, (_, i) => n - i - 1);
  for (let repeat = 0; repeat < 2; repeat++) {
    const result = solveSparseDirectAlignedMany(a, rhs, order);
    result.forEach((r, k) => {
      assert(r.relativeResidual <= 1e-10);
      assert(Math.max(...r.x.map((v, i) => Math.abs(v - exact[k][i]))) < 1e-12);
    });
  }
  assert.deepEqual({ a, rhs }, before);
});

test('a difficult batched RHS retains the ordinary sparse accuracy fallback', () => {
  const { matrix: a, rhs } = JSON.parse(readFileSync(new URL('./fixtures/klu-refined-coupled.json', import.meta.url)));
  const order = Array.from({ length: a.n }, (_, i) => i);
  const exact = Float64Array.from({ length: a.n }, (_, i) => Math.sin(i / 7));
  const results = solveSparseDirectAlignedMany(a, [rhs, sparseProduct(a, exact)], order);
  assert(results.every(r => r.relativeResidual <= 1e-10));
  assert.equal(results[0].ordering, 'colamd');
});

test('batched solve rejects invalid and singular systems without poisoning later factors', () => {
  const a = sparseMatrix([[0, 1], [0, 1]]); a.values.set([1, 1, 2, 2]);
  assert.throws(() => solveSparseDirectAlignedMany(a, [[1, 0]], [0, 1]), /factorization/);
  assert.throws(() => solveSparseDirectAlignedMany(a, [[Infinity, 1]], [0, 1]), /right-hand sides/);
  assert.throws(() => solveSparseDirectAlignedMany(a, [[1, 1]], [0, 0]), /bijection/);
  assert.deepEqual(solveSparseDirectAlignedMany(a, [], [0, 1]), []);
  assert.deepEqual([...solveSparseDirectAlignedMany(a, [[0, 0]], [0, 1])[0].x], [0, 0]);
  a.values.set([0, 2, 3, 0]);
  assert.deepEqual([...solveSparseDirectAlignedMany(a, [[4, 9]], [1, 0])[0].x], [3, 2]);
});
