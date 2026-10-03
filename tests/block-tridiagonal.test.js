import test from 'node:test';
import assert from 'node:assert/strict';
import { solveBlockTridiagonal } from '../src/numerics/block-tridiagonal.js';
import { solveLinear } from '../src/numerics/linear.js';

test('general block Thomas matches known solutions and an independent dense solve with internal row pivots', () => {
  for (const n of [1, 5]) for (const size of [1, 2, 4, 6]) {
    const N = n * size, matrix = new Float64Array(N * N), expected = Array.from({ length: N }, (_, k) => Math.sin(k + .2));
    const diagonal = [], lower = [], upper = [];
    for (let i = 0; i < n; i++) {
      const d = new Float64Array(size ** 2), l = new Float64Array(size ** 2), u = new Float64Array(size ** 2);
      for (let r = 0; r < size; r++) for (let c = 0; c < size; c++) {
        // Permute rows inside every diagonal block, requiring LU pivots.
        const row = size > 1 ? (r + 1) % size : r, scale = 10 ** (2 * r - size);
        d[r * size + c] = scale * (row === c ? -6 : .1 * Math.sin(row + c + 1));
        l[r * size + c] = scale * .06 * Math.cos(row + 2 * c + i);
        u[r * size + c] = scale * .08 * Math.sin(2 * row + c + i);
        matrix[(i * size + r) * N + i * size + c] = d[r * size + c];
        if (i) matrix[(i * size + r) * N + (i - 1) * size + c] = l[r * size + c];
        if (i + 1 < n) matrix[(i * size + r) * N + (i + 1) * size + c] = u[r * size + c];
      }
      diagonal.push(d); lower.push(l); upper.push(u);
    }
    const b = Array.from({ length: N }, (_, r) => expected.reduce((s, x, c) => s + matrix[r * N + c] * x, 0));
    const input = { diagonal, lower, upper, rhs: Array.from({ length: n }, (_, i) => b.slice(i * size, (i + 1) * size)) };
    const before = structuredClone(input), actual = solveBlockTridiagonal(input).flatMap(row => Array.from(row));
    const dense = solveLinear(matrix, b);
    for (let i = 0; i < N; i++) { assert.ok(Math.abs(actual[i] - expected[i]) < 2e-12); assert.ok(Math.abs(actual[i] - dense[i]) < 2e-12); }
    assert.deepEqual(input, before);
  }
});

test('block Thomas rejects nonfinite, inconsistent and unresolved Schur blocks', () => {
  const input = { lower: [[0]], diagonal: [[1]], upper: [[0]], rhs: [[2]] };
  assert.throws(() => solveBlockTridiagonal({ ...input, diagonal: [[NaN]] }), /Invalid/);
  assert.throws(() => solveBlockTridiagonal({ ...input, rhs: [[1, 2]] }), /Invalid/);
  assert.throws(() => solveBlockTridiagonal({ ...input, diagonal: [[0]] }), /Singular/);
  assert.throws(() => solveBlockTridiagonal({ lower: [[0], [1]], diagonal: [[1], [1]], upper: [[1], [0]], rhs: [[1], [1]] }), /Singular/);
});
