import test from 'node:test';
import assert from 'node:assert/strict';
import { compensatedSparseResidual } from '../src/numerics/compensated-sparse-residual.js';

test('compensated residual retains an exactly representable term lost to cancellation', () => {
  const a = { n: 3, rowPtr: [0, 3, 4, 5], colIndex: [0, 1, 2, 1, 2], values: [1e16, 1, -1e16, 1, 1] };
  assert.deepEqual(Array.from(compensatedSparseResidual(a, [1, 1, 1], [0, 1, 1])), [-1, 0, 0]);
});

test('compensated residual recovers the exact product error of binary factors', () => {
  // (1+2^-27)(1-2^-27)-1 = -2^-54, although the first product rounds to 1.
  const a = { n: 2, rowPtr: [0, 2, 3], colIndex: [0, 1, 1], values: [1 + 2 ** -27, -1, 1] };
  assert.deepEqual(Array.from(compensatedSparseResidual(a, [1 - 2 ** -27, 1], [0, 1])), [2 ** -54, 0]);
});
