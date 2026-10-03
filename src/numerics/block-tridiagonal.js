// SPDX-License-Identifier: GPL-2.0-or-later
// Block Thomas with row-pivoted LU inside each Schur block. No inverse,
// positive-definiteness assumption, regularization or pivot across blocks.
import { factorLinear } from './linear.js';

export function solveBlockTridiagonal({ lower, diagonal, upper, rhs }) {
  const n = diagonal?.length, size = rhs?.[0]?.length;
  if (!Number.isInteger(n) || n < 1 || !Number.isInteger(size) || size < 1
    || ![lower, diagonal, upper, rhs].every(a => Array.isArray(a) && a.length === n)
    || [lower, diagonal, upper].some(a => a.some(b => !b || b.length !== size ** 2 || !Array.from(b).every(Number.isFinite)))
    || rhs.some(b => !b || b.length !== size || !Array.from(b).every(Number.isFinite)))
    throw new Error('Invalid block-tridiagonal system.');
  const transformedUpper = [], transformedRhs = [];
  for (let i = 0; i < n; i++) {
    const block = Float64Array.from(diagonal[i]), value = Float64Array.from(rhs[i]);
    if (i) for (let r = 0; r < size; r++) for (let k = 0; k < size; k++) {
      const a = lower[i][r * size + k];
      value[r] -= a * transformedRhs[i - 1][k];
      for (let c = 0; c < size; c++) block[r * size + c] -= a * transformedUpper[i - 1][k * size + c];
    }
    const solve = factorLinear(block, size);
    transformedRhs.push(solve(value));
    const next = new Float64Array(size ** 2);
    if (i + 1 < n) for (let c = 0; c < size; c++) {
      const column = solve(Array.from({ length: size }, (_, r) => upper[i][r * size + c]));
      for (let r = 0; r < size; r++) next[r * size + c] = column[r];
    }
    transformedUpper.push(next);
  }
  const solution = Array(n);
  for (let i = n - 1; i >= 0; i--) {
    const value = Float64Array.from(transformedRhs[i]);
    if (i + 1 < n) for (let r = 0; r < size; r++) for (let c = 0; c < size; c++)
      value[r] -= transformedUpper[i][r * size + c] * solution[i + 1][c];
    if (!value.every(Number.isFinite)) throw new Error('Nonfinite block-tridiagonal solution.');
    solution[i] = value;
  }
  return solution;
}
