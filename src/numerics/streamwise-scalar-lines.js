// SPDX-License-Identifier: GPL-2.0-or-later
// Scalar SLOR: exact tridiagonal Thomas solves along fixed eta lines,
// Gauss--Seidel between lines, and a certificate from the original matrix.
// Line pivots must be positive and resolved. No shifts or matrix repairs.
import { sparseIndex, sparseProduct } from './sparse.js';

const norm = values => {
  let scale = 0;
  for (const value of values) scale = Math.max(scale, Math.abs(value));
  if (scale === 0) return 0;
  let sum = 0;
  for (const value of values) sum += (value / scale) ** 2;
  return scale * Math.sqrt(sum);
};

export function solveStreamwiseScalarLines(matrix, rhs, {
  nx, nt, maxSweeps = 600, tolerance = 1e-10, omega = 1.3,
} = {}) {
  const n = (nx - 1) * (nt - 1);
  if (!Number.isInteger(nx) || nx < 2 || !Number.isInteger(nt) || nt < 2
    || matrix?.n !== n || matrix.rowPtr?.length !== n + 1
    || !matrix.colIndex || !matrix.values
    || matrix.colIndex?.length !== matrix.values?.length
    || matrix.rowPtr[0] !== 0 || matrix.rowPtr[n] !== matrix.values.length)
    throw new Error('Invalid scalar streamwise line dimensions or sparse storage.');
  if (!Number.isInteger(maxSweeps) || maxSweeps < 0 || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isFinite(omega) || omega <= 0 || omega >= 2
    || rhs?.length !== n || !Array.from(rhs).every(Number.isFinite))
    throw new Error('Invalid scalar line solve controls or right-hand side.');
  for (let row = 0; row < n; row++) {
    const start = matrix.rowPtr[row], end = matrix.rowPtr[row + 1];
    if (!Number.isInteger(start) || !Number.isInteger(end) || end < start)
      throw new Error('Invalid scalar line sparse row pointers.');
    let previous = -1;
    for (let p = start; p < end; p++) {
      const col = matrix.colIndex[p], value = matrix.values[p];
      if (!Number.isInteger(col) || col <= previous || col >= n || !Number.isFinite(value))
        throw new Error('Scalar line matrix must have ordered columns and finite entries.');
      if (value !== 0 && row % (nt - 1) === col % (nt - 1)
        && Math.abs(Math.floor(row / (nt - 1)) - Math.floor(col / (nt - 1))) > 1)
        throw new Error('Scalar streamwise lines must be tridiagonal.');
      previous = col;
    }
  }
  const id = (i, j) => (i - 1) * (nt - 1) + j - 1;
  const entry = (row, col) => {
    const p = sparseIndex(matrix, row, col);
    return p < 0 ? 0 : matrix.values[p];
  };
  const lines = Array.from({ length: nt - 1 }, (_, j0) => {
    const j = j0 + 1, diagonal = new Float64Array(nx - 1), upper = new Float64Array(nx - 1);
    const factors = new Float64Array(nx - 1);
    for (let i = 1; i < nx; i++) {
      const k = i - 1, row = id(i, j), original = entry(row, row);
      const lower = k ? entry(row, id(i - 1, j)) : 0;
      upper[k] = i + 1 < nx ? entry(row, id(i + 1, j)) : 0;
      factors[k] = k ? lower / diagonal[k - 1] : 0;
      const correction = k ? factors[k] * upper[k - 1] : 0;
      diagonal[k] = original - correction;
      const scale = Math.max(Math.abs(original), Math.abs(lower), Math.abs(upper[k]), Math.abs(correction));
      if (!Number.isFinite(factors[k]) || !Number.isFinite(diagonal[k])
        || !(diagonal[k] > 32 * Number.EPSILON * scale))
        throw new Error(`Unresolved positive scalar line pivot at station ${i}, streamline ${j}.`);
    }
    return { j, diagonal, upper, factors };
  });

  const rhsNorm = norm(rhs), x = new Float64Array(n), delta = new Float64Array(nx - 1);
  if (!Number.isFinite(rhsNorm)) throw new Error('Nonfinite scalar line right-hand-side norm.');
  let relativeResidual = rhsNorm === 0 ? 0 : 1, sweeps = 0;
  for (; sweeps < maxSweeps && relativeResidual > tolerance; sweeps++) {
    for (const { j, diagonal, upper, factors } of lines) {
      // All entries in this line use the current full residual before any
      // correction in the line; earlier eta lines have already been updated.
      for (let i = 1; i < nx; i++) {
        const row = id(i, j); let value = rhs[row];
        for (let p = matrix.rowPtr[row]; p < matrix.rowPtr[row + 1]; p++)
          value -= matrix.values[p] * x[matrix.colIndex[p]];
        delta[i - 1] = value;
      }
      for (let k = 1; k < nx - 1; k++) delta[k] -= factors[k] * delta[k - 1];
      for (let k = nx - 2; k >= 0; k--) {
        delta[k] = (delta[k] - (k + 1 < nx - 1 ? upper[k] * delta[k + 1] : 0)) / diagonal[k];
        x[id(k + 1, j)] += omega * delta[k];
      }
    }
    if (!x.every(Number.isFinite)) throw new Error('Nonfinite scalar SLOR update.');
    const residual = sparseProduct(matrix, x);
    for (let row = 0; row < n; row++) residual[row] -= rhs[row];
    relativeResidual = norm(residual) / rhsNorm;
    if (!Number.isFinite(relativeResidual)) throw new Error('Nonfinite scalar line solve residual.');
  }
  return { x, converged: relativeResidual <= tolerance, sweeps, relativeResidual };
}
