// SPDX-License-Identifier: GPL-2.0-or-later
// Fixed pairs of transverse and streamwise Gauss--Seidel line corrections,
// as in the DSOLV procedure. Each line uses an exact tridiagonal factorization.
// Assemble from the original sparse matrix, avoiding separately indexed
// copies of the nine-point coefficients. Report the residual of that matrix.
import { sparseIndex, sparseProduct } from './sparse.js';

const norm = values => {
  let largest = 0;
  for (const v of values) largest = Math.max(largest, Math.abs(v));
  if (!largest) return 0;
  let sum = 0; for (const v of values) sum += (v / largest) ** 2;
  return largest * Math.sqrt(sum);
};

export function solveAlternatingScalarLines(matrix, rhs, { nx, nt, pairs = 5, tolerance = 1e-10 } = {}) {
  const n = (nx - 1) * (nt - 1);
  if (!Number.isInteger(nx) || nx < 2 || !Number.isInteger(nt) || nt < 2 || matrix?.n !== n
    || matrix.rowPtr?.length !== n + 1 || matrix.rowPtr[0] !== 0 || matrix.rowPtr[n] !== matrix.values?.length
    || matrix.colIndex?.length !== matrix.values?.length || rhs?.length !== n)
    throw new Error('Invalid alternating scalar line dimensions.');
  if (!Number.isInteger(pairs) || pairs < 0 || !Number.isFinite(tolerance) || tolerance <= 0
    || !Array.from(rhs).every(Number.isFinite)) throw new Error('Invalid alternating scalar line controls.');
  for (let row = 0; row < n; row++) {
    const start = matrix.rowPtr[row], end = matrix.rowPtr[row + 1]; let previous = -1;
    if (!Number.isInteger(start) || !Number.isInteger(end) || end < start) throw new Error('Invalid scalar line row pointers.');
    for (let p = start; p < end; p++) {
      const col = matrix.colIndex[p], v = matrix.values[p];
      if (!Number.isInteger(col) || col <= previous || col >= n || !Number.isFinite(v)) throw new Error('Invalid scalar line matrix entry.');
      if (v && (Math.abs(Math.floor(row / (nt - 1)) - Math.floor(col / (nt - 1))) > 1
        || Math.abs(row % (nt - 1) - col % (nt - 1)) > 1)) throw new Error('Alternating lines require a local nine-point stencil.');
      previous = col;
    }
  }
  const id = (i, j) => (i - 1) * (nt - 1) + j - 1;
  const entry = (row, col) => { const p = sparseIndex(matrix, row, col); return p < 0 ? 0 : matrix.values[p]; };
  const factor = ids => {
    const diagonal = new Float64Array(ids.length), upper = new Float64Array(ids.length), multipliers = new Float64Array(ids.length);
    ids.forEach((row, k) => {
      const original = entry(row, row), lower = k ? entry(row, ids[k - 1]) : 0;
      upper[k] = k + 1 < ids.length ? entry(row, ids[k + 1]) : 0;
      multipliers[k] = k ? lower / diagonal[k - 1] : 0;
      const correction = k ? multipliers[k] * upper[k - 1] : 0;
      diagonal[k] = original - correction;
      const scale = Math.max(Math.abs(original), Math.abs(lower), Math.abs(upper[k]), Math.abs(correction));
      if (!Number.isFinite(diagonal[k]) || !(diagonal[k] > 32 * Number.EPSILON * scale))
        throw new Error(`Unresolved positive alternating line pivot at row ${row}.`);
    });
    return { ids, diagonal, upper, multipliers, delta: new Float64Array(ids.length) };
  };
  const lines = [
    ...Array.from({ length: nx - 1 }, (_, i) => factor(Array.from({ length: nt - 1 }, (_, j) => id(i + 1, j + 1)))),
    ...Array.from({ length: nt - 1 }, (_, j) => factor(Array.from({ length: nx - 1 }, (_, i) => id(i + 1, j + 1)))),
  ];
  const x = new Float64Array(n), rhsNorm = norm(rhs), history = [];
  if (!Number.isFinite(rhsNorm)) throw new Error('Nonfinite alternating line right-hand-side norm.');
  const measure = () => {
    const residual = sparseProduct(matrix, x);
    for (let row = 0; row < n; row++) residual[row] -= rhs[row];
    const residualNorm = norm(residual), relativeResidual = rhsNorm ? residualNorm / rhsNorm : residualNorm;
    if (!Number.isFinite(relativeResidual)) throw new Error('Nonfinite alternating line residual.');
    return { residualNorm, relativeResidual };
  };
  history.push({ pair: 0, ...measure() });
  for (let pair = 1; pair <= pairs; pair++) {
    for (const { ids, diagonal, upper, multipliers, delta } of lines) {
      ids.forEach((row, k) => {
        let r = rhs[row];
        for (let p = matrix.rowPtr[row]; p < matrix.rowPtr[row + 1]; p++) r -= matrix.values[p] * x[matrix.colIndex[p]];
        delta[k] = r;
      });
      for (let k = 1; k < ids.length; k++) delta[k] -= multipliers[k] * delta[k - 1];
      for (let k = ids.length - 1; k >= 0; k--)
        delta[k] = (delta[k] - (k + 1 < ids.length ? upper[k] * delta[k + 1] : 0)) / diagonal[k];
      ids.forEach((row, k) => { x[row] += delta[k]; });
    }
    if (!x.every(Number.isFinite)) throw new Error('Nonfinite alternating line update.');
    history.push({ pair, ...measure() });
  }
  return { x, pairs, ...history.at(-1), converged: history.at(-1).relativeResidual <= tolerance, history };
}
