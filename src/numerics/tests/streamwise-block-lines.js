// SPDX-License-Identifier: GPL-2.0-or-later
// Two-coordinate block SLOR for a quadrilateral nine-point stencil. Each
// streamwise line uses block Thomas elimination, linear work in its length.
import { sparseIndex, sparseProduct } from '../sparse.js';

const product = (a, b) => [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
  a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]];
const apply = (a, x, y) => [a[0] * x + a[1] * y, a[2] * x + a[3] * y];
const inverse = a => {
  const scale = Math.max(...a.map(Math.abs)), b = a.map(v => v / scale), det = b[0] * b[3] - b[1] * b[2];
  if (!(scale > 0 && b[0] > 0 && b[3] > 0 && det > 32 * Number.EPSILON))
    throw new Error('Unresolved positive-definite grid line block.');
  return [b[3], -b[1], -b[2], b[0]].map(v => v / det / scale);
};

export function factorStreamwiseBlockLines(matrix, { nx, nt }) {
  if (!Number.isInteger(nx) || !Number.isInteger(nt) || nx < 2 || nt < 2 || matrix.n !== 2 * (nx - 1) * (nt - 1))
    throw new Error('Invalid two-coordinate line dimensions.');
  for (let r = 0; r < matrix.n; r++) for (let p = matrix.rowPtr[r]; p < matrix.rowPtr[r + 1]; p++) {
    const a = Math.floor(r / 2), b = Math.floor(matrix.colIndex[p] / 2);
    if (!Number.isFinite(matrix.values[p]) || (a % (nt - 1) === b % (nt - 1)
      && Math.abs(Math.floor(a / (nt - 1)) - Math.floor(b / (nt - 1))) > 1 && matrix.values[p] !== 0))
      throw new Error('Grid matrix must have finite, block-tridiagonal streamwise lines.');
  }
  const id = (i, j) => 2 * ((i - 1) * (nt - 1) + j - 1);
  const block = (a, b) => [0, 1].flatMap(r => [0, 1].map(c => {
    const k = sparseIndex(matrix, a + r, b + c); return k < 0 ? 0 : matrix.values[k];
  }));
  const lines = Array.from({ length: nt - 1 }, (_, j0) => {
    const j = j0 + 1, inverses = [], factors = [], upper = [];
    for (let i = 1; i < nx; i++) {
      const k = i - 1, diagonal = block(id(i, j), id(i, j));
      upper[k] = i + 1 < nx ? block(id(i, j), id(i + 1, j)) : [0, 0, 0, 0];
      if (k) {
        factors[k] = product(block(id(i, j), id(i - 1, j)), inverses[k - 1]);
        const correction = product(factors[k], upper[k - 1]);
        for (let q = 0; q < 4; q++) diagonal[q] -= correction[q];
      }
      inverses[k] = inverse(diagonal);
    }
    return { j, inverses, factors, upper };
  });
  const sweep = (x, rhs, omega = 1) => {
    if (x.length !== matrix.n || rhs.length !== matrix.n || !(omega > 0 && omega < 2) || !Number.isFinite(omega))
      throw new Error('Invalid block SLOR state or relaxation.');
    const next = Float64Array.from(x);
    for (const { j, inverses, factors, upper } of lines) {
      const delta = new Float64Array(2 * (nx - 1));
      for (let i = 1; i < nx; i++) for (let c = 0; c < 2; c++) {
        const row = id(i, j) + c; let value = rhs[row];
        for (let p = matrix.rowPtr[row]; p < matrix.rowPtr[row + 1]; p++) value -= matrix.values[p] * next[matrix.colIndex[p]];
        delta[2 * (i - 1) + c] = value;
      }
      for (let k = 1; k < nx - 1; k++) {
        const q = apply(factors[k], delta[2 * k - 2], delta[2 * k - 1]);
        delta[2 * k] -= q[0]; delta[2 * k + 1] -= q[1];
      }
      for (let k = nx - 2; k >= 0; k--) {
        const q = k + 1 < nx - 1 ? apply(upper[k], delta[2 * k + 2], delta[2 * k + 3]) : [0, 0];
        const d = apply(inverses[k], delta[2 * k] - q[0], delta[2 * k + 1] - q[1]);
        delta[2 * k] = d[0]; delta[2 * k + 1] = d[1];
        next[id(k + 1, j)] += omega * d[0]; next[id(k + 1, j) + 1] += omega * d[1];
      }
    }
    if (!next.every(Number.isFinite)) throw new Error('Nonfinite block SLOR update.');
    return next;
  };
  return { sweep };
}

// Standard overrelaxation is applicable to the SPD Gauss-Newton matrix;
// this omega does not overrelax the nonlinear physical-grid update.
export function solveStreamwiseBlockLines(matrix, rhs, dimensions, { maxSweeps = 600, tolerance = 1e-10, omega = 1.5 } = {}) {
  if (!Number.isInteger(maxSweeps) || maxSweeps < 0 || !(tolerance > 0) || !Number.isFinite(tolerance)
    || !(omega > 0 && omega < 2) || !Number.isFinite(omega)
    || rhs.length !== matrix.n || !rhs.every(Number.isFinite)) throw new Error('Invalid block line solve controls.');
  const lines = factorStreamwiseBlockLines(matrix, dimensions), norm = a => Math.sqrt(a.reduce((s, v) => s + v * v, 0));
  const rhsNorm = norm(rhs); let x = new Float64Array(matrix.n), relativeResidual = rhsNorm ? 1 : 0, sweeps = 0;
  for (; sweeps < maxSweeps && relativeResidual > tolerance; sweeps++) {
    x = lines.sweep(x, rhs, omega);
    relativeResidual = norm(sparseProduct(matrix, x).map((v, i) => v - rhs[i])) / rhsNorm;
    if (!Number.isFinite(relativeResidual)) throw new Error('Nonfinite block line solve residual.');
  }
  return { x, relativeResidual, sweeps, converged: relativeResidual <= tolerance };
}
