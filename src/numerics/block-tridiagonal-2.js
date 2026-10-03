// SPDX-License-Identifier: GPL-2.0-or-later
// General nonsymmetric 2x2 block Thomas solve, O(number of blocks).
// Block pivots must be resolved; no SPD assumption or silent regularization.
const product = (a, b) => [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
  a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]];
const apply = (a, b) => [a[0] * b[0] + a[1] * b[1], a[2] * b[0] + a[3] * b[1]];
const inverse = a => {
  const scale = Math.max(...a.map(Math.abs));
  if (!(scale > 0) || !Number.isFinite(scale)) throw new Error('Singular 2x2 block line pivot.');
  const b = a.map(v => v / scale), determinant = b[0] * b[3] - b[1] * b[2];
  if (!(Math.abs(determinant) > 64 * Number.EPSILON)) throw new Error('Unresolved 2x2 block line pivot.');
  return [b[3], -b[1], -b[2], b[0]].map(v => v / determinant / scale);
};
export function solveBlockTridiagonal2({ lower, diagonal, upper, rhs }) {
  const n = diagonal?.length;
  if (!Number.isInteger(n) || n < 1 || ![lower, diagonal, upper, rhs].every(a => Array.isArray(a) && a.length === n)
    || [lower, diagonal, upper].some(a => a.some(b => !Array.isArray(b) || b.length !== 4 || !b.every(Number.isFinite)))
    || rhs.some(b => !Array.isArray(b) || b.length !== 2 || !b.every(Number.isFinite))) throw new Error('Invalid 2x2 block line.');
  const pivots = [], residual = rhs.map(row => row.slice());
  for (let k = 0; k < n; k++) {
    let block = diagonal[k].slice();
    if (k) {
      const factor = product(lower[k], pivots[k - 1]), correction = product(factor, upper[k - 1]), value = apply(factor, residual[k - 1]);
      block = block.map((v, c) => v - correction[c]);
      residual[k] = residual[k].map((v, c) => v - value[c]);
    }
    pivots.push(inverse(block));
  }
  const solution = Array(n);
  for (let k = n - 1; k >= 0; k--) {
    const correction = k + 1 < n ? apply(upper[k], solution[k + 1]) : [0, 0];
    solution[k] = apply(pivots[k], residual[k].map((v, c) => v - correction[c]));
    if (!solution[k].every(Number.isFinite)) throw new Error('Nonfinite 2x2 block line solution.');
  }
  return solution;
}
