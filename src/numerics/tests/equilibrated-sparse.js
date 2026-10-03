// SPDX-License-Identifier: GPL-2.0-or-later
// Change variable units, never the equations: A D y = b, x = D y.
// Both the scaled solve and the original matrix must meet the same tolerance.
import { solveSparseDirect } from '../klu.js';
import { compensatedSparseResidual } from '../compensated-sparse-residual.js';

function norm(values) {
  let scale = 0, sum = 0;
  for (const value of values) {
    if (!Number.isFinite(value)) return Infinity;
    const a = Math.abs(value);
    if (a > scale) { sum = 1 + sum * (scale / a) ** 2; scale = a; }
    else if (a) sum += (a / scale) ** 2;
  }
  return scale * Math.sqrt(sum);
}

export function solveEquilibratedSparse(a, b, { tolerance = 1e-10, maxRefinements = 12,
  ordering = 'auto', preferredOrdering = 'amd' } = {}) {
  if (!Number.isInteger(a.n) || a.n < 1 || a.rowPtr.length !== a.n + 1
    || a.colIndex.length !== a.values.length || a.rowPtr[a.n] !== a.values.length)
    throw new Error('Invalid sparse matrix dimensions.');
  const maxima = new Float64Array(a.n);
  for (let p = 0; p < a.values.length; p++) {
    const j = a.colIndex[p], v = a.values[p];
    if (!Number.isInteger(j) || j < 0 || j >= a.n || !Number.isFinite(v)) throw new Error('Invalid sparse matrix entry.');
    maxima[j] = Math.max(maxima[j], Math.abs(v));
  }
  if (maxima.some(v => !(v > 0))) throw new Error('Cannot equilibrate a zero matrix column.');
  const scaled = { ...a, values: a.values.map((v, p) => v / maxima[a.colIndex[p]]) };
  const linear = solveSparseDirect(scaled, b, { tolerance, maxRefinements, ordering, preferredOrdering });
  const x = linear.x.map((v, i) => v / maxima[i]), residual = compensatedSparseResidual(a, x, b);
  const rhsNorm = norm(b), relativeResidual = rhsNorm ? norm(residual) / rhsNorm : norm(residual);
  if (!x.every(Number.isFinite) || !(relativeResidual <= tolerance)) {
    const error = new Error(`Equilibrated sparse solve missed the original-system residual (${relativeResidual} > ${tolerance}).`);
    error.code = 'EQUILIBRATED_ORIGINAL_RESIDUAL_LIMIT'; error.relativeResidual = relativeResidual;
    error.scaledRelativeResidual = linear.relativeResidual; error.attempts = linear.attempts; throw error;
  }
  let minimum = Infinity, maximum = 0;
  for (const v of maxima) { minimum = Math.min(minimum, v); maximum = Math.max(maximum, v); }
  return { ...linear, x, relativeResidual, scaledRelativeResidual: linear.relativeResidual,
    columnMaximumRange: [minimum, maximum], equilibration: 'column-maximum' };
}
