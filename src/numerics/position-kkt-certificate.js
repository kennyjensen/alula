// SPDX-License-Identifier: GPL-2.0-or-later
// Final equality-system accuracy qualification, independent of objective scale.
import { accurateRowResidual } from './position-kkt-projection.js';

export function certifyPositionKkt({ matrix, rhs, result, solve, xp, zp, scales, rows, active }) {
  const n = result.length, inverse = Array.from({ length: n }, () => new Float64Array(n)), columns = [];
  for (let j = 0; j < n; j++) {
    const unit = new Float64Array(n); unit[j] = 1;
    const column = solve(unit); columns.push(column); column.forEach((v, i) => { inverse[i][j] = v; });
  }
  const eps = Number.EPSILON, sumFactor = 1 + 32 * n * eps;
  let inverseDefectBound = 0;
  for (let i = 0; i < n; i++) {
    let rowSum = 0;
    for (let j = 0; j < n; j++) {
      let absoluteProducts = Number(i === j);
      for (const [k, a] of matrix[i]) absoluteProducts += Math.abs(a * inverse[k][j]);
      const residual = accurateRowResidual(matrix[i], columns[j], Number(i === j));
      const allowance = 64 * (matrix[i].size + 1) ** 2 * eps ** 2 * absoluteProducts
        + 8 * eps * Math.abs(residual) + 64 * (matrix[i].size + 1) * Number.MIN_VALUE;
      rowSum += Math.abs(residual) + allowance;
    }
    inverseDefectBound = Math.max(inverseDefectBound, sumFactor * rowSum);
  }
  if (!(inverseDefectBound < 1)) throw new Error('Position KKT inverse could not be certified.');
  const residual = Float64Array.from(rhs, (b, i) => accurateRowResidual(matrix[i], result, b));
  const residualRounding = Float64Array.from(rhs, (b, i) => {
    let products = Math.abs(b); for (const [j, a] of matrix[i]) products += Math.abs(a * result[j]);
    return 64 * (matrix[i].size + 1) ** 2 * eps ** 2 * products
      + 8 * eps * Math.abs(residual[i]) + 64 * (matrix[i].size + 1) * Number.MIN_VALUE;
  });
  const residualNorm = Math.max(...residual.map((v, i) => Math.abs(v) + residualRounding[i])) * sumFactor;
  const errorBounds = inverse.map(row => {
    const entries = Array.from(row, (a, i) => [i, a]);
    const correction = -accurateRowResidual(entries, residual);
    const absoluteProducts = row.reduce((s, a, i) => s + Math.abs(a * residual[i]), 0);
    const rowNorm = row.reduce((s, a) => s + Math.abs(a), 0) * sumFactor;
    const residualUncertainty = row.reduce((s, a, i) => s + Math.abs(a) * residualRounding[i], 0) * sumFactor;
    return (Math.abs(correction) + 64 * n ** 2 * eps ** 2 * absoluteProducts + 8 * eps * Math.abs(correction)
      + 64 * n * Number.MIN_VALUE + residualUncertainty
      + rowNorm * inverseDefectBound / (1 - inverseDefectBound) * residualNorm) * sumFactor;
  });
  const physicalPositionErrorBound = Math.max(...xp.map((k, i) => errorBounds[k] * scales[i]));
  const inactive = rows.map((row, k) => {
    if (active.includes(k)) return 0;
    const slack = accurateRowResidual(row.terms, xp.map(i => result[i]), row.rhs);
    const products = Math.abs(row.rhs) + row.terms.reduce((s, [i, a]) => s + Math.abs(a * result[xp[i]]), 0);
    const error = (row.terms.reduce((s, [i, a]) => s + Math.abs(a) * errorBounds[xp[i]], 0)
      + 64 * (row.terms.length + 1) ** 2 * eps ** 2 * products + 8 * eps * Math.abs(slack)
      + 64 * (row.terms.length + 1) * Number.MIN_VALUE) * sumFactor;
    return error - slack;
  });
  const maximumInactiveViolationBound = Math.max(0, ...inactive);
  const minimumActiveMultiplierBound = Math.min(Infinity, ...zp.map(k => result[k] - errorBounds[k] - 2 * eps * Math.abs(result[k])));
  if (![inverseDefectBound, physicalPositionErrorBound, maximumInactiveViolationBound].every(Number.isFinite)
    || errorBounds.some(v => !(v >= 0 && Number.isFinite(v)))) throw new Error('Nonfinite position accuracy certificate.');
  return { inverseDefectBound, physicalPositionErrorBound, maximumInactiveViolationBound,
    activeConstraintCount: active.length, minimumActiveMultiplierBound: active.length ? minimumActiveMultiplierBound : null,
    originalResidualNorm: residualNorm,
    method: 'verified approximate inverse of original normalized floating-point active KKT; station error, inactive slacks and multiplier signs' };
}
