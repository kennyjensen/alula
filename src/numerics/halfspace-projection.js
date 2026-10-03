// SPDX-License-Identifier: GPL-2.0-or-later
// Euclidean projection onto sparse halfspaces a_i.z >= b_i. By default the
// halfspaces contain zero. An explicit maximumNorm also permits shifted
// halfspaces, for a bounded minimum-norm second-order correction.
// Dual coordinate minimization (Hildreth/Dykstra); the squared-norm Hessian is
// the identity. No dense normal matrix is formed. A failed KKT check is a
// failed projection, never permission to use an infeasible point.
const norm = x => { let n = 0; for (const v of x) n = Math.hypot(n, v); return n; };
const dotRow = (row, x) => row.entries.reduce((s, [j, a]) => s + a * x[j], 0);
const primalToleranceRatio = (row, x, absoluteTolerance) => {
  let sum = 0, arithmeticScale = Math.abs(row.lower);
  for (const [j, a] of row.entries) { const term = a * x[j]; sum += term; arithmeticScale += Math.abs(term); }
  const violation = Math.max(0, row.lower - sum);
  // A zero-bound homogeneous face has no strictly positive margin to
  // preserve. Its usual absolute projection tolerance still applies.
  return violation / Math.max(Number.MIN_VALUE, row.lower === 0 ? absoluteTolerance : 16 * Number.EPSILON * arithmeticScale);
};

export function projectHalfspaces(point, constraints, { tolerance = 1e-11, maxSweeps = 500, maximumNorm } = {}) {
  if (!point.length || !point.every(Number.isFinite) || !(tolerance > 0) || !Number.isFinite(tolerance)
    || !Number.isInteger(maxSweeps) || maxSweeps < 1) throw new Error('Invalid halfspace projection controls.');
  const pointNorm = norm(point), radius = maximumNorm ?? pointNorm, rows = [];
  if (!Number.isFinite(radius) || radius < 0) throw new Error('Invalid projection norm bound.');
  let toleranceScale = pointNorm;
  for (const { gradient, lower } of constraints) {
    if (!(gradient instanceof Map) || !Number.isFinite(lower) || (maximumNorm === undefined && lower > 0))
      throw new Error('Unbounded projection halfspaces must contain the origin.');
    let length = 0;
    for (const [j, a] of gradient) {
      if (!Number.isInteger(j) || j < 0 || j >= point.length || !Number.isFinite(a)) throw new Error('Invalid projection gradient.');
      length = Math.hypot(length, a);
    }
    if (length === 0 && lower > 0) throw new Error('Inconsistent constant projection constraint.');
    // These halfspaces cannot constrain any point inside the permitted ball.
    // A projected point outside that ball is explicitly rejected below.
    if (length === 0 || -lower / length > radius) continue;
    rows.push({ entries: [...gradient].map(([j, a]) => [j, a / length]), lower: lower / length });
    toleranceScale = Math.max(toleranceScale, lower / length);
  }
  const x = Float64Array.from(point), multipliers = new Float64Array(rows.length);
  const absoluteTolerance = tolerance * Math.max(toleranceScale, Number.MIN_VALUE);
  let sweeps = 0, kktResidual = Infinity, maximumPrimalToleranceRatio = Infinity;
  do {
    for (let i = 0; i < rows.length; i++) {
      const row = rows[i], slack = dotRow(row, x) - row.lower;
      // Compute the increment directly. Subtracting two nearly equal
      // multipliers would discard corrections smaller than their ulp,
      // even when those corrections resolve a small primal margin.
      const change = Math.max(-multipliers[i], -slack);
      for (const [j, a] of row.entries) x[j] += change * a;
      multipliers[i] += change;
    }
    sweeps++;
    kktResidual = 0; maximumPrimalToleranceRatio = 0;
    for (let i = 0; i < rows.length; i++) {
      const slack = dotRow(rows[i], x) - rows[i].lower;
      // Projected dual-gradient residual includes feasibility and complementarity.
      kktResidual = Math.max(kktResidual, Math.abs(multipliers[i] - Math.max(0, multipliers[i] - slack)));
      maximumPrimalToleranceRatio = Math.max(maximumPrimalToleranceRatio, primalToleranceRatio(rows[i], x, absoluteTolerance));
    }
    // A global norm tolerance alone can exceed a small physical margin.
    // Resolve primal feasibility to each row's dot-product roundoff as
    // well as checking the global KKT residual. Do not relax the domain.
  } while ((kktResidual > absoluteTolerance || maximumPrimalToleranceRatio > 1) && sweeps < maxSweeps);
  const stationarity = x.map((v, j) => v - point[j]); let primalViolation = 0;
  rows.forEach((row, i) => {
    primalViolation = Math.max(primalViolation, row.lower - dotRow(row, x));
    for (const [j, a] of row.entries) stationarity[j] -= multipliers[i] * a;
  });
  const stationarityNorm = norm(stationarity), projectedNorm = norm(x);
  const converged = kktResidual <= absoluteTolerance && primalViolation <= absoluteTolerance && maximumPrimalToleranceRatio <= 1
    && stationarityNorm <= absoluteTolerance && projectedNorm <= radius + tolerance * radius;
  return { point: x, converged, sweeps, kktResidual, primalViolation, stationarityNorm,
    projectedNorm, maximumPrimalToleranceRatio, constraints: rows.length, active: multipliers.reduce((n, v) => n + (v > 0 ? 1 : 0), 0) };
}
