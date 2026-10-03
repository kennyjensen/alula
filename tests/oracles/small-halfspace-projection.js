// SPDX-License-Identifier: GPL-2.0-or-later
// Independent bounded oracle: min ||x||²/2 subject to dense a_i.x >= b_i.
// Reduce the row span by twice-reorthogonalized QR, then enumerate independent
// active subsets in that small space. Never form A*A^T. At most 12 faces;
// this exhaustive reference is not a scalable production QP implementation.
const norm = a => { let n = 0; for (const v of a) n = Math.hypot(n, v); return n; };
const dot = (a, b) => {
  let sum = 0, correction = 0;
  for (let i = 0; i < a.length; i++) {
    const v = a[i] * b[i], next = sum + v;
    correction += Math.abs(sum) >= Math.abs(v) ? (sum - next) + v : (v - next) + sum; sum = next;
  }
  return sum + correction;
};
function rowQR(rows) {
  const basis = [], triangle = [];
  for (const row of rows) {
    const z = Float64Array.from(row), coefficients = new Array(basis.length).fill(0);
    for (let pass = 0; pass < 2; pass++) for (let j = 0; j < basis.length; j++) {
      const v = dot(z, basis[j]); coefficients[j] += v;
      for (let k = 0; k < z.length; k++) z[k] -= v * basis[j][k];
    }
    const length = norm(z);
    if (length > 128 * Number.EPSILON * norm(row)) { basis.push(z.map(v => v / length)); coefficients.push(length); }
    triangle.push(coefficients);
  }
  return { basis, triangle: triangle.map(row => Float64Array.from({ length: basis.length }, (_, i) => row[i] ?? 0)) };
}
const primal = (normal, lower, x) => {
  const slack = dot(normal, x) - lower;
  const scale = Math.abs(lower) + normal.reduce((s, v, i) => s + Math.abs(v * x[i]), 0);
  return { slack, tolerance: 64 * Number.EPSILON * scale };
};

export function projectSmallHalfspaces(faces, { maximumNorm = Infinity } = {}) {
  if (!Array.isArray(faces) || !faces.length || faces.length > 12 || !(maximumNorm >= 0))
    throw new Error('Small projection oracle requires 1–12 faces and a nonnegative norm bound.');
  const n = faces[0].normal.length;
  if (!n || faces.some(f => f.normal.length !== n || !f.normal.every(Number.isFinite) || !Number.isFinite(f.lower)))
    throw new Error('Invalid small projection faces.');
  const rows = faces.map(f => {
    const length = norm(f.normal);
    if (!length) throw new Error('Small projection requires nonzero face normals.');
    return { normal: Float64Array.from(f.normal, v => v / length), lower: f.lower / length, length };
  });
  const { basis, triangle: coordinates } = rowQR(rows.map(r => r.normal)), rank = basis.length;
  let subsets = 0;
  for (let mask = 0; mask < 2 ** rows.length; mask++) {
    const ids = rows.map((_, i) => i).filter(i => mask & 2 ** i), count = ids.length;
    if (count > rank) continue;
    subsets++;
    const small = rowQR(ids.map(i => coordinates[i]));
    if (small.basis.length !== count) continue;
    const coefficients = new Float64Array(count), multipliers = new Float64Array(count);
    for (let i = 0; i < count; i++) {
      let value = rows[ids[i]].lower;
      for (let j = 0; j < i; j++) value -= small.triangle[i][j] * coefficients[j];
      coefficients[i] = value / small.triangle[i][i];
    }
    for (let i = count - 1; i >= 0; i--) {
      let value = coefficients[i];
      for (let j = i + 1; j < count; j++) value -= small.triangle[j][i] * multipliers[j];
      multipliers[i] = value / small.triangle[i][i];
    }
    const dualTolerance = 64 * Number.EPSILON * norm(multipliers);
    if (!multipliers.every(v => Number.isFinite(v) && v >= -dualTolerance)) continue;
    const reduced = new Float64Array(rank);
    for (let i = 0; i < count; i++) for (let j = 0; j < rank; j++) reduced[j] += coefficients[i] * small.basis[i][j];
    if (rows.some((r, i) => { const p = primal(coordinates[i], r.lower, reduced); return p.slack < -p.tolerance; })) continue;
    const point = new Float64Array(n);
    for (let i = 0; i < rank; i++) for (let j = 0; j < n; j++) point[j] += reduced[i] * basis[i][j];
    const length = norm(point), stationarity = point.slice();
    const allMultipliers = new Float64Array(rows.length);
    ids.forEach((id, i) => {
      allMultipliers[id] = multipliers[i] / rows[id].length;
      for (let j = 0; j < n; j++) stationarity[j] -= multipliers[i] * rows[id].normal[j];
    });
    let primalViolation = 0, maximumPrimalToleranceRatio = 0, complementarity = 0;
    rows.forEach((r, i) => {
      const p = primal(r.normal, r.lower, point);
      primalViolation = Math.max(primalViolation, -p.slack);
      maximumPrimalToleranceRatio = Math.max(maximumPrimalToleranceRatio, Math.max(0, -p.slack) / Math.max(Number.MIN_VALUE, p.tolerance));
      complementarity += Math.abs(p.slack * allMultipliers[i] * r.length);
    });
    const stationarityNorm = norm(stationarity), tolerance = 1e-11 * Math.max(length, Number.MIN_VALUE);
    const converged = point.every(Number.isFinite) && maximumPrimalToleranceRatio <= 1
      && stationarityNorm <= tolerance && complementarity <= tolerance * Math.max(length, Number.MIN_VALUE)
      && length <= maximumNorm * (1 + 1e-11);
    if (converged) return { point, converged, rank, subsets, active: count, activeIndices: ids,
      multipliers: allMultipliers, projectedNorm: length, stationarityNorm, primalViolation, maximumPrimalToleranceRatio, complementarity };
  }
  return { point: null, converged: false, rank, subsets, reason: 'No independently certified active subset within the oracle limits.' };
}
