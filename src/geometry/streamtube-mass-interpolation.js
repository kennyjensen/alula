// SPDX-License-Identifier: GPL-2.0-or-later
// Refine cumulative mass in physical space. Local quadratic streamfunction
// gradients reproduce affine flows on skew crosslines and quadratic
// stagnation flow. A monotone rational-quadratic inverse locates child mass
// levels on each parent edge; all parent nodes and tube masses are retained.
import { inverseMonotoneRationalQuadratic } from '../numerics/rational-quadratic.js';

// Small column-pivoted Householder QR, avoiding squared condition numbers
// from normal equations. A rank-deficient quadratic fit may fall back to
// its affine subspace; an unresolved affine fit is an explicit failure.
function fit(rows, values, n) {
  const a = rows.map(r => r.slice(0, n)), b = values.slice(), m = a.length;
  const permutation = Array.from({ length: n }, (_, i) => i);
  const reference = Math.max(...Array.from({ length: n }, (_, j) => Math.hypot(...a.map(r => r[j]))));
  for (let k = 0; k < n; k++) {
    let pivot = k, norm = 0;
    for (let j = k; j < n; j++) {
      const candidate = Math.hypot(...a.slice(k).map(r => r[j]));
      if (candidate > norm) { norm = candidate; pivot = j; }
    }
    if (!(norm > 128 * Number.EPSILON * reference)) return null;
    for (const row of a) [row[k], row[pivot]] = [row[pivot], row[k]];
    [permutation[k], permutation[pivot]] = [permutation[pivot], permutation[k]];
    const v = a.slice(k).map(r => r[k]); v[0] += (v[0] >= 0 ? norm : -norm);
    const scale = Math.hypot(...v); for (let i = 0; i < v.length; i++) v[i] /= scale;
    for (let j = k; j < n; j++) {
      let dot = 0; for (let i = k; i < m; i++) dot += v[i - k] * a[i][j];
      for (let i = k; i < m; i++) a[i][j] -= 2 * v[i - k] * dot;
    }
    let dot = 0; for (let i = k; i < m; i++) dot += v[i - k] * b[i];
    for (let i = k; i < m; i++) b[i] -= 2 * v[i - k] * dot;
  }
  const solved = Array(n).fill(0), result = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let rhs = b[i]; for (let j = i + 1; j < n; j++) rhs -= a[i][j] * solved[j];
    solved[i] = rhs / a[i][i]; result[permutation[i]] = solved[i];
  }
  return result.every(Number.isFinite) ? result : null;
}

export function refineStreamtubeMassCoordinates(grid, masses, counts, { lowerStagnation = null, upperStagnation = null } = {}) {
  const nx = grid?.length - 1, nt = grid?.[0]?.length - 1;
  if (!Array.isArray(grid) || nx < 2 || nt < 2 || !grid.every(row => Array.isArray(row) && row.length === nt + 1
    && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))
    || !Array.isArray(masses) || masses.length !== nt || masses.some(m => !(m > 0) || !Number.isFinite(m))
    || !Array.isArray(counts) || counts.length !== nt || counts.some(n => !Number.isInteger(n) || n < 1 || n > 4)
    || [lowerStagnation, upperStagnation].some(i => i !== null && (!Number.isInteger(i) || i < 0 || i > nx)))
    throw new Error('Invalid streamtube mass interpolation data.');
  const diagnostics = { affineFallbacks: 0, limitedEndpointSlopes: 0, maximumFractionChange: 0 };
  const gradient = (i, j) => {
    if (j === 0 && i === lowerStagnation || j === nt && i === upperStagnation) return { x: 0, y: 0 };
    const startI = Math.min(nx - 2, Math.max(0, i - 1)), startJ = Math.min(nt - 2, Math.max(0, j - 1)), origin = grid[i][j];
    const data = [];
    for (let a = startI; a <= startI + 2; a++) for (let b = startJ; b <= startJ + 2; b++) if (a !== i || b !== j) {
      // Differences of nearby levels avoid cancellation against a large
      // arbitrary absolute passage streamfunction.
      let psi = 0; for (let k = Math.min(b, j); k < Math.max(b, j); k++) psi += masses[k];
      data.push({ x: grid[a][b].x - origin.x, y: grid[a][b].y - origin.y, psi: b < j ? -psi : psi });
    }
    const length = Math.max(...data.map(p => Math.hypot(p.x, p.y))), mass = Math.max(...data.map(p => Math.abs(p.psi)));
    if (!(length > 0 && mass > 0)) throw new Error('Degenerate mass interpolation stencil.');
    const rows = data.map(p => { const x = p.x / length, y = p.y / length; return [x, y, .5 * x * x, x * y, .5 * y * y]; });
    const rhs = data.map(p => p.psi / mass);
    let coefficients = fit(rows, rhs, 5);
    if (!coefficients) { coefficients = fit(rows, rhs, 2); diagnostics.affineFallbacks++; }
    if (!coefficients) throw new Error('Unresolved streamfunction gradient stencil.');
    return { x: coefficients[0] * mass / length, y: coefficients[1] * mass / length };
  };
  const nodes = grid.map((row, i) => {
    const result = [{ ...row[0] }], gradients = new Map();
    const slope = j => { if (!gradients.has(j)) gradients.set(j, gradient(i, j)); return gradients.get(j); };
    for (let j = 0; j < nt; j++) {
      const p = row[j], q = row[j + 1], dx = q.x - p.x, dy = q.y - p.y;
      if (!(Math.hypot(dx, dy) > 0)) throw new Error('Degenerate crossline interval.');
      let a, b;
      if (counts[j] > 1) [a, b] = [slope(j), slope(j + 1)].map(g => {
        const d = (g.x * dx + g.y * dy) / masses[j];
        if (d < 0) diagnostics.limitedEndpointSlopes++;
        return Math.max(0, d);
      });
      for (let k = 1; k < counts[j]; k++) {
        const fraction = k / counts[j], t = inverseMonotoneRationalQuadratic(fraction, a, b);
        diagnostics.maximumFractionChange = Math.max(diagnostics.maximumFractionChange, Math.abs(t - fraction));
        result.push({ x: p.x + t * dx, y: p.y + t * dy });
      }
      result.push({ ...q });
    }
    return result;
  });
  return { nodes, diagnostics };
}
