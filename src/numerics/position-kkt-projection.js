// SPDX-License-Identifier: GPL-2.0-or-later
import { certifyPositionKkt } from './position-kkt-certificate.js';
// Internal sparse banded equality-constrained position projection.
export function accurateRowResidual(terms, x, rhs = 0) {
  let sum = rhs, correction = 0;
  for (const [i, a] of terms) {
    const v = -a, y = x[i], p = v * y, vc = 134217729 * v, yc = 134217729 * y;
    const vh = vc - (vc - v), vl = v - vh, yh = yc - (yc - y), yl = y - yh;
    const pe = ((vh * yh - p) + vh * yl + vl * yh) + vl * yl;
    const next = sum + p, t = next - sum;
    correction += pe + (sum - (next - t)) + (p - t); sum = next;
  }
  return sum + correction;
}
export function projectPositionKkt({ active, rows, scales, certify = false }) {
  const n = scales.length, byCenter = Array.from({ length: n }, () => []);
  if (!n || scales.some(v => !(v > 0 && Number.isFinite(v))) || new Set(active).size !== active.length)
    throw new Error('Invalid position KKT dimensions or weights.');
  active.forEach((r, k) => {
    if (!Number.isInteger(r) || r < 0 || r >= rows.length) throw new Error('Invalid position KKT row index.');
    const row = rows[r];
    if (!Number.isFinite(row.rhs) || !row.terms.length || row.terms.length > 3
      || new Set(row.terms.map(t => t[0])).size !== row.terms.length
      || row.terms.some(([i, a]) => !Number.isInteger(i) || i < 0 || i >= n || !Number.isFinite(a))
      || Math.max(...row.terms.map(t => t[0])) - Math.min(...row.terms.map(t => t[0])) > 2)
      throw new Error('Invalid position KKT local row.');
    const center = Math.round(rows[r].terms.reduce((s, t) => s + t[0], 0) / rows[r].terms.length);
    byCenter[center].push(k);
  });
  const xp = [], zp = []; let size = 0;
  for (let i = 0; i < n; i++) { xp[i] = size++; for (const k of byCenter[i]) zp[k] = size++; }
  const matrix = Array.from({ length: size }, () => new Map()), rhs = new Float64Array(size);
  let bandwidth = 0;
  const add = (i, j, v) => { matrix[i].set(j, v); bandwidth = Math.max(bandwidth, Math.abs(i - j)); };
  for (const i of xp) add(i, i, 1);
  active.forEach((r, k) => {
    rhs[zp[k]] = rows[r].rhs;
    for (const [i, a] of rows[r].terms) { add(xp[i], zp[k], a); add(zp[k], xp[i], a); }
  });
  const a = matrix.map(row => new Map(row)), operations = [];
  for (let k = 0; k < size; k++) {
      let pivot = k, magnitude = Math.abs(a[k].get(k) ?? 0);
      for (let i = k + 1; i < Math.min(size, k + bandwidth + 1); i++) {
        const v = Math.abs(a[i].get(k) ?? 0); if (v > magnitude) { pivot = i; magnitude = v; }
      }
      if (!(magnitude > 0 && Number.isFinite(magnitude))) throw new Error('Unresolved position KKT pivot.');
      if (pivot !== k) [a[k], a[pivot]] = [a[pivot], a[k]];
      const operation = { pivot, updates: [] }; operations.push(operation);
      for (let i = k + 1; i < Math.min(size, k + bandwidth + 1); i++) {
        const v = a[i].get(k); if (!v) continue;
        const q = v / a[k].get(k); a[i].delete(k);
        for (const [j, u] of a[k]) if (j > k) {
          const updated = (a[i].get(j) ?? 0) - q * u;
          if (updated === 0) a[i].delete(j); else a[i].set(j, updated);
        }
        operation.updates.push([i, q]);
      }
  }
  const solve = input => {
    const b = Float64Array.from(input);
    operations.forEach(({pivot, updates}, k) => {
      if (pivot !== k) [b[k], b[pivot]] = [b[pivot], b[k]];
      for (const [i, q] of updates) b[i] -= q * b[k];
    });
    const result = new Float64Array(size);
    for (let k = size - 1; k >= 0; k--) {
      let v = b[k]; for (const [j, aij] of a[k]) if (j > k) v -= aij * result[j];
      result[k] = v / a[k].get(k);
    }
    return result;
  };
  const result = solve(rhs);
  for (let pass = 0; pass < 2; pass++) {
    const residual = Float64Array.from(rhs, (b, i) => accurateRowResidual(matrix[i], result, b));
    const update = solve(residual); result.forEach((_, i) => { result[i] += update[i]; });
    if (![...residual, ...update, ...result].every(Number.isFinite)) throw new Error('Nonfinite position KKT refinement.');
  }
  const target = xp.map(i => result[i]), multipliers = zp.map(i => result[i]);
  const forwardCertificate = certify ? certifyPositionKkt({ matrix, rhs, result, solve, xp, zp, scales, rows, active }) : undefined;
  return { target, multipliers, ...(forwardCertificate ? { forwardCertificate } : {}) };
}
