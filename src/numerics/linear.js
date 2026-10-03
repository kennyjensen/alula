// SPDX-License-Identifier: GPL-2.0-or-later
// Dense Float64 reference backend. Row-major buffers form the future WASM ABI.
export function normInf(vector) {
  let norm = 0;
  for (const value of vector) norm = Math.max(norm, Math.abs(value));
  return norm;
}

// Reusable row-equilibrated LU for an influence matrix with many right-hand
// sides. Pivot whole rows, including the already formed lower triangle.
export function factorLinear(matrix, n = Math.sqrt(matrix.length)) {
  if (!Number.isInteger(n) || n < 1 || matrix.length !== n * n) throw new Error('Invalid linear system dimensions.');
  const a = Float64Array.from(matrix); const scales = new Float64Array(n); const pivots = new Int32Array(n);
  if (!a.every(Number.isFinite)) throw new Error('Linear system contains nonfinite values.');
  for (let i = 0; i < n; i++) {
    scales[i] = normInf(a.subarray(i*n,(i+1)*n));
    if (scales[i] === 0) throw new Error('Singular linear system: zero row.');
    for (let j = 0; j < n; j++) a[i*n+j] /= scales[i];
  }
  for (let k = 0; k < n; k++) {
    let pivot = k;
    for (let i = k+1; i < n; i++) if (Math.abs(a[i*n+k]) > Math.abs(a[pivot*n+k])) pivot = i;
    if (Math.abs(a[pivot*n+k]) < 64*Number.EPSILON) throw new Error('Singular or numerically unresolved linear system.');
    pivots[k] = pivot;
    if (pivot !== k) for (let j = 0; j < n; j++) [a[k*n+j],a[pivot*n+j]] = [a[pivot*n+j],a[k*n+j]];
    for (let i = k+1; i < n; i++) {
      const f = a[i*n+k] /= a[k*n+k];
      if (f !== 0) for (let j = k+1; j < n; j++) a[i*n+j] -= f*a[k*n+j];
    }
  }
  return rhs => {
    if (rhs.length !== n || !Array.from(rhs).every(Number.isFinite)) throw new Error('Invalid linear right-hand side.');
    const x = Float64Array.from(rhs, (v,i) => v/scales[i]);
    for (let k = 0; k < n; k++) if (pivots[k] !== k) [x[k],x[pivots[k]]] = [x[pivots[k]],x[k]];
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) x[i] -= a[i*n+j]*x[j];
    for (let i = n-1; i >= 0; i--) {
      for (let j = i+1; j < n; j++) x[i] -= a[i*n+j]*x[j];
      x[i] /= a[i*n+i];
    }
    if (!x.every(Number.isFinite)) throw new Error('Nonfinite linear solution.');
    return x;
  };
}

export function solveLinear(matrix, rhs) {
  const n = rhs.length;
  if (!n || matrix.length !== n * n) throw new Error('Invalid linear system dimensions.');
  const a = Float64Array.from(matrix);
  const b = Float64Array.from(rhs);
  if (!a.every(Number.isFinite) || !b.every(Number.isFinite)) throw new Error('Linear system contains nonfinite values.');
  // Row equilibration prevents units from deciding which pivots are acceptable.
  for (let i = 0; i < n; i++) {
    const scale = normInf(a.subarray(i * n, (i + 1) * n));
    if (scale === 0) throw new Error('Singular linear system: zero row.');
    for (let j = 0; j < n; j++) a[i * n + j] /= scale;
    b[i] /= scale;
  }
  for (let k = 0; k < n; k++) {
    let pivot = k;
    for (let i = k + 1; i < n; i++) {
      if (Math.abs(a[i * n + k]) > Math.abs(a[pivot * n + k])) pivot = i;
    }
    if (Math.abs(a[pivot * n + k]) < 64 * Number.EPSILON) {
      throw new Error('Singular or numerically unresolved linear system.');
    }
    if (pivot !== k) {
      for (let j = k; j < n; j++) {
        [a[k * n + j], a[pivot * n + j]] = [a[pivot * n + j], a[k * n + j]];
      }
      [b[k], b[pivot]] = [b[pivot], b[k]];
    }
    for (let i = k + 1; i < n; i++) {
      const factor = a[i * n + k] / a[k * n + k];
      a[i * n + k] = 0;
      if (factor === 0) continue;
      for (let j = k + 1; j < n; j++) a[i * n + j] -= factor * a[k * n + j];
      b[i] -= factor * b[k];
    }
  }
  const x = new Float64Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let value = b[i];
    for (let j = i + 1; j < n; j++) value -= a[i * n + j] * x[j];
    x[i] = value / a[i * n + i];
  }
  if (!x.every(Number.isFinite)) throw new Error('Nonfinite linear solution.');
  return x;
}

export function linearResidual(matrix, x, rhs) {
  const n = x.length;
  const residual = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    residual[i] = -rhs[i];
    for (let j = 0; j < n; j++) residual[i] += matrix[i * n + j] * x[j];
  }
  return residual;
}
