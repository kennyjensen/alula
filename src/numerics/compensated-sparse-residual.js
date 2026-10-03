// SPDX-License-Identifier: GPL-2.0-or-later
// Evaluate the original b-Ax using error-free product and sum transforms.
// The low-order terms are accumulated separately and rounded into each row
// once. Fine-grid rows can otherwise lose the correction through cancellation
// even when sparse LU has a small componentwise backward error. This changes
// residual arithmetic for refinement and certification, not A, b or tolerance.
export function compensatedSparseResidual(a, x, b) {
  if (x.length !== a.n || b.length !== a.n) throw new Error('Residual dimensions differ.');
  const residual = new Float64Array(a.n), splitter = 134217729;
  for (let i = 0; i < a.n; i++) {
    let sum = b[i], correction = 0;
    for (let k = a.rowPtr[i]; k < a.rowPtr[i + 1]; k++) {
      const av = -a.values[k], xv = x[a.colIndex[k]], product = av * xv;
      const ac = splitter * av, xc = splitter * xv;
      const ah = ac - (ac - av), al = av - ah, xh = xc - (xc - xv), xl = xv - xh;
      const productError = ((ah * xh - product) + ah * xl + al * xh) + al * xl;
      const next = sum + product, bv = next - sum;
      const sumError = (sum - (next - bv)) + (product - bv);
      correction += productError + sumError; sum = next;
    }
    residual[i] = sum + correction;
    if (!Number.isFinite(residual[i])) throw new Error('Compensated residual overflow or nonfinite input.');
  }
  return residual;
}
