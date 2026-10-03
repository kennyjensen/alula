// SPDX-License-Identifier: GPL-2.0-or-later
// Nine-node tensor-product Q2 scalar basis on a bilinear quadrilateral.
// Geometry corners are counterclockwise; basis node (a,b) has index 3*a+b
// and reference coordinates (a/2,b/2). Geometry remains bilinear, not curved.
import { gaussUnitRule } from '../gauss-unit.js';

export function quadraticQuadLaplaceMatrix(vertices, { quadratureOrder = 5 } = {}) {
  if (!Array.isArray(vertices) || vertices.length !== 4 || vertices.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))
    throw new Error('A quadratic harmonic reference requires four finite vertices.');
  const p = vertices.map(q => ({ x: q.x - vertices[0].x, y: q.y - vertices[0].y }));
  for (let k = 0; k < 4; k++) {
    const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4];
    if (!((b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x) > 0))
      throw new Error('Invalid quadrilateral in quadratic harmonic reference.');
  }
  const rule = gaussUnitRule(quadratureOrder);
  const basis = rule.nodes.map(s => ({ values: [2 * s * s - 3 * s + 1, 4 * s * (1 - s), 2 * s * s - s],
    derivatives: [4 * s - 3, 4 - 8 * s, 4 * s - 1] }));
  const matrix = new Float64Array(81), gx = new Float64Array(9), gy = new Float64Array(9);
  for (let q = 0; q < rule.nodes.length; q++) for (let r = 0; r < rule.nodes.length; r++) {
    const s = rule.nodes[q], t = rule.nodes[r], ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
    let xs = 0, ys = 0, xt = 0, yt = 0;
    for (let k = 0; k < 4; k++) { xs += ds[k] * p[k].x; ys += ds[k] * p[k].y; xt += dt[k] * p[k].x; yt += dt[k] * p[k].y; }
    const det = xs * yt - ys * xt;
    if (!(det > 0) || !Number.isFinite(det)) throw new Error('Invalid quadrilateral in quadratic harmonic reference.');
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
      const du = basis[q].derivatives[a] * basis[r].values[b], dv = basis[q].values[a] * basis[r].derivatives[b], k = 3 * a + b;
      gx[k] = (yt * du - ys * dv) / det; gy[k] = (-xt * du + xs * dv) / det;
    }
    const w = rule.weights[q] * rule.weights[r] * det;
    for (let a = 0; a < 9; a++) for (let b = 0; b < 9; b++) matrix[9 * a + b] += w * (gx[a] * gx[b] + gy[a] * gy[b]);
  }
  return matrix;
}
