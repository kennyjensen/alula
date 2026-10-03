// SPDX-License-Identifier: GPL-2.0-or-later
// Bilinear quadrilateral weak Laplace matrix. Existing audits default to
// 3x3 integration; the Giles STIFF comparison uses its 2x2 Gauss rule.
import { requirePositiveSimpleQuad } from '../geometry/simple-quadrilateral.js';
const rules = {
  2: { points: [.5 - 1 / Math.sqrt(12), .5 + 1 / Math.sqrt(12)], weights: [.5, .5] },
  3: { points: [.5 - Math.sqrt(3 / 5) / 2, .5, .5 + Math.sqrt(3 / 5) / 2], weights: [5 / 18, 4 / 9, 5 / 18] },
};

export function quadLaplaceMatrix(vertices, { quadratureOrder = 3, quadratureDomain = 'convex' } = {}) {
  if (!['convex', 'sampled-positive'].includes(quadratureDomain)) throw new Error('Unknown quadrilateral quadrature domain.');
  if (![2, 3].includes(quadratureOrder)) throw new Error('Unsupported quadrilateral Laplace quadrature order.');
  const { points, weights } = rules[quadratureOrder];
  if (!Array.isArray(vertices) || vertices.length !== 4 || vertices.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))
    throw new Error('A harmonic audit requires four finite quadrilateral vertices.');
  // Translate before differentiating to retain invariance at remote origins.
  const p = vertices.map(q => ({ x: q.x - vertices[0].x, y: q.y - vertices[0].y }));
  // Sampled positivity reproduces the domain inspected by Giles's 2x2
  // STIFF formula for intermediate iterates. It is not a globally valid Q1
  // map certificate; final mesh acceptance must still require convexity.
  if (quadratureDomain === 'sampled-positive') requirePositiveSimpleQuad(vertices);
  else for (let k = 0; k < 4; k++) {
    const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4];
    if (!((b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x) > 0))
      throw new Error('Invalid quadrilateral in harmonic audit.');
  }
  const matrix = new Float64Array(16);
  for (let q = 0; q < quadratureOrder; q++) for (let r = 0; r < quadratureOrder; r++) {
    const s = points[q], t = points[r], ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
    let xs = 0, ys = 0, xt = 0, yt = 0;
    for (let k = 0; k < 4; k++) { xs += ds[k] * p[k].x; ys += ds[k] * p[k].y; xt += dt[k] * p[k].x; yt += dt[k] * p[k].y; }
    const det = xs * yt - ys * xt;
    if (!(det > 0)) throw new Error('Invalid quadrilateral in harmonic audit.');
    const gx = ds.map((v, k) => (yt * v - ys * dt[k]) / det);
    const gy = ds.map((v, k) => (-xt * v + xs * dt[k]) / det);
    const w = weights[q] * weights[r] * det;
    for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) matrix[4 * a + b] += w * (gx[a] * gx[b] + gy[a] * gy[b]);
  }
  return matrix;
}
