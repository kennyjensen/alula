// SPDX-License-Identifier: GPL-2.0-or-later
// A dimensionless corner shape measure. For e=b-a, f=c-b and prescribed
// aspect r>0, E=(|e|²/r+r|f|²)/(2 cross(e,f)) >= 1 on positive corners.
// Equality means perpendicular edges with |e|/|f|=r. With one vertex free,
// the numerator is quadratic and the denominator affine. The derivatives
// below are for that vertex's two physical coordinates.
export function quadrilateralCornerShape(points, aspect, vertex) {
  if (!Array.isArray(points) || points.length !== 3 || !points.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))
    || !Number.isFinite(aspect) || aspect <= 0 || !Number.isInteger(vertex) || vertex < 0 || vertex > 2)
    throw new Error('Invalid quadrilateral corner shape inputs.');
  const [a, b, c] = points, e = { x: b.x - a.x, y: b.y - a.y }, f = { x: c.x - b.x, y: c.y - b.y };
  const determinant = e.x * f.y - e.y * f.x;
  if (!(determinant > 0) || !Number.isFinite(determinant)) throw new Error('Corner shape requires positive geometry.');
  const numerator = (e.x * e.x + e.y * e.y) / aspect + aspect * (f.x * f.x + f.y * f.y);
  const alpha = [-1, 1, 0][vertex], beta = [0, -1, 1][vertex];
  const dn = [2 * (alpha * e.x / aspect + beta * aspect * f.x), 2 * (alpha * e.y / aspect + beta * aspect * f.y)];
  const dd = [alpha * f.y - beta * e.y, -alpha * f.x + beta * e.x];
  const value = numerator / (2 * determinant), curvature = alpha * alpha / aspect + beta * beta * aspect;
  const gradient = dn.map((n, k) => (n / 2 - value * dd[k]) / determinant);
  const hessian = Array.from({ length: 4 }, (_, k) => {
    const i = Math.floor(k / 2), j = k % 2;
    return ((i === j ? curvature : 0) - (dn[i] * dd[j] + dd[i] * dn[j]) / (2 * determinant)
      + 2 * value * dd[i] * dd[j] / determinant) / determinant;
  });
  if (![value, ...gradient, ...hessian].every(Number.isFinite)) throw new Error('Unresolved quadrilateral corner shape.');
  return { value, gradient, hessian, determinant, determinantGradient: dd };
}
