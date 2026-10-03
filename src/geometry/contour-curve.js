// SPDX-License-Identifier: GPL-2.0-or-later
// Parametric C2 natural cubic through an open TE-to-TE contour. The knot
// parameter is cumulative input chord length; it approaches arc length under
// geometry refinement. The repeated TE has independent one-sided tangents.
// Linear coordinate interpolation preserves rigid-transform invariance.
export function createContourCurve(points, { allowOpenEndpoints = false } = {}) {
  if (typeof allowOpenEndpoints !== 'boolean') throw new Error('Invalid open contour endpoint option.');
  if (!Array.isArray(points) || points.length < 4 || !points.every(p => Number.isFinite(p.x) && Number.isFinite(p.y)))
    throw new Error('Supply finite TE-to-TE contour points.');
  points = points.map(p => ({ ...p }));
  const knots = [0];
  for (let i = 1; i < points.length; i++) {
    const h = Math.hypot(points[i].x - points[i - 1].x, points[i].y - points[i - 1].y);
    if (!(h > 0)) throw new Error('Duplicate contour-curve point.'); knots.push(knots.at(-1) + h);
  }
  const length = knots.at(-1), count = points.length;
  if (!allowOpenEndpoints && Math.hypot(points[0].x - points.at(-1).x, points[0].y - points.at(-1).y) > 1e-10 * length)
    throw new Error('A sharp trailing edge requires a repeated closing point.');
  const second = key => {
    const diagonal = new Float64Array(count), rhs = new Float64Array(count), result = new Float64Array(count);
    diagonal[0] = diagonal[count - 1] = 1;
    for (let i = 1; i < count - 1; i++) {
      const left = knots[i] - knots[i - 1], right = knots[i + 1] - knots[i];
      diagonal[i] = 2 * (left + right);
      rhs[i] = 6 * ((points[i + 1][key] - points[i][key]) / right - (points[i][key] - points[i - 1][key]) / left);
      if (i > 1) { const factor = left / diagonal[i - 1]; diagonal[i] -= factor * left; rhs[i] -= factor * rhs[i - 1]; }
    }
    for (let i = count - 2; i > 0; i--) result[i] = (rhs[i] - (knots[i + 1] - knots[i]) * result[i + 1]) / diagonal[i];
    return result;
  };
  const curvature = { x: second('x'), y: second('y') };
  const evaluate = s => {
    if (!Number.isFinite(s) || s < 0 || s > length) throw new Error('Contour parameter is outside the curve.');
    let lo = 0, hi = count - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (knots[mid] > s) hi = mid; else lo = mid; }
    if (lo === count - 1) lo--;
    const h = knots[lo + 1] - knots[lo], b = (s - knots[lo]) / h, a = 1 - b;
    // Both coordinates share these basis values. Keep the original operation
    // order (including h * h / 6) to preserve thin-cell geometry bit for bit.
    const ac = a ** 3 - a, bc = b ** 3 - b, ad = 1 - 3 * a * a, bd = 3 * b * b - 1;
    const p0 = points[lo], p1 = points[lo + 1];
    const x0 = curvature.x[lo], x1 = curvature.x[lo + 1], y0 = curvature.y[lo], y1 = curvature.y[lo + 1];
    const point = { x: a * p0.x + b * p1.x + (ac * x0 + bc * x1) * h * h / 6,
      y: a * p0.y + b * p1.y + (ac * y0 + bc * y1) * h * h / 6 };
    const derivative = { x: (p1.x - p0.x) / h + (ad * x0 + bd * x1) * h / 6,
      y: (p1.y - p0.y) / h + (ad * y0 + bd * y1) * h / 6 };
    const secondDerivative = { x: a * x0 + b * x1, y: a * y0 + b * y1 };
    return { point, derivative, secondDerivative };
  };
  const branch = (side, fraction, stagnation) => {
    if (!['upper', 'lower'].includes(side) || !Number.isFinite(fraction) || fraction < 0 || fraction > 1
      || !(stagnation > 0 && stagnation < length)) throw new Error('Invalid stagnation/branch parameter.');
    const s = fraction === 1 ? (side === 'upper' ? 0 : length)
      : side === 'upper' ? (1 - fraction) * stagnation : stagnation + fraction * (length - stagnation);
    const value = evaluate(s), factor = 1 - fraction;
    return { ...value, parameter: s, stagnationDerivative: { x: factor * value.derivative.x, y: factor * value.derivative.y } };
  };
  return { length, knots, evaluate, branch };
}
