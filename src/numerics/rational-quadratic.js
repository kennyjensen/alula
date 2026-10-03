// SPDX-License-Identifier: GPL-2.0-or-later
// Normalized monotone rational-quadratic Hermite segment (Gregory and
// Delbourgo): endpoints (0,0),(1,1), nonnegative endpoint derivatives a,b.
// Its positive Bernstein denominator and derivative numerator preserve
// ordering even when a monotone cubic with the same derivatives would not.
const validate = (t, a, b) => {
  if (![t, a, b].every(Number.isFinite) || t < 0 || t > 1 || a < 0 || b < 0)
    throw new Error('Invalid monotone rational-quadratic segment.');
};
export function monotoneRationalQuadratic(t, a, b) {
  validate(t, a, b);
  if (t === 0 || t === 1) return t;
  const u = 1 - t, scale = Math.max(1, a, b);
  const numerator = t * t / scale + (a / scale) * t * u;
  const remainder = u * u / scale + (b / scale) * t * u;
  return numerator / (numerator + remainder);
}
export function inverseMonotoneRationalQuadratic(y, a, b) {
  validate(y, a, b);
  if (y === 0 || y === 1) return y;
  // A*t^2+B*t-y=0, A+B=1. The discriminant identity avoids subtracting
  // nearly equal squares; choose a root formula without cancellation.
  const B = a * (1 - y) - b * y + 2 * y;
  const root = Math.hypot(B - 2 * y, 2 * Math.sqrt(y * (1 - y)));
  const t = B >= 0 ? y / (.5 * B + .5 * root) : (.5 * root - .5 * B) / (1 - B);
  if (!(t > 0 && t < 1) || !Number.isFinite(t)) throw new Error('Unresolved rational-quadratic inverse.');
  return t;
}
