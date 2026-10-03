// SPDX-License-Identifier: GPL-2.0-or-later
// Physical arc coordinate on the C2 parametric contour. Integrate separately
// on each spline interval; the supplied chord-length parameter is not arc.
const abscissae = [.1834346424956498, .525532409916329, .7966664774136267, .9602898564975363];
const weights = [.362683783378362, .3137066458778873, .2223810344533745, .1012285362903763];
export function createContourArc(curve) {
  const speed = s => { const d = curve.evaluate(s).derivative; return Math.hypot(d.x, d.y); };
  const integral = (a, b) => {
    const mid = .5 * (a + b), half = .5 * (b - a); let sum = 0;
    for (let i = 0; i < 4; i++) sum += weights[i] * (speed(mid - half * abscissae[i]) + speed(mid + half * abscissae[i]));
    return half * sum;
  };
  const arc = [0]; for (let i = 1; i < curve.knots.length; i++) arc.push(arc.at(-1) + integral(curve.knots[i - 1], curve.knots[i]));
  const at = s => {
    if (!Number.isFinite(s) || s < 0 || s > curve.length) throw new Error('Contour arc parameter is outside the curve.');
    let lo = 0, hi = curve.knots.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (curve.knots[mid] > s) hi = mid; else lo = mid; }
    return arc[lo] + integral(curve.knots[lo], s);
  };
  return { at, speed, length: arc.at(-1) };
}
