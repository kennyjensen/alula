// SPDX-License-Identifier: GPL-2.0-or-later
// Resolve contour turning before matching multielement grid cross-lines.
// The derivative of each cubic spline piece is a quadratic Bezier curve.
// Its three control vectors bound the tangent cone by convexity; checking
// that cone avoids missing a high-curvature nose between endpoint samples.
export function resolveSurfaceTurns({ curve, side, stagnation, fractions, maxTurn = Math.PI / 12, maxPoints = 4097 }) {
  if (!curve || !['upper', 'lower'].includes(side) || !Array.isArray(fractions) || fractions[0] !== 0 || fractions.at(-1) !== 1
    || fractions.some((f, i) => !Number.isFinite(f) || i && f <= fractions[i - 1])
    || !(stagnation > 0 && stagnation < curve.length) || !Number.isFinite(maxTurn) || !(maxTurn > 0 && maxTurn < Math.PI / 2)
    || !Number.isInteger(maxPoints) || maxPoints < fractions.length) throw new Error('Invalid surface-turn resolution controls.');
  const parameter = f => side === 'upper' ? (1 - f) * stagnation : stagnation + f * (curve.length - stagnation);
  const cone = (f0, f1) => {
    const a = Math.min(parameter(f0), parameter(f1)), b = Math.max(parameter(f0), parameter(f1));
    const knots = [a, ...curve.knots.filter(s => s > a && s < b), b], base = curve.evaluate(a).derivative;
    let lo = 0, hi = 0;
    for (let k = 1; k < knots.length; k++) {
      const u = curve.evaluate(knots[k - 1]), v = curve.evaluate(knots[k]), half = .5 * (knots[k] - knots[k - 1]);
      const middle = { x: u.derivative.x + half * u.secondDerivative.x, y: u.derivative.y + half * u.secondDerivative.y };
      for (const d of [u.derivative, middle, v.derivative]) {
        if (!(Math.hypot(d.x, d.y) > 0)) return Infinity;
        const angle = Math.atan2(base.x * d.y - base.y * d.x, base.x * d.x + base.y * d.y);
        lo = Math.min(lo, angle); hi = Math.max(hi, angle);
      }
    }
    return hi - lo;
  };
  const resolved = [0]; let worst = 0;
  const interval = (a, b, depth) => {
    const turn = cone(a, b);
    if (turn <= maxTurn) { resolved.push(b); worst = Math.max(worst, turn); }
    else {
      if (depth >= 40 || resolved.length >= maxPoints - 1) throw new Error('Surface-turn refinement budget exhausted.');
      const mid = .5 * (a + b); interval(a, mid, depth + 1); interval(mid, b, depth + 1);
    }
    if (resolved.length > maxPoints) throw new Error('Surface-turn refinement budget exhausted.');
  };
  for (let i = 1; i < fractions.length; i++) interval(fractions[i - 1], fractions[i], 0);
  return { fractions: resolved, addedPoints: resolved.length - fractions.length, maximumTangentCone: worst, maxTurn };
}
