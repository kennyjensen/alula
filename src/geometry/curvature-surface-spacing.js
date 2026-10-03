// SPDX-License-Identifier: GPL-2.0-or-later
// MSES manual 2.2.2: ds ~ 1/(1 + a |kappa|^b), calibrated by
// physical LE spacing / (half perimeter / number of side points).
// Equidistribution, adaptive quadrature and the optional compact TE density
// bumps are reconstruction choices; no MSET source parity is claimed.
import { createContourArc } from './contour-arc.js';

const gx = [.1834346424956498, .525532409916329, .7966664774136267, .9602898564975363];
const gw = [.362683783378362, .3137066458778873, .2223810344533745, .1012285362903763];
const gauss = (fn, a, b) => {
  const mid = .5 * (a + b), half = .5 * (b - a), sum = [0, 0];
  for (let k = 0; k < 4; k++) for (const sign of [-1, 1]) {
    const v = fn(mid + sign * half * gx[k]);
    for (let j = 0; j < 2; j++) sum[j] += gw[k] * v[j];
  }
  return sum.map(v => half * v);
};
const integrate = (fn, knots, tolerance) => {
  const piece = (a, b) => {
    const coarse = gauss(fn, a, b), mid = .5 * (a + b), left = gauss(fn, a, mid), right = gauss(fn, mid, b);
    const value = left.map((v, k) => v + right[k]);
    return { a, b, value, error: Math.max(...value.map((v, k) => Math.abs(v - coarse[k]))) };
  };
  const pieces = knots.slice(1).map((b, k) => piece(knots[k], b));
  let sum = [0, 0], error = 0;
  for (;;) {
    sum = [0, 0]; error = 0; let worst = 0;
    pieces.forEach((p, k) => { p.value.forEach((v, j) => { sum[j] += v; }); error += p.error; if (p.error > pieces[worst].error) worst = k; });
    if (error <= tolerance * Math.max(...sum)) break;
    if (pieces.length >= 8192) throw new Error('Curvature spacing quadrature budget exhausted.');
    const p = pieces[worst], mid = .5 * (p.a + p.b);
    if (mid === p.a || mid === p.b) throw new Error('Curvature spacing quadrature reached roundoff.');
    pieces.splice(worst, 1, piece(p.a, mid), piece(mid, p.b));
  }
  let cumulative = [0, 0];
  for (const p of pieces) { p.cumulative = cumulative; cumulative = cumulative.map((v, k) => v + p.value[k]); }
  const at = s => {
    if (s === knots.at(-1)) return sum.slice();
    let lo = 0, hi = pieces.length - 1;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (pieces[mid].b <= s) lo = mid + 1; else hi = mid; }
    const p = pieces[lo];
    if (s === p.a) return p.cumulative.slice();
    const value = gauss(fn, p.a, s);
    return value.map((v, k) => v + p.cumulative[k]);
  };
  return { at, total: sum, error, pieces: pieces.length };
};

export function distributeCurvatureSurface({ curve, side, stagnation, count = 65, exponent = .5,
  leadingSpacingRatio, trailingSpacingRatio, curvatureWeight, artificialLeadingCurvature = true, tolerance = 1e-11 }) {
  if (!curve || !Array.isArray(curve.knots) || typeof curve.evaluate !== 'function' || !['upper', 'lower'].includes(side)
    || !Number.isInteger(count) || count < 4 || count > 4097 || !Number.isFinite(exponent) || exponent < 0 || exponent > 4
    || !(stagnation > 0 && stagnation < curve.length) || !Number.isFinite(tolerance) || tolerance <= 0 || tolerance > 1e-4)
    throw new Error('Invalid curvature surface-spacing inputs.');
  const prescribed = curvatureWeight !== undefined;
  if (prescribed ? !Number.isFinite(curvatureWeight) || curvatureWeight < 0 || leadingSpacingRatio !== undefined || trailingSpacingRatio !== undefined
    : !Number.isFinite(leadingSpacingRatio) || leadingSpacingRatio <= 0
      || (trailingSpacingRatio !== undefined && (!Number.isFinite(trailingSpacingRatio) || trailingSpacingRatio <= 0)))
    throw new Error('Supply a nonnegative curvature weight or positive endpoint spacing ratios.');
  const physicalArc = createContourArc(curve), referenceLength = .5 * physicalArc.length, averageSpacing = referenceLength / count;
  const first = side === 'upper' ? 0 : stagnation, last = side === 'upper' ? stagnation : curve.length;
  const metric = s => {
    const { derivative: d, secondDerivative: dd } = curve.evaluate(s), speed = Math.hypot(d.x, d.y);
    const curvature = Math.abs(d.x * dd.y - d.y * dd.x) / speed ** 3;
    const weight = (referenceLength * curvature) ** exponent;
    if (!(speed > 0) || ![speed, weight].every(Number.isFinite)) throw new Error('Singular curvature surface metric.');
    return [speed, speed * weight];
  };
  const knots = [first, ...curve.knots.filter(s => s > first && s < last), last];
  const measure = integrate(metric, knots, tolerance), [length, curvatureIntegral] = measure.total;
  const intervals = count - 1;
  const parameter = f => f === 1 ? (side === 'upper' ? 0 : curve.length)
    : side === 'upper' ? (1 - f) * stagnation : stagnation + f * (curve.length - stagnation);
  const at = f => {
    const v = measure.at(parameter(f));
    return side === 'upper' ? v.map((value, k) => measure.total[k] - value) : v;
  };
  const invertArc = distance => {
    let lo = 0, hi = 1;
    for (let k = 0; k < 55; k++) { const mid = .5 * (lo + hi); if (at(mid)[0] < distance) lo = mid; else hi = mid; }
    return .5 * (lo + hi);
  };
  let a = curvatureWeight ?? 0, teWeight = 0, teWidth = 0, leWeight = 0, leWidth = 0;
  // A C3 compact endpoint bump: B(t)=(1-t)^4(1+4t), 0<=t<=1.
  // Its exact physical-arc integral avoids quadrature through artificial
  // curvature. Width and shape are documented reconstruction choices.
  const primitive = t => t - (10 / 3) * t ** 3 + 5 * t ** 4 - 3 * t ** 5 + (2 / 3) * t ** 6;
  const bumpIntegral = distance => !teWidth || distance <= length - teWidth ? 0
    : distance >= length ? teWidth / 3 : teWidth * (1 / 3 - primitive((length - distance) / teWidth));
  const leadingIntegral = distance => !leWidth ? 0 : leWidth * primitive(Math.min(1, Math.max(0, distance / leWidth)));
  if (!prescribed) {
    const le = leadingSpacingRatio * averageSpacing, te = (trailingSpacingRatio ?? 0) * averageSpacing;
    if (!(le > 0 && le < length / 2) || (trailingSpacingRatio !== undefined && !(te > 0 && te < length / 2)))
      throw new Error('Endpoint spacing consumes too much of the surface branch.');
    const kle = at(invertArc(le))[1], left = intervals * kle - curvatureIntegral, rhs = length - intervals * le;
    if (trailingSpacingRatio === undefined) {
      if (Math.abs(rhs) <= tolerance * length) a = 0;
      else if (Math.abs(left) > tolerance * Math.max(length, curvatureIntegral)) a = rhs / left;
      else a = NaN;
    } else {
      teWidth = Math.min(.25 * length, 4 * te);
      const kte = curvatureIntegral - at(invertArc(length - te))[1], total = teWidth / 3;
      const right = intervals * bumpIntegral(le) - total, lowerLeft = intervals * kte - curvatureIntegral;
      const lowerRight = intervals * (total - bumpIntegral(length - te)) - total, lowerRhs = length - intervals * te;
      const determinant = left * lowerRight - right * lowerLeft;
      if (Math.abs(determinant) > tolerance * (Math.abs(left * lowerRight) + Math.abs(right * lowerLeft))) {
        a = (rhs * lowerRight - right * lowerRhs) / determinant;
        teWeight = (left * lowerRhs - rhs * lowerLeft) / determinant;
      } else a = NaN;
    }
    if (![a, teWeight].every(v => Number.isFinite(v) && v >= 0)) {
      if (!artificialLeadingCurvature) throw new Error('Endpoint spacing requires artificial leading-edge curvature.');
      // Low curvature at stagnation can make the pure curvature calibration
      // impossible. The manual permits a local artificial LE contribution.
      // Keep the physical-curvature monitor, and solve nonnegative endpoint
      // bump weights exactly. The remaining degree of freedom chooses equal
      // integrated base/curvature density, clipped to the feasible interval.
      leWidth = Math.min(.25 * length, 4 * le);
      const ltotal = leWidth / 3, lleft = intervals * leadingIntegral(le) - ltotal;
      let base, slope;
      if (trailingSpacingRatio === undefined) {
        base = [rhs / lleft, 0]; slope = [-left / lleft, 0];
      } else {
        const ttotal = teWidth / 3, tleft = intervals * bumpIntegral(le) - ttotal;
        const lright = intervals * (ltotal - leadingIntegral(length - te)) - ltotal;
        const tright = intervals * (ttotal - bumpIntegral(length - te)) - ttotal;
        const det = lleft * tright - tleft * lright;
        if (!(Math.abs(det) > tolerance * (Math.abs(lleft * tright) + Math.abs(tleft * lright))))
          throw new Error('Endpoint density-spacing constraints are singular.');
        const solve = (u, v) => [(u * tright - tleft * v) / det, (lleft * v - u * lright) / det];
        base = solve(rhs, length - intervals * te);
        slope = solve(-left, curvatureIntegral - intervals * (curvatureIntegral - at(invertArc(length - te))[1]));
      }
      let low = 0, high = Infinity;
      for (let k = 0; k < 2; k++) {
        if (slope[k] > 0) low = Math.max(low, -base[k] / slope[k]);
        else if (slope[k] < 0) high = Math.min(high, -base[k] / slope[k]);
        else if (base[k] < 0) high = -1;
      }
      if (!(low <= high) || !Number.isFinite(low)) throw new Error('Requested endpoint spacing has no nonnegative density distribution.');
      a = Math.max(low, Math.min(high, curvatureIntegral > 0 ? length / curvatureIntegral : 0));
      [leWeight, teWeight] = base.map((v, k) => Math.max(0, v + a * slope[k]));
      if (![a, leWeight, teWeight].every(Number.isFinite)) throw new Error('Unresolved artificial endpoint curvature.');
    }
  }
  const total = length + a * curvatureIntegral + teWeight * teWidth / 3 + leWeight * leWidth / 3;
  const weightAt = f => { const [arc, curvature] = at(f); return arc + a * curvature + teWeight * bumpIntegral(arc) + leWeight * leadingIntegral(arc); };
  const fractions = Array.from({ length: count }, (_, i) => {
    if (!i) return 0; if (i === intervals) return 1;
    let lo = 0, hi = 1; const target = total * i / intervals;
    for (let k = 0; k < 55; k++) { const mid = .5 * (lo + hi); if (weightAt(mid) < target) lo = mid; else hi = mid; }
    return .5 * (lo + hi);
  });
  if (!(total > 0) || !Number.isFinite(total) || fractions.some((f, i) => i && f <= fractions[i - 1]))
    throw new Error('Unresolved curvature surface distribution.');
  const arcLengths = fractions.map(f => at(f)[0]), spacings = arcLengths.slice(1).map((v, i) => v - arcLengths[i]);
  return { fractions, parameters: fractions.map(parameter), arcLengths, spacings,
    diagnostics: { method: 'curvature equidistribution', side, count, exponent, length, referenceLength, averageSpacing,
      curvatureWeight: a, dimensionalCurvatureCoefficient: a * referenceLength ** exponent, leadingSpacingRatio, trailingSpacingRatio,
      actualLeadingSpacingRatio: spacings[0] / averageSpacing, actualTrailingSpacingRatio: spacings.at(-1) / averageSpacing,
      artificialLeadingCurvature: leWeight > 0, leadingBump: { weight: leWeight, width: leWidth, shape: 'compact C3 (1-t)^4(1+4t)' },
      trailingBump: { weight: teWeight, width: teWidth, shape: 'compact C3 (1-t)^4(1+4t)' },
      quadratureErrorEstimate: measure.error, quadraturePieces: measure.pieces, totalMetric: total } };
}

// Preserve the point-density distribution of the supplied contour, splitting
// the continuous source index at the panel stagnation point. Refine the index
// linearly, not with another cosine map in arclength: the source may already
// cluster its points at the nose. Geometry stays on the original spline.
export function distributeSourceSurface({ curve, side, stagnation, count = 65 }) {
  const knots = curve?.knots;
  if (!Array.isArray(knots) || knots.length < 4 || knots[0] !== 0
    || knots.at(-1) !== curve.length || !knots.every(Number.isFinite)
    || knots.some((s, i) => i && s <= knots[i - 1])
    || !['upper', 'lower'].includes(side) || !Number.isFinite(stagnation)
    || !(stagnation > 0 && stagnation < curve.length)
    || !Number.isInteger(count) || count < 4 || count > 4097)
    throw new Error('Invalid source surface-spacing inputs.');
  let k = 0;
  while (knots[k + 1] < stagnation) k++;
  const rank = k + (stagnation - knots[k]) / (knots[k + 1] - knots[k]);
  return Array.from({ length: count }, (_, i) => {
    if (!i) return 0;
    if (i === count - 1) return 1;
    const f = i / (count - 1);
    const u = side === 'upper' ? rank * (1 - f) : rank + (knots.length - 1 - rank) * f;
    const j = Math.floor(u), s = knots[j] + (u - j) * (knots[j + 1] - knots[j]);
    return side === 'upper' ? (stagnation - s) / stagnation : (s - stagnation) / (curve.length - stagnation);
  });
}
