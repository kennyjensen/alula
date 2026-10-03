// SPDX-License-Identifier: GPL-2.0-or-later
// Vinokur, NASA-CR-3313, printed p. 15, Eqs. (54)-(61): one-sided
// stretching with prescribed initial slope and zero final curvature.
// This is the published spacing function, not a recovered MSET X law.
const MIN_NORMAL = 2 ** -1022;

function logShapeRatio(a, hyperbolic) {
  const x = 2 * a, x2 = x * x;
  if (x === 0) return 0;
  if (x < 0.1) {
    // log(sinh(x)/x), or -log(sin(x)/x), without loss near x=0.
    const sign = hyperbolic ? -1 : 1;
    return x2 * (1 / 6 + sign * x2 * (1 / 180 + sign * x2 * (1 / 2835
      + sign * x2 * (1 / 37800 + sign * x2 / 467775))));
  }
  if (!hyperbolic) return -Math.log(Math.sin(x) / x);
  return x < 20 ? Math.log(Math.sinh(x) / x)
    : x + Math.log(-Math.expm1(-2 * x)) - Math.log(2 * x);
}

function logCosh(x) {
  return x < 20 ? Math.log1p(2 * Math.sinh(x / 2) ** 2)
    : x + Math.log1p(Math.exp(-2 * x)) - Math.LN2;
}
function logCos(x) {
  return x < 0.25 ? Math.log1p(-2 * Math.sin(x / 2) ** 2) : Math.log(Math.cos(x));
}
function logSinh(x) {
  return x < 20 ? Math.log(Math.sinh(x)) : x + Math.log(-Math.expm1(-2 * x)) - Math.LN2;
}
function scaledPositive(scale, logFactor) {
  const factor = Math.exp(logFactor), direct = scale * factor;
  const value = factor >= MIN_NORMAL && direct > 0 && Number.isFinite(direct)
    ? direct : Math.exp(Math.log(scale) + logFactor);
  if (!(value > 0) || !Number.isFinite(value))
    throw new Error('Vinokur stretch value or derivative is not representable.');
  return value;
}

export function createVinokurStretch({ length, logicalWidth, initialDerivative } = {}) {
  if (![length, logicalWidth, initialDerivative].every(v => Number.isFinite(v) && v > 0))
    throw new Error('Vinokur length, logical width, and initial derivative must be finite and positive.');
  const mean = length / logicalWidth;
  if (!(mean > 0) || !Number.isFinite(mean)) throw new Error('Vinokur mean derivative is not representable.');
  const branch = initialDerivative === mean ? 'affine' : initialDerivative < mean ? 'hyperbolic' : 'trigonometric';
  const hyperbolic = branch === 'hyperbolic';
  let parameter = 0;
  if (branch !== 'affine') {
    const larger = Math.max(mean, initialDerivative), smaller = Math.min(mean, initialDerivative);
    const relative = (larger - smaller) / smaller;
    const target = relative < 0.5 ? Math.log1p(relative) : Math.log(larger) - Math.log(smaller);
    let lo = 0, hi = hyperbolic ? 1 : Math.PI / 2 - Number.EPSILON;
    // Hyperbolic log ratios grow without bound. Sixteen doublings exceed
    // every ratio of positive finite doubles; trigonometric a must stay <pi/2.
    if (hyperbolic) for (let i = 0; i < 16 && logShapeRatio(hi, true) < target; i++) hi *= 2;
    if (!(target > 0) || logShapeRatio(hi, hyperbolic) < target)
      throw new Error('Vinokur requested slope has no resolved representable parameter.');
    let solved = false;
    for (let i = 0; i < 128; i++) {
      const mid = lo + (hi - lo) / 2;
      if (mid === lo || mid === hi) {
        parameter = Math.abs(logShapeRatio(lo, hyperbolic) - target) <= Math.abs(logShapeRatio(hi, hyperbolic) - target) ? lo : hi;
        solved = true; break;
      }
      const residual = logShapeRatio(mid, hyperbolic) - target;
      if (residual === 0) { parameter = mid; solved = true; break; }
      if (residual < 0) lo = mid; else hi = mid;
    }
    if (!solved || !(parameter > 0)
      || Math.abs(logShapeRatio(parameter, hyperbolic) - target) > 64 * Number.EPSILON * Math.max(1, target))
      throw new Error('Vinokur parameter solve cannot resolve the requested initial slope.');
  }
  const a = parameter;
  const derivativeAt = t => {
    if (t === 0 || branch === 'affine') return initialDerivative;
    const u = a * (1 - t);
    const logRatio = hyperbolic ? 2 * (logCosh(a) - logCosh(u)) : 2 * (logCos(a) - logCos(u));
    return scaledPositive(initialDerivative, logRatio);
  };
  const secondAt = t => {
    if (t === 1 || branch === 'affine') return 0;
    const u = a * (1 - t), factor = 2 * a * (hyperbolic ? Math.tanh(u) : Math.tan(u));
    const derivative = derivativeAt(t), product = derivative * factor, direct = product / logicalWidth;
    const magnitude = product >= MIN_NORMAL && direct > 0 && Number.isFinite(direct)
      ? direct : Math.exp(Math.log(derivative) + Math.log(factor) - Math.log(logicalWidth));
    if (!(magnitude > 0) || !Number.isFinite(magnitude)) throw new Error('Vinokur second derivative is not representable.');
    return hyperbolic ? magnitude : -magnitude;
  };
  // The derivative is monotone. Hyperbolic curvature can have an interior
  // maximum, at tanh(a*(1-t))=1/sqrt(3); trigonometric curvature peaks at 0.
  derivativeAt(1); secondAt(0);
  const critical = Math.atanh(1 / Math.sqrt(3));
  if (hyperbolic && a > critical) secondAt(1 - critical / a);
  const value = t => {
    if (!Number.isFinite(t) || t < 0 || t > 1) throw new Error('Vinokur argument must lie in [0,1].');
    if (t === 0) return 0;
    if (t === 1) return length;
    if (branch === 'affine') {
      const affine = length * t;
      return affine > 0 && Number.isFinite(affine) ? affine : scaledPositive(length, Math.log(t));
    }
    const u = a * (1 - t), at = a * t;
    let fraction, logFraction;
    if (hyperbolic && a >= 20) {
      logFraction = logSinh(at) - logCosh(u) - logSinh(a);
    } else {
      // Multiplying t by sinc/sinhc factors avoids both cancellation at the
      // near endpoint and division of two underflowed a*t factors.
      const numerator = at === 0 ? 1 : (hyperbolic ? Math.sinh(at) : Math.sin(at)) / at;
      const denominator = (hyperbolic ? Math.sinh(a) : Math.sin(a)) / a;
      const cosine = hyperbolic ? Math.cosh(u) : Math.cos(u);
      const ratio = numerator / denominator / cosine;
      fraction = t * ratio;
      logFraction = Math.log(t) + Math.log(ratio);
    }
    const direct = length * fraction;
    return fraction >= MIN_NORMAL && direct > 0 && Number.isFinite(direct) ? direct : scaledPositive(length, logFraction);
  };
  const evaluate = t => ({ value: value(t), derivative: derivativeAt(t), secondDerivative: secondAt(t) });
  // Positive analytic derivatives imply strict monotonicity. Consumers must
  // still check that their sampled intervals are distinguishable as doubles.
  return { parameter, branch, value, evaluate };
}
