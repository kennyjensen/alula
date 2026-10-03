// SPDX-License-Identifier: GPL-2.0-or-later
// C1 increasing Hermite interpolation. The default harmonic rule retains
// positive derivatives and secant endpoint slopes for existing inverse maps.
// Optional PCHIP uses the increasing-data subset of SLATEC DPCHIM:
// https://www.netlib.org/slatec/pchip/dpchim.f (Fritsch/Butland derivatives).
// Its endpoint derivatives can be zero. This is generic interpolation, not
// a recovered MSET X-spacing law. Prescribed positive slopes are accepted
// only in the sufficient monotonicity cone: at most twice either adjacent
// secant. These interpolants are not generally C2.
export function createMonotoneCubicMap(x, y, options = {}) {
  if (!options || typeof options !== 'object' || Array.isArray(options)
    || Object.keys(options).some(key => !['derivatives', 'slopes'].includes(key))
    || !['harmonic', 'pchip', 'prescribed'].includes(options.derivatives === undefined ? 'harmonic' : options.derivatives))
    throw new Error('Invalid monotone cubic derivative options.');
  const derivatives = options.derivatives === undefined ? 'harmonic' : options.derivatives;
  if (derivatives !== 'prescribed' && Object.hasOwn(options, 'slopes'))
    throw new Error('Prescribed slopes require the prescribed derivative option.');
  if (!Array.isArray(x) || !Array.isArray(y) || x.length < 2 || x.length !== y.length
    || x.some((v, i) => !Number.isFinite(v) || i && v <= x[i - 1])
    || y.some((v, i) => !Number.isFinite(v) || i && v <= y[i - 1]))
    throw new Error('Monotone cubic data must be finite and strictly increasing.');
  x = x.slice(); y = y.slice();
  const secants = x.slice(1).map((v, i) => (y[i + 1] - y[i]) / (v - x[i]));
  if (secants.some(v => !(v > 0) || !Number.isFinite(v))) throw new Error('Unresolved monotone cubic secant.');
  const slopes = x.map((_, i) => !i ? secants[0] : i === x.length - 1 ? secants.at(-1)
    : Math.min(secants[i - 1], secants[i]) * (2 / (1 + Math.min(secants[i - 1], secants[i]) / Math.max(secants[i - 1], secants[i]))));
  if (derivatives === 'prescribed') {
    if (!Array.isArray(options.slopes) || options.slopes.length !== x.length
      || Array.from(options.slopes).some((slope, i) => !Number.isFinite(slope) || !(slope > 0)
        || (i > 0 && slope / secants[i - 1] > 2) || (i < secants.length && slope / secants[i] > 2)))
      throw new Error('Prescribed cubic slopes must be finite, positive, and inside the two-secant monotonicity cone.');
    options.slopes.forEach((slope, i) => { slopes[i] = slope; });
  }
  if (derivatives === 'pchip' && x.length > 2) {
    const widths = x.slice(1).map((v, i) => v - x[i]);
    for (let i = 1; i < x.length - 1; i++) {
      const left = secants[i - 1], right = secants[i];
      if (left === right) { slopes[i] = left; continue; }
      const hmax = Math.max(widths[i - 1], widths[i]), h1 = widths[i - 1] / hmax, h2 = widths[i] / hmax;
      const w1 = (2 * h1 + h2) / (3 * (h1 + h2)), w2 = (h1 + 2 * h2) / (3 * (h1 + h2));
      const dmax = Math.max(left, right), dmin = Math.min(left, right);
      // DPCHIM's scaled form avoids reciprocals of tiny secants. It equals
      // (w1+w2)/(w2/left+w1/right), with interval-dependent weights.
      slopes[i] = dmin / (w1 * (left / dmax) + w2 * (right / dmax));
    }
    const endpoint = (h1, h2, d1, d2) => {
      const largest = Math.max(h1, h2), ratio = (h1 / largest) / (h1 / largest + h2 / largest);
      // One-sided three-point derivative. All input secants are positive,
      // so DPCHIM's sign-switch 3*d1 limiter is inapplicable to this subset.
      return Math.max(0, d1 + ratio * (d1 - d2));
    };
    slopes[0] = endpoint(widths[0], widths[1], secants[0], secants[1]);
    slopes[slopes.length - 1] = endpoint(widths.at(-1), widths.at(-2), secants.at(-1), secants.at(-2));
    if (slopes.some(v => !Number.isFinite(v) || v < 0)) throw new Error('Unresolved monotone cubic PCHIP derivative.');
  }
  const evaluate = value => {
    if (!Number.isFinite(value) || value < x[0] || value > x.at(-1)) throw new Error('Monotone cubic argument lies outside its knots.');
    let lo = 0, hi = x.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (x[mid] <= value) lo = mid; else hi = mid; }
    const width = x[hi] - x[lo], t = (value - x[lo]) / width, a = slopes[lo] / secants[lo], b = slopes[hi] / secants[lo];
    const f = t * (a + t * (3 - 2 * a - b + t * (a + b - 2)));
    const derivative = secants[lo] * (a + t * (2 * (3 - 2 * a - b) + 3 * t * (a + b - 2)));
    // At an interior knot, use the interval on its right, as for the value
    // and first derivative. The final endpoint uses the interval on its left.
    const secondDerivative = secants[lo] / width * (2 * (3 - 2 * a - b) + 6 * t * (a + b - 2));
    return { value: value === x[lo] ? y[lo] : value === x[hi] ? y[hi] : y[lo] + (y[hi] - y[lo]) * f, derivative, secondDerivative };
  };
  return { knots: x, values: y, slopes, evaluate, value: x => evaluate(x).value };
}
