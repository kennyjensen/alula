// SPDX-License-Identifier: GPL-2.0-or-later
// Second-order local derivative inside a known scalar domain. Near a domain
// edge, a one-sided stencil retains the nominal step without probing another
// active equation branch or subtracting at a vanishingly small central step.
export function boundedResidualPartial(evaluate, value, { step, lower = -Infinity, upper = Infinity, base } = {}) {
  if (!Number.isFinite(value) || !Number.isFinite(step) || !(step > 0) || !(lower < upper) || value < lower || value > upper)
    throw new Error('Invalid bounded residual derivative domain.');
  const left = value - lower, right = upper - value;
  let a, b, result;
  if (left >= step && right >= step) {
    a = evaluate(value + step); b = evaluate(value - step);
    result = a.map((v, i) => (v - b[i]) / (2 * step));
  } else {
    const sign = right >= left ? 1 : -1, room = sign > 0 ? right : left;
    const h = sign * Math.min(step, room / 3);
    if (value + h === value || value + 2 * h === value + h)
      throw new Error('Unresolved bounded residual derivative step.');
    const zero = base ?? evaluate(value);
    a = evaluate(value + h); b = evaluate(value + 2 * h);
    result = a.map((v, i) => (4 * (v - zero[i]) - (b[i] - zero[i])) / (2 * h));
  }
  if (!result.length || !result.every(Number.isFinite)) throw new Error('Nonfinite bounded residual derivative.');
  return result;
}
