// SPDX-License-Identifier: GPL-2.0-or-later
// In local linear stagnation flow, phi = K(t^2-n^2)/2 and psi = K*t*n.
// Along the cross-line phi-phi_s = slope*(psi-psi_s), the distance R from
// stagnation gives rho*|delta psi| = rho*K*R^2/(2*sqrt(1+slope^2)).
export function stagnationStreamtubeMass({ strainRate, distance, potentialSlope = 0, density = 1 }) {
  if (![strainRate, distance, potentialSlope, density].every(Number.isFinite)
    || Math.min(strainRate, distance, density) <= 0) throw new Error('Invalid stagnation-flow spacing controls.');
  const mass = .5 * density * strainRate * distance * distance / Math.hypot(1, potentialSlope);
  if (!(mass > 0) || !Number.isFinite(mass)) throw new Error('Unresolved stagnation streamtube mass.');
  return mass;
}

// Positive mass fractions with prescribed first and/or last fraction. One
// wall uses a geometric progression; two walls use a quadratic log profile.
// The one remaining scalar determines the total, leaving endpoint masses
// fixed. This is a spacing construction, not an equation in the Euler solve.
export function createNormalMassWeights({ count, first, last }) {
  const constrained = [first, last].filter(v => v !== undefined);
  if (!Number.isInteger(count) || count < 3 || !constrained.length
    || constrained.some(v => !Number.isFinite(v) || v <= 0 || v >= 1)
    || constrained.reduce((a, b) => a + b, 0) >= 1) throw new Error('Incompatible normal streamtube endpoint fractions.');
  const both = first !== undefined && last !== undefined;
  const terms = Array.from({ length: count }, (_, i) => {
    const t = i / (count - 1);
    return both ? { a: (1 - t) * Math.log(first) + t * Math.log(last), b: 4 * t * (1 - t) }
      : { a: Math.log(first ?? last), b: first === undefined ? count - 1 - i : i };
  });
  const logarithms = lambda => terms.map(p => p.a + lambda * p.b);
  const logTotal = lambda => {
    const a = logarithms(lambda), largest = Math.max(...a);
    return largest + Math.log(a.reduce((sum, v) => sum + Math.exp(v - largest), 0));
  };
  let lo = -1024, hi = 1024;
  if (!(logTotal(lo) < 0 && logTotal(hi) > 0)) throw new Error('Normal mass spacing cannot resolve the requested endpoints.');
  for (let i = 0; i < 80; i++) { const mid = .5 * (lo + hi); if (logTotal(mid) < 0) lo = mid; else hi = mid; }
  const weights = logarithms(.5 * (lo + hi)).map(Math.exp);
  const free = weights.map((_, i) => i).filter(i => !(i === 0 && first !== undefined) && !(i === count - 1 && last !== undefined));
  const remaining = 1 - constrained.reduce((a, b) => a + b, 0), sum = free.reduce((a, i) => a + weights[i], 0);
  for (const i of free) weights[i] *= remaining / sum;
  if (first !== undefined) weights[0] = first; if (last !== undefined) weights[count - 1] = last;
  if (weights.some(v => !(v > 0) || !Number.isFinite(v))) throw new Error('Unresolved normal mass fraction.');
  return weights;
}

// Add resolution to meet a bound on adjacent mass ratios while preserving
// both measured wall fractions and total channel mass. No endpoint clipping.
export function resolveNormalMassWeights({ count, first, last, maximumAdjacentRatio = Infinity, maximumTubes = 63 }) {
  if (!(maximumAdjacentRatio > 1) || !Number.isInteger(maximumTubes) || maximumTubes < count)
    throw new Error('Invalid normal-spacing refinement controls.');
  for (let n = count; n <= maximumTubes; n++) {
    const weights = createNormalMassWeights({ count: n, first, last });
    const maxAdjacentRatio = Math.max(...weights.slice(1).map((v, i) => Math.max(v / weights[i], weights[i] / v)));
    if (maxAdjacentRatio <= maximumAdjacentRatio) return { weights, maxAdjacentRatio, requestedTubes: count, addedTubes: n - count };
  }
  throw new Error('Stagnation aspect target requires more transverse resolution than the grid budget.');
}
