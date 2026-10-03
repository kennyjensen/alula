// SPDX-License-Identifier: GPL-2.0-or-later
// Specialize Spekreijse NLR TP 96735, Eqs. (12)--(18), to a prescribed
// harmonic map (S(xi,eta), eta). The mass coordinate remains harmonic.
// C_ab = -S_ab/S_xi; metric contraction gives J^2*Laplacian(xi).
// This is an inspectable elliptic control-map formulation, not MSET's X law.
export function harmonicMapCoefficients({ xi, xiXi, xiEta, etaEta }) {
  if (![xi, xiXi, xiEta, etaEta].every(Number.isFinite) || !(xi > 0))
    throw new Error('A streamwise harmonic map needs finite derivatives and positive S_xi.');
  const result = { xiXi: -xiXi / xi, xiEta: -xiEta / xi, etaEta: -etaEta / xi };
  if (!Object.values(result).every(Number.isFinite)) throw new Error('Unresolved harmonic-map coefficients.');
  return result;
}

export function createDiscreteHarmonicMapControl({ values, xi, eta, discretization = 'giles-1985' }) {
  const ordered = a => Array.isArray(a) && a.length >= 3 && a[0] === 0 && a.at(-1) === 1
    && a.every((v, i) => Number.isFinite(v) && (!i || v > a[i - 1]));
  if (!ordered(xi) || !ordered(eta) || !['giles-1985', 'quadratic'].includes(discretization)
    || !Array.isArray(values) || values.length !== xi.length || values.some(row => !Array.isArray(row) || row.length !== eta.length || !row.every(Number.isFinite))
    || values.some((row, i) => i && row.some((v, j) => !(v > values[i - 1][j]))))
    throw new Error('Harmonic-map samples require increasing coordinates and streamwise-ordered scalar values.');
  const first = (a, b, lo, center, hi) => discretization === 'giles-1985' ? (hi - lo) / (a + b)
    : (b * (center - lo) / a + a * (hi - center) / b) / (a + b);
  const second = (a, b, lo, center, hi) => 2 * ((hi - center) / b - (center - lo) / a) / (a + b);
  const coefficients = values.map(row => row.map(() => ({ xiXi: 0, xiEta: 0, etaEta: 0 })));
  for (let i = 1; i < xi.length - 1; i++) for (let j = 1; j < eta.length - 1; j++) {
    const a = xi[i] - xi[i - 1], b = xi[i + 1] - xi[i], c = eta[j] - eta[j - 1], d = eta[j + 1] - eta[j];
    const transverse = k => first(c, d, values[k][j - 1], values[k][j], values[k][j + 1]);
    coefficients[i][j] = harmonicMapCoefficients({
      xi: first(a, b, values[i - 1][j], values[i][j], values[i + 1][j]),
      xiXi: second(a, b, values[i - 1][j], values[i][j], values[i + 1][j]),
      xiEta: first(a, b, transverse(i - 1), transverse(i), transverse(i + 1)),
      etaEta: second(c, d, values[i][j - 1], values[i][j], values[i][j + 1]),
    });
  }
  return coefficients;
}
