// SPDX-License-Identifier: GPL-2.0-or-later
// Independent smooth, isentropic Ringleb solution, gamma=7/5 only.
// Sources and branch/derivative derivation: docs/RINGLEB_LOCAL_VERIFICATION.md.
// No production imports, coordinate inversion, grid solve or fitted gas state.
export function ringlebPoint(q, psi) {
  if (![q, psi].every(Number.isFinite) || !(q > 0 && psi > 0 && q * psi < 1))
    throw new Error('Ringleb upper branch requires q>0, psi>0 and q*psi<1.');
  const gamma = 1.4, c2 = 1 - q * q / 5;
  if (!(c2 > 0)) throw new Error('Ringleb static enthalpy must be positive.');
  const c = Math.sqrt(c2), rho = c ** 5, z = Math.sqrt(1 - q * q * psi * psi);
  const J = 1 / c + 1 / (3 * c ** 3) + 1 / (5 * c ** 5) - Math.atanh(c);
  const x = .5 * ((1 / (q * q) - 2 * psi * psi) / rho + J), y = psi * z / (rho * q);
  const F = q * q * z * z - c2;
  // This patch uses F<0, for which q decreases in the downstream direction.
  if (!(F < 0)) throw new Error('Ringleb patch must remain before the hodograph envelope.');
  const derivatives = {
    x: { q: F / (rho * q ** 3 * c2), psi: -2 * psi / rho },
    y: { q: psi * F / (rho * q * q * z * c2), psi: (1 - 2 * q * q * psi * psi) / (rho * q * z) },
    rho: { q: -rho * q / c2, psi: 0 },
    p: { q: -rho * q, psi: 0 },
    u: { q: (1 - 2 * q * q * psi * psi) / z, psi: -(q ** 3) * psi / z },
    v: { q: 2 * q * psi, psi: q * q }
  };
  const determinant = derivatives.x.q * derivatives.y.psi - derivatives.x.psi * derivatives.y.q;
  if (![x, y, rho, c, determinant, ...Object.values(derivatives).flatMap(d => Object.values(d))].every(Number.isFinite))
    throw new Error('Nonfinite Ringleb state; use a regular finite patch away from q=0.');
  const spatialGradient = d => ({
    x: (d.q * derivatives.y.psi - d.psi * derivatives.y.q) / determinant,
    y: (-d.q * derivatives.x.psi + d.psi * derivatives.x.q) / determinant
  });
  return { x, y, q, psi, rho, p: rho * c2 / gamma, u: q * z, v: q * q * psi,
    c, mach: q / c, enthalpy: c2 / (gamma - 1), h0: 2.5, gamma,
    entropyOverR: Math.log(c2) / (gamma - 1) - Math.log(rho),
    envelopeMargin: -F, determinant, derivatives, spatialGradient };
}

// Five stations give the two-history upstream stencil of the central cell.
// Its reference point r(qCenter, psiCenter) is fixed; finite-width physical
// bank midpoints differ from it and move by O(width^2) under refinement.
// Both parameter widths halve together in a local consistency study.
export function ringlebLocalTube({ qCenter, psiCenter = .7, qHalfStep, psiWidth }) {
  if (!(qHalfStep > 0 && psiWidth > 0)) throw new Error('Positive Ringleb refinement widths required.');
  const qs = [3, 2, 1, 0, -1].map(i => qCenter + i * qHalfStep);
  const lower = qs.map(q => ringlebPoint(q, psiCenter - .5 * psiWidth));
  const upper = qs.map(q => ringlebPoint(q, psiCenter + .5 * psiWidth));
  const parameterSections = qs.slice(1).map((q, i) => ringlebPoint(.5 * (q + qs[i]), psiCenter));
  const exactSections = parameterSections.map((s, i) => {
    const x = (lower[i].x + lower[i + 1].x + upper[i].x + upper[i + 1].x) / 4;
    const y = (lower[i].y + lower[i + 1].y + upper[i].y + upper[i + 1].y) / 4;
    return ringlebAtPoint(x, y, { guess: { q: s.q, psi: s.psi } });
  });
  return { lower, upper, densities: exactSections.map(s => s.rho), massFlow: psiWidth,
    stagnationEnthalpy: 2.5, gamma: 1.4, exactSections, parameterSections,
    qCenter, psiCenter, qHalfStep, psiWidth };
}

// Bounded two-variable inverse hodograph mapping, not an Euler flow solve.
// Initial guesses are supplied from the known analytic patch. Stop at a
// machine-scaled coordinate residual; never alter density to fit cell mass.
export function ringlebAtPoint(x, y, { guess, maxIterations = 8 } = {}) {
  if (![x, y, guess?.q, guess?.psi].every(Number.isFinite)
    || !Number.isInteger(maxIterations) || maxIterations < 0 || maxIterations > 16)
    throw new Error('Supply a finite Ringleb inverse-map target and initial coordinates.');
  let q = guess.q, psi = guess.psi, maxCondition = 0;
  const tolerance = 64 * Number.EPSILON * Math.max(1, Math.abs(x), Math.abs(y));
  for (let iteration = 0; iteration <= maxIterations; iteration++) {
    const s = ringlebPoint(q, psi), d = s.derivatives;
    const normJ = Math.max(Math.abs(d.x.q) + Math.abs(d.x.psi), Math.abs(d.y.q) + Math.abs(d.y.psi));
    const normInverse = Math.max(Math.abs(d.y.psi) + Math.abs(d.x.psi), Math.abs(d.y.q) + Math.abs(d.x.q)) / Math.abs(s.determinant);
    const condition = normJ * normInverse;
    if (!Number.isFinite(condition) || condition > 1e8) throw new Error('Ill-conditioned Ringleb inverse map.');
    maxCondition = Math.max(maxCondition, condition);
    const rx = s.x - x, ry = s.y - y, residual = Math.hypot(rx, ry);
    if (residual <= tolerance) return { ...s, mapping: { iterations: iteration, residual, tolerance, maxCondition } };
    if (iteration === maxIterations) throw new Error(`Ringleb inverse map exceeded ${maxIterations} iterations; residual=${residual}.`);
    q -= (d.y.psi * rx - d.x.psi * ry) / s.determinant;
    psi -= (-d.y.q * rx + d.x.q * ry) / s.determinant;
  }
}
