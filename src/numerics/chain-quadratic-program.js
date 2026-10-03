// SPDX-License-Identifier: GPL-2.0-or-later
// Mehrotra predictor-corrector for min 1/2 ||x-1||^2 with one equality
// and inequalities coupling adjacent variables. The eliminated Hessian is
// SPD tridiagonal; a scalar Schur complement handles the total-length row.
import { polishChainQuadratic } from './chain-quadratic-polish.js';
export function solveChainQuadratic({ constraints, tolerance = 1e-10, maxIterations = 100 }) {
  const equality = constraints?.[0], rows = constraints?.slice(1), n = equality?.terms?.length;
  if (!n || !equality.equality || !rows.length || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isInteger(maxIterations) || maxIterations < 1
    || constraints.some(c => !Number.isFinite(c.rhs) || !Array.isArray(c.terms) || c.terms.some(([i, a]) => !Number.isInteger(i) || i < 0 || i >= n || !Number.isFinite(a)))
    || rows.some(c => c.equality || c.terms.length > 2 || c.terms.length < 1
      || c.terms.length === 2 && Math.abs(c.terms[0][0] - c.terms[1][0]) !== 1)) throw new Error('Invalid chain quadratic program.');
  const dot = (c, x) => c.terms.reduce((v, [i, a]) => v + a * x[i], 0);
  const e = Array(n).fill(0); equality.terms.forEach(([i, a]) => { e[i] = a; });
  const x = Array(n).fill(1), slack = rows.map(c => Math.max(1, c.rhs - dot(c, x))), dual = rows.map(() => 1);
  let lambda = 0, iterations = 0, primalResidual = Infinity, dualResidual = Infinity, complementarity = Infinity;
  const stepLimit = (values, step) => Math.min(1, ...values.map((v, k) => step[k] < 0 ? -v / step[k] : Infinity));
  for (; iterations < maxIterations; iterations++) {
    const rd = x.map((v, i) => v - 1 + e[i] * lambda), rp = rows.map((c, k) => dot(c, x) + slack[k] - c.rhs), re = dot(equality, x) - equality.rhs;
    rows.forEach((c, k) => c.terms.forEach(([i, a]) => { rd[i] += a * dual[k]; }));
    primalResidual = Math.max(Math.abs(re), ...rp.map(Math.abs)); dualResidual = Math.max(...rd.map(Math.abs));
    complementarity = Math.max(...slack.map((s, k) => s * dual[k]));
    if (Math.max(primalResidual, dualResidual, complementarity) <= tolerance) break;
    // A near-optimal barrier Hessian can lose its unit diagonal to very
    // large active-row weights. Solve a candidate active chain directly;
    // its full original KKT certificate must pass the same tolerance.
    if (Math.max(primalResidual, dualResidual, complementarity) <= Math.sqrt(tolerance)) {
      const polished = polishChainQuadratic({ constraints, slack, dual, equalityMultiplier: lambda, tolerance });
      if (polished) return { ...polished, iterations, converged: true,
        method: 'primal-dual predictor-corrector with certified active-chain KKT polishing' };
    }
    const diagonal = Array(n).fill(1), off = Array(Math.max(0, n - 1)).fill(0);
    rows.forEach((c, k) => {
      const weight = dual[k] / slack[k];
      c.terms.forEach(([i, a]) => { diagonal[i] += weight * a * a; });
      if (c.terms.length === 2) { const [[i, a], [j, b]] = c.terms; off[Math.min(i, j)] += weight * a * b; }
    });
    const pivots = [...diagonal], multipliers = [];
    for (let i = 0; i < n; i++) {
      if (i) { multipliers[i - 1] = off[i - 1] / pivots[i - 1]; pivots[i] -= multipliers[i - 1] * off[i - 1]; }
      if (!(pivots[i] > 0) || !Number.isFinite(pivots[i])) throw new Error('Unresolved chain quadratic Hessian.');
    }
    const solve = rhs => {
      const v = [...rhs]; for (let i = 1; i < n; i++) v[i] -= multipliers[i - 1] * v[i - 1];
      for (let i = n - 1; i >= 0; i--) v[i] = (v[i] - (i < n - 1 ? off[i] * v[i + 1] : 0)) / pivots[i];
      return v;
    };
    const he = solve(e), schur = e.reduce((s, v, i) => s + v * he[i], 0);
    if (!(schur > 0) || !Number.isFinite(schur)) throw new Error('Unresolved chain total-length Schur complement.');
    const direction = center => {
      const rhs = rd.map(v => -v);
      rows.forEach((c, k) => c.terms.forEach(([i, a]) => { rhs[i] += a * (center[k] / slack[k] - dual[k] / slack[k] * rp[k]); }));
      const hb = solve(rhs), dl = (e.reduce((s, v, i) => s + v * hb[i], 0) + re) / schur;
      const dx = hb.map((v, i) => v - he[i] * dl), ds = rows.map((c, k) => -rp[k] - dot(c, dx));
      const dz = ds.map((v, k) => (-center[k] - dual[k] * v) / slack[k]);
      return { dx, ds, dz, dl };
    };
    const products = slack.map((s, k) => s * dual[k]), mu = products.reduce((a, b) => a + b, 0) / rows.length;
    const affine = direction(products), ap = stepLimit(slack, affine.ds), ad = stepLimit(dual, affine.dz);
    const muAffine = slack.reduce((sum, s, k) => sum + (s + ap * affine.ds[k]) * (dual[k] + ad * affine.dz[k]), 0) / rows.length;
    const sigma = Math.max(0, Math.min(1, (muAffine / mu) ** 3));
    const step = direction(products.map((v, k) => v + affine.ds[k] * affine.dz[k] - sigma * mu));
    const primalStep = Math.min(1, .995 * stepLimit(slack, step.ds)), dualStep = Math.min(1, .995 * stepLimit(dual, step.dz));
    x.forEach((_, i) => { x[i] += primalStep * step.dx[i]; });
    slack.forEach((_, k) => { slack[k] += primalStep * step.ds[k]; dual[k] += dualStep * step.dz[k]; });
    lambda += dualStep * step.dl;
    if (![...x, ...slack, ...dual, lambda].every(Number.isFinite) || slack.some(s => !(s > 0)) || dual.some(z => !(z > 0)))
      throw new Error('Nonfinite or nonpositive chain quadratic iterate.');
  }
  return { x, equalityMultiplier: lambda, inequalityMultipliers: dual, iterations,
    converged: Math.max(primalResidual, dualResidual, complementarity) <= tolerance,
    primalResidual, dualResidual, complementarity, method: 'primal-dual predictor-corrector with tridiagonal Hessian and scalar equality Schur complement' };
}
