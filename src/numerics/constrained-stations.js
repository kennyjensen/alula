// SPDX-License-Identifier: GPL-2.0-or-later
// Project requested interval lengths onto a fixed-length, fixed-endpoint
// distribution with bounded adjacent growth. This is a convex spacing
// reconstruction, not a recovered MSET routine or a grid-flow smoother.
// Dykstra projections solve the diagonal least-squares problem; no global
// flow or mesh iteration is involved.
import { solveChainQuadratic } from './chain-quadratic-program.js';
export function stationGrowthFeasibility({ intervals, firstSpacing, lastSpacing, maximumGrowth = 1.5 }) {
  if (!Number.isInteger(intervals) || intervals < 3 || ![firstSpacing, lastSpacing].every(v => Number.isFinite(v) && v > 0)
    || !Number.isFinite(maximumGrowth) || maximumGrowth <= 1) throw new Error('Invalid station growth feasibility inputs.');
  const logR = Math.log(maximumGrowth), logA = Math.log(firstSpacing), logB = Math.log(lastSpacing);
  const lower = Array.from({ length: intervals }, (_, i) => Math.exp(Math.max(logA - i * logR, logB - (intervals - 1 - i) * logR)));
  const upper = Array.from({ length: intervals }, (_, i) => Math.exp(Math.min(0, logA + i * logR, logB + (intervals - 1 - i) * logR)));
  const lowerSum = lower.reduce((s, h) => s + h, 0), upperSum = upper.reduce((s, h) => s + h, 0);
  const feasible = !lower.some((h, i) => h > upper[i] * (1 + 64 * Number.EPSILON))
    && lowerSum <= 1 + 64 * Number.EPSILON && upperSum >= 1 - 64 * Number.EPSILON;
  return { feasible, lowerSum, upperSum };
}

export function fitConstrainedStations({ positions, firstSpacing, lastSpacing, maximumGrowth = 1.5, tolerance = 1e-10, maxSweeps = 20000, method = 'dykstra' }) {
  if (!Array.isArray(positions) || positions.length < 4 || positions[0] !== 0 || positions.at(-1) !== 1
    || positions.some((s, i) => !Number.isFinite(s) || i && !(s > positions[i - 1]))
    || ![firstSpacing, lastSpacing].every(h => Number.isFinite(h) && h > 0) || !(firstSpacing + lastSpacing < 1)
    || !Number.isFinite(maximumGrowth) || maximumGrowth <= 1 || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isInteger(maxSweeps) || maxSweeps < 1 || !['dykstra', 'primal-dual'].includes(method))
    throw new Error('Constrained stations require normalized ordered nodes, feasible endpoint lengths and growth above one.');
  const count = positions.length - 1, natural = positions.slice(1).map((s, i) => s - positions[i]);
  const { feasible, lowerSum, upperSum } = stationGrowthFeasibility({ intervals: count, firstSpacing, lastSpacing, maximumGrowth });
  // These bounds are sufficient as well as necessary: each envelope obeys
  // the growth inequalities, so their convex combination spans every
  // admissible total length between lowerSum and upperSum.
  if (!feasible)
    throw new Error(`Endpoint/growth constraints are infeasible for ${count} intervals (total bounds ${lowerSum}, ${upperSum}).`);
  const n = count - 2, weights = natural.slice(1, -1), x = Array(n).fill(1), constraints = [];
  const add = (terms, rhs, equality = false) => {
    const norm = Math.hypot(...terms.map(t => t[1]));
    constraints.push({ terms: terms.map(([i, a]) => [i, a / norm]), rhs: rhs / norm, equality, correction: 0 });
  };
  add(weights.map((h, i) => [i, h]), 1 - firstSpacing - lastSpacing, true);
  for (let i = 1; i < count; i++) for (const direction of [1, -1]) {
    const a = direction === 1 ? i : i - 1, b = direction === 1 ? i - 1 : i;
    let rhs = 0; const terms = [];
    for (const [k, factor] of [[a, 1], [b, -maximumGrowth]]) {
      if (k === 0 || k === count - 1) rhs -= factor * (k === 0 ? firstSpacing : lastSpacing);
      else terms.push([k - 1, factor * weights[k - 1]]);
    }
    add(terms, rhs);
  }
  let sweeps = 0, change = Infinity, primalResidual = Infinity;
  if (method === 'primal-dual') {
    const r = solveChainQuadratic({ constraints, tolerance, maxIterations: Math.min(maxSweeps, 100) });
    if (!r.converged) throw new Error(`Constrained station projection did not converge in ${r.iterations} primal-dual iterations (primal ${r.primalResidual}, dual ${r.dualResidual}, complementarity ${r.complementarity}).`);
    r.x.forEach((v, i) => { x[i] = v; });
    constraints[0].correction = r.equalityMultiplier;
    r.inequalityMultipliers.forEach((v, i) => { constraints[i + 1].correction = v; });
    sweeps = r.iterations; change = 0; primalResidual = r.primalResidual;
  } else for (; sweeps < maxSweeps; sweeps++) {
    const before = x.slice();
    for (const c of constraints) {
      for (const [i, a] of c.terms) x[i] += c.correction * a;
      const residual = c.terms.reduce((s, [i, a]) => s + a * x[i], 0) - c.rhs;
      const correction = c.equality ? residual : Math.max(0, residual);
      for (const [i, a] of c.terms) x[i] -= correction * a;
      c.correction = correction;
    }
    change = Math.max(...x.map((v, i) => Math.abs(v - before[i])));
    primalResidual = Math.max(...constraints.map(c => {
      const r = c.terms.reduce((s, [i, a]) => s + a * x[i], 0) - c.rhs;
      return c.equality ? Math.abs(r) : Math.max(0, r);
    }));
    if (change <= tolerance && primalResidual <= tolerance) { sweeps++; break; }
  }
  if (change > tolerance || primalResidual > tolerance) throw new Error(`Constrained station projection did not converge in ${maxSweeps} sweeps.`);
  const intervals = [firstSpacing, ...x.map((v, i) => v * weights[i]), lastSpacing];
  const output = [0]; intervals.forEach(h => output.push(output.at(-1) + h));
  const totalResidual = output.at(-1) - 1;
  if (intervals.some(h => !(h > 0) || !Number.isFinite(h)) || Math.abs(totalResidual) > 10 * tolerance)
    throw new Error('Constrained stations failed their positive-length or total check.');
  output[output.length - 1] = 1;
  const stationarity = x.map(v => v - 1);
  for (const c of constraints) for (const [i, a] of c.terms) stationarity[i] += c.correction * a;
  const complementarity = Math.max(0, ...constraints.filter(c => !c.equality).map(c => Math.abs(c.correction
    * (c.terms.reduce((s, [i, a]) => s + a * x[i], 0) - c.rhs))));
  return { positions: output, sweeps, primalResidual, totalResidual, complementarity,
    ...(method === 'primal-dual' ? { method } : {}),
    stationarityResidual: Math.max(...stationarity.map(Math.abs)), lowerSum, upperSum,
    objective: x.reduce((s, v) => s + .5 * (v - 1) ** 2, 0), maximumGrowth,
    policy: 'minimum squared relative changes of interior intervals subject to total length, endpoint lengths and adjacent growth' };
}
