// SPDX-License-Identifier: GPL-2.0-or-later
import { normInf, solveLinear } from './linear.js';

// All unknowns are updated together. Residuals and variables must be scaled by
// the caller. An analytic/AD Jacobian may replace this verification backend.
export function solveNewton({ residual, initial, jacobian, admissible = () => true,linearSolve=solveLinear,
  tolerance = 1e-10, maxIterations = 40, onIteration, onState, activeSet }) {
  if (activeSet && !['snapshot', 'restore', 'prepare'].every(k => typeof activeSet[k] === 'function'))
    throw new Error('Invalid Newton active-set callbacks.');
  let x = Float64Array.from(initial);
  const n = x.length;
  if (!n || !x.every(Number.isFinite) || !admissible(x)) throw new Error('Invalid Newton initial state.');
  if (!(tolerance > 0) || !Number.isFinite(tolerance) || !Number.isInteger(maxIterations) || maxIterations < 0) {
    throw new Error('Invalid Newton convergence controls.');
  }
  const evaluate = state => {
    const r = Float64Array.from(residual(state));
    if (r.length !== n || !r.every(Number.isFinite)) throw new Error('Invalid nonlinear residual.');
    return r;
  };
  let r = evaluate(x);
  const history = [{ iteration: 0, residual: normInf(r), step: 0 }];
  onIteration?.(history[0]);
  // Isolated snapshots of accepted states only; observers cannot mutate the
  // solve, and rejected line-search proposals are never published.
  onState?.({ x: x.slice(), iteration: { ...history[0] } });
  for (let iteration = 0; iteration < maxIterations; iteration++) {
    if (normInf(r) <= tolerance) return { converged: true, x, history, reason: 'residual' };
    let j;
    if (jacobian) j = jacobian(x);
    else {
      j = new Float64Array(n * n);
      for (let col = 0; col < n; col++) {
        const h = Math.cbrt(Number.EPSILON) * Math.max(1, Math.abs(x[col]));
        const plus = x.slice(); const minus = x.slice();
        plus[col] += h; minus[col] -= h;
        const plusOK = admissible(plus); const minusOK = admissible(minus);
        if (!plusOK && !minusOK) throw new Error('No admissible Jacobian perturbation.');
        const rp = plusOK ? evaluate(plus) : r;
        const rm = minusOK ? evaluate(minus) : r;
        const width = h * (plusOK && minusOK ? 2 : 1);
        for (let row = 0; row < n; row++) j[row * n + col] = (rp[row] - rm[row]) / width;
      }
    }
    let delta;
    try { delta = linearSolve(j, r.map(v => -v)); }
    catch (error) { return { converged: false, x, history, reason: error.message }; }
    let accepted = false;
    const saved = activeSet?.snapshot();
    try {
      for (let step = 1; step >= 2 ** -20; step *= 0.5) {
        if (activeSet) activeSet.restore(saved);
        const candidate = x.map((value, i) => value + step * delta[i]);
        let next, event;
        try {
          event = activeSet?.prepare(candidate, x);
          if (!candidate.every(Number.isFinite) || !admissible(candidate)) continue;
          next = evaluate(candidate);
        } catch { continue; }
        // A bounded phase event can change an auxiliary variable's physical
        // meaning. As in dogleg, accept its admissible transfer and rebuild
        // the Jacobian; an old-phase merit comparison would be misleading.
        // The caller's prepare callback must bound physical event travel.
        if (event?.changed || normInf(next) <= (1 - 1e-4 * step) * normInf(r)) {
          x = candidate; r = next; accepted = true;
          history.push({ iteration: iteration + 1, residual: normInf(r), step,
            ...(event?.changed ? { activeChange: true, changes: event.changes, meritComparable: false } : {}) });
          onIteration?.(history.at(-1));
          onState?.({ x: x.slice(), iteration: { ...history.at(-1) } });
          break;
        }
      }
    } finally { if (activeSet && !accepted) activeSet.restore(saved); }
    if (!accepted) return { converged: false, x, history, reason: 'line search failed' };
  }
  return { converged: normInf(r) <= tolerance, x, history,
    reason: normInf(r) <= tolerance ? 'residual' : 'iteration limit' };
}
