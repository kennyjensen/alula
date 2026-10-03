// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded monitoring of accepted physical updates. Thresholds are recovery
// policies, not convergence tolerances. Paired merits use identical equations.
export function repeatedGridStepLimit(history, tolerance) {
  if (history.length < 4) return null;
  const before = history.at(-4), recent = history.slice(-3), last = recent.at(-1);
  if (!(last.residual > tolerance) || recent.some(h => !(h.step > 0 && h.step < 1e-5))
    || last.residual < .999 * before.residual) return null;
  const law = h => JSON.stringify([h.dissipation?.mucon, h.dissipation?.mcrit]);
  if (recent.some(h => law(h) !== law(before))) return null;
  const cells = recent.map(h => h.rejections?.find(r => r.code === 'streamtube-grid-step')?.diagnostics?.cell);
  if (!cells[0] || cells.some(c => JSON.stringify(c) !== JSON.stringify(cells[0]))) return null;
  return { cell: cells[0], iterations: recent.map(h => h.iteration), steps: recent.map(h => h.step),
    beforeResidual: before.residual, residual: last.residual };
}

// A passage chart can stall well before a corner becomes nearly flat.
// Compare each paired merit under its own law; a single pressure peak or a
// substantial dissipation-law change is not evidence that the Newton solve has stalled.
export function repeatedPassageStepLimit(history, tolerance) {
  const recent = history.slice(-6);
  if (recent.length < 6 || !(recent.at(-1).residual > tolerance)
    || recent.some(h => !(h.step > 0 && h.step < .02))) return null;
  const first = recent[0].dissipation;
  if (recent.some(h => h.dissipation?.mucon !== first?.mucon
    || Math.abs((h.dissipation?.mcrit ?? 0) - (first?.mcrit ?? 0)) > 1e-5)) return null;
  const ratios = recent.map(h => h.residualDecrease?.afterSquaredNorm / h.residualDecrease?.beforeSquaredNorm);
  if (ratios.some(r => !Number.isFinite(r) || r <= 0) || ratios.reduce((a, b) => a * b, 1) < .9) return null;
  return { iterations: recent.map(h => h.iteration), steps: recent.map(h => h.step),
    pairedSquaredMeritRatio: ratios.reduce((a, b) => a * b, 1), residual: recent.at(-1).residual };
}

export const residualMerit = residual => {
  let scale = 0;
  for (const r of residual) scale = Math.max(scale, Math.abs(r));
  if (scale === 0) return 0;
  let sum = 0;
  for (const r of residual) sum += (r / scale) ** 2;
  return scale * Math.sqrt(sum / residual.length);
};
export function relativeStateChange(a, b) {
  if (!a || !b || a.length !== b.length) return Infinity;
  let change = 0;
  for (let i = 0; i < a.length; i++) change = Math.max(change, Math.abs(a[i] - b[i]) / Math.max(1, Math.abs(a[i]), Math.abs(b[i])));
  return change;
}
export function physicalIterationVector(state, nodes) {
  const vector = Array.from(state);
  for (const group of nodes) for (const row of group) for (const p of row) vector.push(p.x, p.y);
  return vector;
}
export function createIterationProgress() {
  let states = [], ratios = [], tiny = 0, cycles = 0;
  return {
    seed(vector, phase, mcrit) { states = [{ vector, phase, mcrit }]; ratios = []; tiny = 0; cycles = 0; },
    observe({ vector, before, after, residual, tolerance, phase, mcrit, phaseChanged = false, maintenance = false }) {
      const previous = states.at(-1), twoAgo = states.at(-2);
      const movement = relativeStateChange(vector, previous?.vector), returnDistance = relativeStateChange(vector, twoAgo?.vector);
      // A transition event changes the row interpretation: never use its
      // merit ratio to declare progress or stagnation.
      const ratio = !phaseChanged && before > 0 ? after / before : null;
      ratios.push(ratio); if (ratios.length > 6) ratios.shift();
      tiny = movement <= 32 * Number.EPSILON && ratio !== null && ratio >= .99 ? tiny + 1 : 0;
      const repeated = phase === twoAgo?.phase && movement > 1e-9 && returnDistance <= Math.max(1e-12, .05 * movement)
        && twoAgo?.residual !== undefined && residual >= .9 * twoAgo.residual;
      cycles = repeated ? cycles + 1 : 0;
      const plateau = ratios.length === 6 && ratios.every(r => r !== null && Number.isFinite(r))
        && ratios.reduce((p, r) => p * r, 1) >= .9;
      const cause = residual <= tolerance ? null : tiny >= 3 ? 'negligible-state-change'
        : cycles >= 2 ? 'two-state-cycle' : plateau ? 'residual-plateau' : null;
      const diagnostics = { cause, stateChange: movement, twoStepChange: Number.isFinite(returnDistance) ? returnDistance : null,
        pairedMeritRatio: ratio, transitionChanged: phaseChanged || previous?.phase !== phase,
        mcritChanged: previous?.mcrit !== mcrit, gridMaintenance: maintenance };
      states.push({ vector, phase, mcrit, residual }); if (states.length > 2) states.shift();
      return diagnostics;
    },
    canExtend() { return ratios.length >= 3 && ratios.slice(-3).every(r => r !== null && r < .8)
      && ratios.slice(-3).reduce((p, r) => p * r, 1) < .2; },
  };
}

// Optional line-search acceptance on the existing scaled governing rows.
// This function does not evaluate a flow, change a state or relax a domain.
export function requireCoupledResidualDecrease(current, candidate, { step, tolerance, weights } = {}) {
  const vector = r => (Array.isArray(r) || ArrayBuffer.isView(r)) && r.length > 0
    && Array.from(r).every(Number.isFinite);
  if (!vector(current) || !vector(candidate) || current.length !== candidate.length
    || weights !== undefined && (!vector(weights) || weights.length !== current.length || weights.some(w => w <= 0))
    || !Number.isFinite(step) || step <= 0 || step > 1 || !Number.isFinite(tolerance) || tolerance <= 0)
    throw new Error('Invalid coupled residual-decrease controls or governing rows.');
  const squared = r => r.reduce((sum, value, i) => sum + (value * (weights?.[i] ?? 1)) ** 2, 0);
  const beforeSquaredNorm = squared(current), afterSquaredNorm = squared(candidate);
  if (!Number.isFinite(beforeSquaredNorm) || !Number.isFinite(afterSquaredNorm))
    throw new Error('Nonfinite coupled squared residual norm.');
  const coefficient = 1e-4, allowedSquaredNorm = (1 - coefficient * step) * beforeSquaredNorm;
  const maximumResidual = candidate.reduce((value, r) => Math.max(value, Math.abs(r)), 0);
  const converged = maximumResidual <= tolerance;
  const diagnostics = { method: 'armijo-squared-residual', coefficient, step, beforeSquaredNorm,
    afterSquaredNorm, allowedSquaredNorm, maximumResidual, tolerance, converged };
  if (!converged && afterSquaredNorm > allowedSquaredNorm)
    throw Object.assign(new Error(`Maintained coupled squared residual ${afterSquaredNorm} exceeds the allowed ${allowedSquaredNorm}.`), {
      code: 'COUPLED_RESIDUAL_DECREASE', diagnostics,
    });
  return diagnostics;
}
