// SPDX-License-Identifier: GPL-2.0-or-later
// A failed elliptic initializer may retain a well-conditioned initial guess.
// This is deliberately separate from elliptic and Euler/BL convergence.
import { streamtubeGridConvexity } from './streamtube-convex-step.js';
export const PARTIAL_SLOR_MINIMUM_CORNER_SINE = .1;
const allowed = new Set(['no-admissible-convex-correction', 'no-admissible-decreasing-line-step', 'sweep-limit', 'roundoff-stagnation']);
export const isRecoverableSlorTermination = termination => termination?.origin === 'solver' && allowed.has(termination.termination);
const finiteMeasure = h => Number.isFinite(h?.residual) && h.residual >= 0 && Number.isFinite(h?.merit) && h.merit >= 0;

export function createPartialSlorSeedRecorder({ nodes, stage }) {
  if (typeof stage !== 'string' || !stage) throw new Error('A fixed named smoothing stage is required.');
  const seed = structuredClone(nodes), seedQuality = streamtubeGridConvexity([seed]);
  // Reuse the existing grid-spacing audit's 0.1 obliqueness screen only as
  // failed-process recovery eligibility. The physical 1e-12 domain and the
  // requested residual tolerance are unchanged. Successful SLOR ignores it.
  const eligible = seedQuality.valid && seedQuality.minCornerSine >= PARTIAL_SLOR_MINIMUM_CORNER_SINE;
  let first = null, best = null, lastIteration = -1;
  const observe = (history, current) => {
    if (!eligible || history?.stage !== stage || !Number.isInteger(history.iteration)
      || history.iteration <= lastIteration || !finiteMeasure(history)) return false;
    lastIteration = history.iteration;
    if (!Array.isArray(current) || current.length !== seed.length
      || current.some((row, i) => !Array.isArray(row) || row.length !== seed[i].length)) return false;
    for (let i = 0; i < seed.length; i++) for (let j = 0; j < seed[i].length; j++)
      if (!i || i === seed.length - 1 || !j || j === seed[i].length - 1)
        if (current[i][j].x !== seed[i][j].x || current[i][j].y !== seed[i][j].y) return false;
    const quality = streamtubeGridConvexity([current]);
    if (!quality.valid || quality.minCornerSine < seedQuality.minCornerSine) return false;
    if (!first) {
      if (history.iteration !== 0) return false;
      for (let i = 0; i < seed.length; i++) for (let j = 0; j < seed[i].length; j++)
        if (current[i][j].x !== seed[i][j].x || current[i][j].y !== seed[i][j].y) return false;
      first = { residual: history.residual, merit: history.merit }; return true;
    }
    if (!(history.residual < first.residual && history.merit < first.merit)
      || best && !(history.residual < best.residual)) return false;
    best = { nodes: structuredClone(current), residual: history.residual, merit: history.merit, iteration: history.iteration, quality };
    return true;
  };
  const afterRequestedFailure = failure => {
    if (!eligible || !first || failure?.completedRequestedStages !== true
      || failure.converged !== false || !isRecoverableSlorTermination(failure)) return null;
    return { nodes: structuredClone(best?.nodes ?? seed), report: { smoothingConverged: false,
      selected: best ? 'lower-residual accepted iterate with nondegraded seed corner conditioning' : 'original admissible seed',
      stage, iteration: best?.iteration ?? 0, residual: best?.residual ?? first.residual,
      merit: best?.merit ?? first.merit, seedMeasure: structuredClone(first), seedQuality: structuredClone(seedQuality),
      selectedQuality: structuredClone(best?.quality ?? seedQuality), eligibilityCornerSine: PARTIAL_SLOR_MINIMUM_CORNER_SINE,
      eligibilityMeaning: 'Existing spacing-audit conditioning screen; not a physical-validity or convergence threshold.',
      requestedFailure: structuredClone(failure),
      scope: 'Initial guess only; requested elliptic tolerance remains unmet and Euler/BL convergence remains separate.' } };
  };
  return { observe, afterRequestedFailure };
}
