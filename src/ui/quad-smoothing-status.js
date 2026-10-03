// SPDX-License-Identifier: GPL-2.0-or-later
// Initial-grid smoothing and final flow convergence are separate outcomes.
export function quadSmoothingLabel(smoothing, { summary = false } = {}) {
  if (smoothing?.initialGuessAccepted === true && smoothing.converged === false)
    return summary ? 'incomplete · admissible initial grid retained' : 'SLOR incomplete · admissible initial grid retained';
  if (smoothing?.retainedOriginal)
    return summary ? 'rejected · original mesh retained' : 'SLOR rejected; original mesh retained';
  if (smoothing?.converged)
    return summary ? 'completed' : `SLOR smoothed${smoothing.fixedSpacingPassages?.length ? ' · fixed spacing used' : ''}${smoothing.harmonicPassages?.length ? ' · harmonic SLOR used' : ''}`;
  if (smoothing?.attempted) return summary ? 'in progress' : 'SLOR smoothing in progress';
  return summary ? 'not run' : 'SLOR not run';
}
