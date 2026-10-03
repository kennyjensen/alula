// SPDX-License-Identifier: GPL-2.0-or-later
// Temporary, user-authorized comparison policy; no numerical solver changes.
export const originalGates = Object.freeze({ cpMax: .04, ueMax: .02, thetaRelativeMax: .05, deltaStarRelativeMax: .08 });
export const currentEffective = Object.freeze({ ...originalGates, thetaRelativeMax: .12, deltaStarRelativeMax: .15 });
export const residualLimit = 1e-10;
export const coupledAcceptancePolicy = Object.freeze({
  id: 'temporary-thickness-screen-2026-09-12', temporary: true,
  authorization: 'User requested temporarily relaxing thickness gates so the current converged baseline passes, while prioritizing coupled multielement convergence.',
  scope: 'Provisional single-grid native profile screen. No claim of force/refinement acceptance, multielement validation, or complete solver capability.',
  originalGates, currentEffective, residualLimit, requireFullCoverage: true
});

// Every readiness flag is mandatory. In particular, a small residual alone
// cannot promote a partial, inadmissible, stale, or unresolved saved result.
export function evaluateCoupledAcceptance(evidence = {}) {
  const failures = [], checks = {};
  for (const key of ['complete', 'converged', 'admissible', 'provenanceValid', 'exactReplay']) {
    checks[key] = evidence[key] === true;
    if (!checks[key]) failures.push(`${key} must be true.`);
  }
  checks.residual = Number.isFinite(evidence.residual) && evidence.residual >= 0 && evidence.residual < residualLimit;
  if (!checks.residual) failures.push(`Residual ${evidence.residual} must be finite, nonnegative, and below ${residualLimit}.`);
  const coverage = evidence.coverage ?? {};
  for (const field of ['BL', 'Pressure']) {
    const requested = coverage[`requested${field}`], covered = coverage[`covered${field}`];
    checks[`full${field}Coverage`] = Number.isSafeInteger(requested) && requested > 0 && covered === requested;
    if (!checks[`full${field}Coverage`]) failures.push(`Complete ${field} coverage is required (${covered}/${requested}).`);
  }
  checks.noUncoveredStations = coverage.uncoveredStations === 0;
  if (!checks.noUncoveredStations) failures.push('Uncovered native stations are not allowed.');
  const errors = { ...evidence.errors }, gateResults = {};
  for (const [key, bound] of Object.entries(currentEffective)) {
    const value = errors[key], valid = Number.isFinite(value) && value >= 0;
    const passed = valid && value <= bound;
    gateResults[key] = { value, original: originalGates[key], currentEffective: bound,
      originalPassed: valid && value <= originalGates[key], passed };
    if (!passed) failures.push(`${key}: ${value} must be finite, nonnegative, and at most ${bound}.`);
  }
  return { accepted: failures.length === 0, policy: coupledAcceptancePolicy,
    originalGates: { ...originalGates }, currentEffective: { ...currentEffective },
    errors, residual: evidence.residual, coverage: { ...coverage }, checks, gateResults, failures,
    physicalAcceptance: false, fullSolverComplete: false, runtimeChanged: false };
}
