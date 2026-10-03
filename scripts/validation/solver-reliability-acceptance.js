// SPDX-License-Identifier: GPL-2.0-or-later
// Numerical acceptance is deliberately distinct from physical validation.
import { coupledConvergenceSatisfied, COUPLED_CHANGE_MAXIMUM } from '../../src/euler/streamtube-coupled-convergence.js';
export function assessSolverResult(caseData, result) {
  const quad = caseData.flowModel === 'streamtube-grid';
  const bl = caseData.flowModel === 'coupled' || caseData.quadBoundaryLayers === true;
  const checks = {};
  checks.terminalConvergence = result?.converged !== false
    && !['unconverged', 'failed', 'outside-model'].includes(result?.status)
    && (result?.status === 'solved' || result?.status === 'research-converged' || result?.converged === true);
  const residuals = [result?.diagnostics?.equationResidual, result?.diagnostics?.linearResidual,
    ...Object.values(result?.families ?? {})].filter(v => v !== undefined);
  const residual = residuals.length ? Math.max(...residuals) : undefined;
  const tolerance = quad ? 1e-10 : bl ? 1e-8 : 1e-9;
  const residualConverged = residuals.length > 0 && residuals.every(v => Number.isFinite(v) && v <= tolerance && v >= 0);
  checks.equations = residualConverged || quad && bl && result?.solverSettings?.convergence === 'mses'
    && residuals.every(v => Number.isFinite(v) && v >= 0 && v <= COUPLED_CHANGE_MAXIMUM)
    && coupledConvergenceSatisfied(result, tolerance);
  checks.requestedMach = Number.isFinite(result?.mach) && Math.abs(result.mach - caseData.mach) <= 1e-12
    && result?.targetReached !== false && result?.requestedConditionReached !== false
    && result?.machContinuation?.reachedTarget !== false && result?.continuation?.reachedTarget !== false;
  checks.requestedIncidence = Number.isFinite(result?.alpha) && Math.abs(result.alpha - caseData.alpha) <= 1e-12;
  checks.coefficientsFinite = ['cl', 'cm', ...(bl ? ['cd'] : [])].every(k => Number.isFinite(result?.[k]));
  if (quad) checks.convexGrid = result?.mesh?.quality?.valid === true;
  if (quad && bl) {
    if (caseData.reynolds !== undefined) checks.requestedReynolds = result?.referenceReynolds === caseData.reynolds;
    if (caseData.ncrit !== undefined) checks.requestedNcrit = result?.conditions?.ncrit === caseData.ncrit
      && result?.ncritContinuation?.reachedTarget !== false;
    if (caseData.materialTrips !== undefined) checks.requestedTrips = JSON.stringify(result?.materialTrips)
      === JSON.stringify(caseData.materialTrips);
  }
  if (bl) {
    const layers = result?.boundaryLayer;
    checks.completeBoundaryLayers = Array.isArray(layers?.surfaces) && layers.surfaces.length === 2 * caseData.elements.length
      && Array.isArray(layers?.wakes) && layers.wakes.length === caseData.elements.length
      && [...layers.surfaces, ...layers.wakes].every(s => Array.isArray(s.stations) && s.stations.length > 1
        && s.stations.every(p => Number.isFinite(p.theta) && p.theta > 0 && Number.isFinite(p.ue) && p.ue > 0
          && Number.isFinite(p.deltaStar) && Number.isFinite(p.wakeGap ?? 0) && (p.wakeGap ?? 0) >= 0
          && p.deltaStar - (p.wakeGap ?? 0) > 0))
      && caseData.elements.every((_, element) => ['upper', 'lower'].every(side => layers.surfaces
        .filter(s => s.element === element && s.side === side).length === 1)
        && layers.wakes.filter(w => w.element === element).length === 1);
    if (!quad) checks.wakeConvergence = Number.isFinite(result?.diagnostics?.wakeResidual)
      && result.diagnostics.wakeResidual <= (result.diagnostics.wakeTolerance ?? 1e-6);
  }
  return { passed: Object.values(checks).every(Boolean), checks, residual, tolerance, residualConverged,
    failures: Object.keys(checks).filter(k => !checks[k]), physicalValidation: 'not evaluated by this numerical gate' };
}
