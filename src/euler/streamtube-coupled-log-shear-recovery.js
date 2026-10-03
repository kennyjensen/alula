// SPDX-License-Identifier: GPL-2.0-or-later
// Work selection and an explicit numerical-coordinate transfer only.
// No equation, physical state, mesh, or existing checkpoint is modified.
import { coupledCheckpointShearCoordinate } from './streamtube-coupled-shear-policy.js';

const families = ['euler', 'boundaryLayer', 'edgeMatching'];
const vector = x => (Array.isArray(x) || ArrayBuffer.isView(x)) && x.length > 0 && x.every(Number.isFinite);
const sameVector = (a, b) => vector(a) && vector(b) && a.length === b.length && a.every((v, i) => Object.is(v, b[i]));
const finiteFamilies = f => f && families.every(k => Number.isFinite(f[k]) && f[k] >= 0);
const sameFamilies = (a, b) => finiteFamilies(a) && finiteFamilies(b) && families.every(k => Object.is(a[k], b[k]));
const points = nodes => Array.isArray(nodes) && nodes.length > 0 && nodes.every(group => Array.isArray(group)
  && group.length > 1 && group.every(row => Array.isArray(row) && row.length > 1
    && row.every(p => p && Number.isFinite(p.x) && Number.isFinite(p.y))));
const sameNodes = (a, b) => points(a) && points(b) && a.length === b.length && a.every((g, i) => g.length === b[i].length
  && g.every((row, j) => row.length === b[i][j].length && row.every((p, k) => Object.is(p.x, b[i][j][k].x) && Object.is(p.y, b[i][j][k].y))));
const terminalTrips = x => Array.isArray(x) && x.length > 0 && x.every(pair => Array.isArray(pair)
  && pair.length === 2 && pair.every(v => v === 1));

export function planCoupledLogarithmicShearRecovery(result, { startupPlan, startupAttempt, maxIterations, tolerance } = {}) {
  const cp = result?.checkpoint, f = cp?.restart, c = cp?.continuation, history = result?.history;
  if (startupAttempt !== 2 || startupPlan?.kind !== 'ncrit-transition-startup'
    || startupPlan.sourceNcrit !== 4 || !Number.isFinite(startupPlan.targetNcrit) || !(startupPlan.targetNcrit > startupPlan.sourceNcrit)
    || startupPlan.shearCoordinate !== 'linear' || startupPlan.shearRecovery
    || !Number.isInteger(maxIterations) || maxIterations < 1 || !Number.isFinite(tolerance) || tolerance <= 0
    || result?.converged !== false || ((result.solverStopReason ?? result.reason) !== 'iteration limit'
      && !(result.progressControl?.stopReason && result.progressControl.stopReason === result.reason))
    || result.lastRejectedStep != null || result.shearRecovery || result.automaticRefinement || result.gridSequence
    || result.initialRedistribution?.accepted !== true || result.mesh?.quality?.valid !== true
    || !Array.isArray(result.mesh.quality.invalidCells) || result.mesh.quality.invalidCells.length
    || cp?.version !== 1 || !f?.input || !f.options || !c || typeof c !== 'object' || Array.isArray(c)
    || !sameFamilies(result.families, cp.families) || Math.max(...families.map(k => result.families[k])) <= tolerance
    || f.options.ncrit !== startupPlan.sourceNcrit || f.options.transitionMode !== 'automatic'
    || (f.options.hkFloorLinearization ?? 'exact') !== startupPlan.hkFloorLinearization
    || !terminalTrips(f.options.tripFractions) || f.input.wakeGeometry !== 'independent-banks'
    || f.input.wakeDisplacementMotion !== 'te-center' || !Array.isArray(f.input.bodies)
    || f.input.bodies.length !== f.options.tripFractions.length
    || !Array.isArray(f.options.transitionState) || f.options.transitionState.length !== 2 * f.input.bodies.length
    || !f.options.transitionState.every(x => Number.isInteger(x) && x >= 0)
    || !vector(c.lastRedistributedStagnation) || c.lastRedistributedStagnation.length !== f.input.bodies.length
    || !Array.isArray(c.fractions) || c.fractions.length !== f.input.bodies.length
    || c.fractions.some((row, b) => !vector(row) || row.length !== f.input.bodies[b].leadingIndex + 1
      || row[0] !== 0 || row.at(-1) !== 1 || row.some((v, i) => i > 0 && !(v > row[i - 1])))
    || !Number.isFinite(f.input.mach) || !(f.input.mach > 0 && f.input.mach < 1)
    || !Number.isFinite(f.options.reynolds) || f.options.reynolds <= 0
    || result.conditions?.mach !== f.input.mach || result.conditions?.reynolds !== f.options.reynolds
    || result.conditions?.ncrit !== f.options.ncrit || result.conditions?.transitionMode !== 'automatic'
    || c.blUpdate !== 'xfoil' || c.projectionGeometry !== 'boundary-increment'
    || c.iterationGeometry !== 'ises-sampled' || c.stepAcceptance !== 'admissible'
    || !['listing', 'prose'].includes(c.stagnationLimiter)
    || c.eventProfile !== undefined && c.eventProfile !== 'none'
    || !vector(f.initialEuler?.x) || !vector(f.initialBL) || !vector(result.x) || !vector(result.residual)
    || result.x.length !== f.initialEuler.x.length + f.initialBL.length || result.residual.length !== result.x.length
    || result.residual.reduce((m, v) => Math.max(m, Math.abs(v)), 0) !== Math.max(...families.map(k => result.families[k]))
    || !sameVector(result.x.slice(0, f.initialEuler.x.length), f.initialEuler.x)
    || !sameVector(result.x.slice(f.initialEuler.x.length), f.initialBL)
    || !sameNodes(result.flow?.nodes, f.initialEuler.nodes) || !sameNodes(result.flow?.undisplacedNodes, f.initialEuler.undisplacedNodes)
    || !Array.isArray(history) || (result.progressControl?.stopReason
      ? history.length < 2 || history.length > maxIterations + 1 : history.length !== maxIterations + 1)
    || history.some((h, i) => h?.iteration !== i) || !sameFamilies(history.at(-1), result.families)) return null;
  // Invalid explicit policies reject; omission continues to mean linear.
  if (coupledCheckpointShearCoordinate(cp) !== 'linear') return null;
  const resume = structuredClone(cp);
  resume.continuation.shearCoordinate = 'logarithmic';
  const additionalIterations = Math.min(40, maxIterations);
  return { kind: 'coupled-logarithmic-shear-recovery', sourceShearCoordinate: 'linear', shearCoordinate: 'logarithmic',
    selection: result.progressControl?.stopReason ? 'One bounded logarithmic continuation after detected linear-shear stagnation.'
      : 'One bounded logarithmic continuation after the complete linear Ncrit startup reaches its iteration limit.',
    originalIterations: history.length - 1, additionalIterations, maximumTotalIterations: history.length - 1 + additionalIterations,
    sourceNcrit: f.options.ncrit, targetNcrit: startupPlan.targetNcrit, sourceFamilies: { ...result.families },
    tolerance, equationsChanged: false, physicalStateChanged: false, resume };
}
