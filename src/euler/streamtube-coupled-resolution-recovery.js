// SPDX-License-Identifier: GPL-2.0-or-later
import { recoverCoupledTransition } from './streamtube-transition-recovery.js';

// A bounded change of initial grid, not a change of operating conditions.
// On coarse requests this may increase resolution; larger requests must solve
// again after nested surface refinement. Never accept the coarse root as the
// final answer to a finer request, or inherit its convergence flag.
export function coupledResolutionRecoveryCase(input, result, { enabled = true, maxIterations = 0 } = {}) {
  if (!enabled || maxIterations < 1 || result?.converged || !result?.checkpoint
    || !result.mesh?.quality?.valid || result.resolutionRecovery
    || !(input.mach >= .7) || input.transitionMode !== 'automatic'
    || input.materialTrips?.some(pair => pair.some(t => t !== 1))
    || !Number.isInteger(input.gridIntervals) || input.gridIntervals < 8 || input.gridIntervals > 64
    || input.gridIntervals === 16
    || !['armijo', 'event-armijo'].includes(result.checkpoint.continuation.stepAcceptance)) return null;
  const recent = result.history?.filter(h => h.iteration > 0).slice(-8) ?? [];
  if (recent.length < 4 || !recent.some(h => h.step > 0 && h.step < .02)) return null;
  return { ...structuredClone(input), gridIntervals: 16 };
}

export function mergeCoupledSolveWork(phases) {
  let offset = 0;
  const history = [], iterations = [], diagnostics = { solves: 0, refinements: 0, pivotRecoveries: 0, maxRelativeResidual: 0 };
  for (const { result, label } of phases) {
    const rows = result.history;
    history.push(...rows.slice(history.length ? 1 : 0).map(h => ({ ...h, iteration: h.iteration + offset,
      ...(label ? { resolutionRecovery: label } : {}) })));
    const d = result.linearDiagnostics;
    if (d) {
      for (const key of ['solves', 'refinements', 'pivotRecoveries']) diagnostics[key] += d[key];
      diagnostics.maxRelativeResidual = Math.max(diagnostics.maxRelativeResidual, d.maxRelativeResidual);
      iterations.push(...d.iterations.map(h => ({ ...h, iteration: h.iteration + offset })));
    }
    offset += rows.length - 1;
  }
  return { history, linearDiagnostics: { ...diagnostics, iterations } };
}

export function recoverCoupledResolution(result, input, { sourceInput, solve, maxIterations, eulerMaxIterations,
  tolerance, onIteration, onMesh, onStage, onCheckpoint, onIterationCheckpoint, onFlow } = {}) {
  const diagnostics = { method: 'coupled-resolution-restart', maximumAttempts: 1, accepted: false,
    sourceGridIntervals: sourceInput.gridIntervals, requestedGridIntervals: input.gridIntervals,
    finalNominalGridIntervals: Math.max(input.gridIntervals, sourceInput.gridIntervals),
    requestedGridRetained: false, operatingConditionsChanged: false, equationsChanged: false,
    originalFamilies: result.families, originalIterations: result.history.length - 1 };
  let observerFailure = false, source, refined, offset = result.history.length - 1;
  const observe = callback => (...args) => {
    try { return callback?.(...args); } catch (error) { observerFailure = true; throw error; }
  };
  const root = value => value?.converged === true && value.checkpoint && value.mesh?.quality?.valid === true
    && ['euler', 'boundaryLayer', 'edgeMatching'].every(k => Number.isFinite(value.families?.[k]) && value.families[k] >= 0 && value.families[k] <= tolerance);
  const label = event => ({ ...event, resolutionRecovery: diagnostics, requestedGridIntervals: input.gridIntervals });
  try {
    source = solve(sourceInput, { direct: true, resolutionRecovery: false, machRecovery: false, maxIterations, eulerMaxIterations, tolerance,
      onStage: observe(event => onStage?.(label({ ...event, gridLevel: sourceInput.gridIntervals }))),
      onIteration: observe(h => onIteration?.({ ...h, iteration: h.stage === 'coupled' ? h.iteration + offset : h.iteration,
        gridLevel: sourceInput.gridIntervals, resolutionRecovery: 'startup' })),
      onMesh: observe(onMesh), onCheckpoint: observe(onCheckpoint),
      onIterationCheckpoint: observe((cp, detail) => onIterationCheckpoint?.(cp, { ...detail,
        iterationOffset: (detail.iterationOffset ?? 0) + offset, resolutionRecovery: 'startup' })),
      onFlow: onFlow ? observe(flow => onFlow({ ...flow, iteration: { ...flow.iteration,
        iteration: flow.iteration.iteration + offset }, gridLevel: sourceInput.gridIntervals, resolutionRecovery: 'startup' })) : undefined });
    diagnostics.source = { converged: source.converged, reason: source.reason, families: source.families,
      iterations: source.history?.length - 1, euler: source.initialization?.euler };
    if (!root(source)) return { diagnostics };
    if (input.gridIntervals > sourceInput.gridIntervals) {
      offset += source.history.length - 1;
      const f = source.checkpoint.restart.input, nx = f.outerLower.length - 1;
      const factor = Math.ceil(input.gridIntervals / sourceInput.gridIntervals);
      const counts = Array.from({ length: nx }, (_, i) => f.bodies.some(b => i >= b.leadingIndex && i < b.trailingIndex) ? factor : 1);
      const plan = { reason: 'converged-coupled-resolution-startup', streamwiseSubdivisions: counts,
        normalFactor: 1, maxNodes: 50000, sourceGridIntervals: sourceInput.gridIntervals,
        requestedGridIntervals: input.gridIntervals, requestedGridRetained: false };
      refined = recoverCoupledTransition(source, { plan, maxIterations, tolerance,
        normalization: Object.fromEntries(['referenceChord', 'referenceReynolds', 'solverLength', 'kernelReynolds'].map(k => [k, source[k]])), onStage: observe(event => onStage?.(label({ ...event,
          stage: event.stage === 'transition-refinement' ? 'coupled-grid-refinement' : event.stage,
          transitionRecovery: undefined, refinement: true, gridLevel: input.gridIntervals }))),
        onIteration: observe(h => onIteration?.({ ...h, stage: 'coupled', iteration: h.iteration + offset,
          gridLevel: input.gridIntervals, resolutionRecovery: 'refinement' })),
        onMesh: observe(onMesh),
        onIterationCheckpoint: observe((cp, detail) => onIterationCheckpoint?.(cp, { ...detail,
          iterationOffset: offset, resolutionRecovery: 'refinement' })),
        onFlow: onFlow ? observe(flow => onFlow({ ...flow, iteration: { ...flow.iteration,
          iteration: flow.iteration.iteration + offset }, gridLevel: input.gridIntervals, resolutionRecovery: 'refinement' })) : undefined });
      diagnostics.refinement = { plan, converged: refined.converged, reason: refined.reason,
        families: refined.families, iterations: refined.history.length - 1 };
      if (!root(refined)) return { diagnostics };
    }
  } catch (error) {
    if (observerFailure) throw error;
    return { diagnostics: { ...diagnostics, reason: error.message } };
  }
  const accepted = refined ?? source;
  diagnostics.accepted = true;
  diagnostics.cells = accepted.mesh.cells.length;
  accepted.mesh.initialization ??= {};
  accepted.mesh.initialization.resolutionRecovery = diagnostics;
  observe(onMesh)(accepted.mesh, 'solving', 'coupled');
  const work = mergeCoupledSolveWork([{ result }, { result: source, label: 'startup' },
    ...(refined ? [{ result: refined, label: 'refinement' }] : [])]);
  return { result: { ...accepted, ...work, resolutionRecovery: diagnostics }, diagnostics,
    initialization: source.initialization };
}
