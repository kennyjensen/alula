// SPDX-License-Identifier: GPL-2.0-or-later
import { initializeCoupledStreamtubeFromFlow } from './tests/streamtube-coupled-flow-restart.js';
import { solveCoupledStreamtubeIses } from './streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { observableFlow, observableBL } from './streamtube-flow-preview.js';
import { mergeCoupledSolveWork } from './streamtube-coupled-resolution-recovery.js';

export function coupledMachRecoveryCase(input, result, { enabled = true, maxIterations = 0 } = {}) {
  if (!enabled || maxIterations < 1 || result?.converged || !result?.checkpoint
    || !result.mesh?.quality?.valid || result.machRecovery
    || !(input.mach >= .72 && input.mach <= .85)
    || result.checkpoint.restart.options.edgeMatching !== 'section-velocity'
    || input.transitionMode !== 'automatic'
    || input.materialTrips?.some(pair => pair.some(t => t !== 1))) return null;
  return { ...structuredClone(input), mach: Number((input.mach - .02).toFixed(12)) };
}

// At most one cold source and two complete coupled transfers. Intermediate
// roots are initial guesses only; a failed target leaves the original result.
export function recoverCoupledMach(original, input, { sourceInput, solve, maxIterations, eulerMaxIterations,
  tolerance, normalization, onStage, onIteration, onMesh, onFlow, onIterationCheckpoint } = {}) {
  const diagnostics = { method: 'coupled-mach-restart', accepted: false, maximumTransfers: 2,
    sourceMach: sourceInput.mach, targetMach: input.mach, attempts: [] };
  let observerFailed = false, offset = original.history.length - 1;
  const observe = fn => (...args) => { try { return fn?.(...args); }
    catch (error) { observerFailed = true; throw error; } };
  const label = mach => ({ mach, actualMach: mach, targetMach: input.mach,
    actualAlpha: input.alpha ?? 0, targetAlpha: input.alpha ?? 0, machRecovery: diagnostics });
  const root = r => r?.converged === true && r.checkpoint && r.mesh?.quality?.valid === true
    && ['euler', 'boundaryLayer', 'edgeMatching'].every(k => Number.isFinite(r.families?.[k]) && r.families[k] >= 0 && r.families[k] <= tolerance);
  const phases = [{ result: original }];
  let retained;
  try {
    retained = solve(sourceInput, { direct: true, machRecovery: false, resolutionRecovery: false,
      maxIterations, eulerMaxIterations, tolerance,
      onStage: observe(e => onStage?.({ ...e, ...label(sourceInput.mach) })),
      onIteration: observe(h => onIteration?.({ ...h, iteration: h.stage === 'coupled' ? h.iteration + offset : h.iteration,
        ...label(sourceInput.mach) })), onMesh: observe(onMesh),
      onFlow: onFlow ? observe(f => onFlow({ ...f, iteration: { ...f.iteration, iteration: f.iteration.iteration + offset }, ...label(sourceInput.mach) })) : undefined,
      onIterationCheckpoint: observe((cp, d) => onIterationCheckpoint?.(cp, { ...d,
        iterationOffset: (d.iterationOffset ?? 0) + offset, ...label(sourceInput.mach) })) });
    diagnostics.source = { converged: retained.converged, families: retained.families, reason: retained.reason,
      iterations: retained.history.length - 1 };
    if (!root(retained)) return { diagnostics };
    phases.push({ result: retained, label: 'mach-startup' });
    for (let stage = 1; stage <= 2; stage++) {
      offset += retained.history.length - 1;
      const mach = stage === 2 ? input.mach : Number(((sourceInput.mach + input.mach) / 2).toFixed(12));
      const seed = initializeCoupledStreamtubeFromFlow(mach, retained.checkpoint, { tolerance });
      const f = seed.checkpoint.restart;
      const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
      let checkpoint;
      observe(onStage)({ stage: 'coupled', ...label(mach) });
      const candidate = solveCoupledStreamtubeIses(undefined, { ...seed.checkpoint.continuation,
        resume: seed.checkpoint, maxIterations, tolerance,
        onIteration: observe(h => onIteration?.({ ...h, stage: 'coupled', iteration: h.iteration + offset, ...label(mach) })),
        onCheckpoint: observe((cp, d) => { checkpoint = cp;
          onIterationCheckpoint?.(cp, { ...d, stage: 'coupled', iterationOffset: offset, ...label(mach) }); }),
        onMesh: observe(state => {
          const mesh = streamtubeMeshSnapshot(state);
          onMesh?.(mesh, 'solving', 'coupled');
          if (onFlow && checkpoint) onFlow({ checkpoint, flow: observableFlow(state.flow), bl: observableBL(system.bl),
            bodies: state.system.layout.bodies, normalization, stage: 'coupled', ...label(mach),
            iteration: { ...state.iteration, iteration: state.iteration.iteration + offset } });
        }) });
      diagnostics.attempts.push({ mach, converged: candidate.converged, families: candidate.families,
        iterations: candidate.history.length - 1, reason: candidate.reason });
      if (!root(candidate)) return { diagnostics };
      retained = candidate;
      phases.push({ result: candidate, label: 'mach-continuation' });
    }
  } catch (error) {
    if (observerFailed) throw error;
    diagnostics.reason = error.message;
    return { diagnostics };
  }
  diagnostics.accepted = true;
  return { diagnostics, result: { ...retained, ...mergeCoupledSolveWork(phases), machRecovery: diagnostics } };
}
