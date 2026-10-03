// SPDX-License-Identifier: GPL-2.0-or-later
// Complete-state nested refinement, including a retained off-target root
// when the caller supplies operating-point continuation on the finer grid.
import { coupledAssemblyConditions } from '../streamtube-coupled-assembly.js';
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { prepareCoupledGridSequence } from '../streamtube-coupled-grid-sequence.js';
import { solveCoupledGridLevel } from '../streamtube-coupled-grid-levels.js';
import { checkpointDataEqual } from '../streamtube-nested-checkpoint.js';

export function targetMachGridPlan(input, { parentResult, preparedEuler, enabled = true } = {}) {
  if (typeof enabled !== 'boolean') throw new Error('Invalid target-Mach grid sequencing policy.');
  if (!enabled || preparedEuler !== undefined || !(input?.mach > .3)
    || input.flowModel !== 'streamtube-grid' || input.quadBoundaryLayers !== true
    || input.transitionMode !== 'automatic' || !Number.isInteger(input.gridIntervals)
    || input.gridIntervals < 32 || input.gridIntervals > 128 || !Array.isArray(input.elements)) return null;
  if (input.materialTrips !== undefined && (!Array.isArray(input.materialTrips)
    || input.materialTrips.length !== input.elements.length
    || !input.materialTrips.every(pair => Array.isArray(pair) && pair.length === 2 && pair.every(v => v === 1)))) return null;
  const coarseIntervals = parentResult?.sourceCase?.gridIntervals ?? (input.gridIntervals === 32 ? 16 : 32);
  if (parentResult !== undefined && (![16, 32, 64].includes(coarseIntervals)
    || coarseIntervals >= input.gridIntervals || !parentResult.checkpoint?.restart
    || !(parentResult.converged || parentResult.stateConverged)
    || !checkpointDataEqual({ ...parentResult.sourceCase, mach: input.mach, alpha: input.alpha,
      gridIntervals: input.gridIntervals }, input))) return null;
  return { method: 'target-mach-before-refinement', coarseIntervals,
    requestedIntervals: input.gridIntervals, targetMach: input.mach };
}

export function solveCoupledTargetGrid(input, { plan, solveCoarse, parentResult, maxIterations = 40,
  tolerance = 1e-10, resumeOperatingPoint, dissipationEnhancement = false, iterationRecovery,
  onStage, onIteration, onMesh, onCheckpoint, onPrepared, onFlow } = {}) {
  if (!plan || plan.method !== 'target-mach-before-refinement' || plan.targetMach !== input.mach
    || plan.requestedIntervals !== input.gridIntervals || ![16, 32, 64].includes(plan.coarseIntervals) || plan.coarseIntervals >= input.gridIntervals
    || !Number.isInteger(maxIterations) || maxIterations < 0 || !(tolerance > 0 && Number.isFinite(tolerance)))
    throw new Error('Invalid target-Mach refinement plan.');
  const coarseCase = { ...structuredClone(input), gridIntervals: plan.coarseIntervals };
  const labels = { targetMach: input.mach, requestedGridIntervals: input.gridIntervals, gridLevel: plan.coarseIntervals,
    gridStrategy: plan.method };
  let observerError, observerFailed = false;
  const observe = (fn, ...args) => {
    try { return fn?.(...args); } catch (e) { observerFailed = true; observerError = e; throw e; }
  };
  const decorate = event => ({ ...event, ...labels });
  let current = parentResult;
  if (current === undefined) {
    if (typeof solveCoarse !== 'function') throw new Error('Missing coarse Mach solver.');
    current = solveCoarse(coarseCase, {
      stopAtShock: true,
      onStage: event => observe(onStage, decorate(event)), onIteration: event => observe(onIteration, decorate(event)),
      onMesh: (mesh, phase, state) => observe(onMesh, mesh, phase, typeof state === 'object' ? decorate(state) : state),
      onFlow: onFlow ? frame => observe(onFlow, decorate(frame)) : undefined, onPrepared: p => observe(onPrepared, p),
      // The worker must never cache a coarse root under the requested fine controls.
      onCheckpoint: (cp, details) => observe(onCheckpoint, cp, { ...decorate(details),
        ...(details.kind === 'accepted' ? { kind: 'coarse-accepted', reachedTarget: false } : {}) }),
    });
  } else if (!checkpointDataEqual(current.sourceCase, resumeOperatingPoint
    ? { ...coarseCase, mach: current.checkpoint?.restart.input.mach,
      alpha: current.checkpoint?.restart.input.alpha ?? 0 } : coarseCase)) {
    throw new Error('Target-grid source controls differ from the requested coarse case.');
  }
  const levels = [], sequence = () => ({ ...plan, levels: structuredClone(levels),
    actualGridIntervals: nominal, requestedGridIntervals: input.gridIntervals,
    reachedTarget: nominal >= input.gridIntervals && current.converged === true && atTarget(current) });
  let nominal = plan.coarseIntervals;
  const retained = reason => ({ ...current, converged: false, stateConverged: current.stateConverged ?? current.converged === true,
    requestedCase: structuredClone(input), gridSequence: sequence(), automaticRefinement: sequence(),
    status: 'research-coupled-target-not-reached', reason });
  const rootAvailable = value => value?.mesh?.quality?.valid && value.checkpoint?.restart
    && (value.converged === true || value.stateConverged === true)
    && ['euler', 'boundaryLayer', 'edgeMatching'].every(k =>
      Number.isFinite(value.checkpoint.families?.[k]) && value.checkpoint.families[k] >= 0 && value.checkpoint.families[k] <= tolerance);
  const atTarget = value => value.checkpoint?.restart.input.mach === input.mach
    && (value.checkpoint.restart.input.alpha ?? 0) === (input.alpha ?? 0);
  if ((!current.converged || current.mach !== input.mach || !current.mesh?.quality?.valid)
    && !(resumeOperatingPoint && rootAvailable(current)))
    return retained(`Coarse target-Mach solve did not converge: ${current.reason}`);
  while (nominal < input.gridIntervals) {
    const next = Math.min(2 * nominal, input.gridIntervals), cp = current.checkpoint, f = cp?.restart;
    if (!f) throw new Error('Target-grid source requires a complete converged checkpoint.');
    const actualMach = f.input.mach, actualAlpha = f.input.alpha ?? 0;
    const sourceCase = current.sourceCase ?? input;
    if (f.input.mach !== (sourceCase.mach ?? input.mach) || (f.input.alpha ?? 0) !== (sourceCase.alpha ?? 0))
      throw new Error('Target-grid checkpoint Mach or incidence does not match the requested operating point.');
    const settings = coupledAssemblyConditions(input, f.input.bodies, current.solverLength, { checkpointOptions: f.options });
    for (const k of ['reynolds', 'ncrit', 'edgeMatching', 'tripFractions'])
      if (!checkpointDataEqual(settings.options[k], f.options[k])) throw new Error(`Target-grid source changed ${k}.`);
    for (const [k, v] of Object.entries(settings.normalization))
      if (v !== current[k]) throw new Error(`Target-grid source changed ${k}.`);
    const source = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
    const level = { sourceIntervals: nominal, requestedIntervals: next, mach: actualMach, alpha: actualAlpha, converged: false };
    levels.push(level);
    const gridStrategy = current.refinementRequested ? 'shock-triggered-refinement' : atTarget(current) ? plan.method : 'refine-retained-operating-point';
    try {
      observe(onStage, { stage: 'coupled-grid-refinement', ...labels, gridLevel: next, mach: actualMach,
        actualMach, actualAlpha, targetAlpha: input.alpha,
        gridStrategy, dissipationEnhancement });
      const prepared = prepareCoupledGridSequence({ sourceSystem: source, sourceResult: { ...current, converged: true },
        input: f.input, options: f.options, requestedGridIntervals: next, sourceGridIntervals: nominal },
      { tolerance, blPredictor: source.bl.hasFiniteBase ? 'interpolate' : 'xfoil-mrchdu' });
      const normalization = Object.fromEntries(['referenceChord', 'referenceReynolds', 'solverLength', 'kernelReynolds'].map(k => [k, current[k]]));
      observe(onPrepared, { system: prepared.system, settings: { normalization }, checkpoint: undefined,
        parentResult: { ...current, sourceCase: { ...structuredClone(input), mach: actualMach, alpha: actualAlpha, gridIntervals: next } } });
      const solved = solveCoupledGridLevel(prepared, { nominalGridIntervals: prepared.transfer.nominalGridIntervals,
        targetGridIntervals: input.gridIntervals }, { maxIterations, tolerance, normalization, startupAttempt: 1, dissipationEnhancement, iterationRecovery,
        onStage: e => observe(onStage, { ...e, mach: actualMach, actualAlpha, targetAlpha: input.alpha, targetMach: input.mach, gridStrategy }),
        onIteration: e => observe(onIteration, { ...e, mach: actualMach, actualAlpha, targetAlpha: input.alpha, targetMach: input.mach }),
        onMesh: (mesh, phase, stage, snapshot) => {
          mesh.initialization.gridSmoothing = current.initialization?.euler?.gridSmoothing;
          const { bodyPressure, ...flow } = snapshot?.flow ?? {};
          observe(onMesh, mesh, phase, snapshot ? { flow, iteration: snapshot.iteration,
            coupledFamilies: snapshot.coupledFamilies, mach: actualMach, actualAlpha, targetAlpha: input.alpha, targetMach: input.mach,
            stage, startupAttempt: 1, gridLevel: next } : stage);
        },
        onFlow: onFlow ? frame => observe(onFlow, { ...frame, actualAlpha, targetAlpha: input.alpha, targetMach: input.mach }) : undefined,
        onIterationCheckpoint: (value, details) => observe(onCheckpoint, value, { ...details, mach: actualMach, actualAlpha, targetAlpha: input.alpha,
          targetMach: input.mach, kind: 'iterate' }),
      });
      if (observerFailed) throw observerError;
      if (!solved.checkpoint) throw new Error(solved.failure?.reason ?? solved.reason ?? 'Refined state has no checkpoint.');
      Object.assign(level, { actualIntervals: prepared.transfer.nominalGridIntervals, converged: solved.converged, families: solved.families,
        reason: solved.reason, iterations: solved.linearDiagnostics.solves, cells: solved.mesh.cells.length });
      if (!solved.converged) {
        observe(onMesh, current.mesh, 'retained', { stage: 'coupled-grid-refinement', retained: true,
          mach: actualMach, actualAlpha, targetMach: input.mach, targetAlpha: input.alpha });
        return retained(`Target-Mach grid refinement did not converge: ${solved.reason}`);
      }
      nominal = prepared.transfer.nominalGridIntervals;
      current = { ...current, ...solved, refinementRequested: undefined, model: 'research-streamtube-euler-bl', ...normalization,
        sourceCase: { ...structuredClone(input), mach: actualMach, alpha: actualAlpha, gridIntervals: nominal }, requestedCase: structuredClone(input),
        mach: actualMach, actualMach, alpha: actualAlpha, actualAlpha, targetAlpha: input.alpha, targetMach: input.mach, stateConverged: solved.converged,
        bodies: solved.solverInput.bodies, restart: solved.checkpoint.restart,
        solverSettings: { ...current.solverSettings, maxIterations, tolerance, dissipationEnhancement },
        initialization: { ...current.initialization, targetMachGridSequence: { ...plan, levels: structuredClone(levels) } } };
      if (!atTarget(current) && resumeOperatingPoint) {
        const resumed = resumeOperatingPoint(current, {
          stopAtShock: next < input.gridIntervals,
          onStage: e => observe(onStage, { ...e, gridLevel: nominal }),
          onIteration: e => observe(onIteration, { ...e, gridLevel: nominal }),
          onMesh: (mesh, phase, info) => observe(onMesh, mesh, phase, info),
          onFlow: onFlow ? frame => observe(onFlow, frame) : undefined,
          onCheckpoint: (cp, e) => observe(onCheckpoint, cp, { ...e,
            kind: e.kind === 'accepted' && nominal < input.gridIntervals ? 'coarse-accepted' : e.kind,
            reachedTarget: nominal >= input.gridIntervals && e.reachedTarget }),
        });
        if (!rootAvailable(resumed)) throw new Error('Operating-point continuation did not retain a converged root.');
        current = resumed;
        level.operatingPointContinuation = resumed.operatingPointContinuation;
      }
      observe(onCheckpoint, current.checkpoint, { stage: 'coupled-grid-refinement', mach: current.mach, actualAlpha: current.actualAlpha,
        targetMach: input.mach, gridLevel: nominal, requestedGridIntervals: input.gridIntervals,
        kind: nominal >= input.gridIntervals && current.converged && atTarget(current) ? 'accepted' : 'coarse-accepted', reachedTarget: nominal >= input.gridIntervals && current.converged && atTarget(current) });
    } catch (error) {
      if (observerFailed) throw observerError;
      level.reason = error.message;
      return retained(`Requested ${input.gridIntervals}-interval grid was not reached; retained ${nominal}-interval state. ${error.message}`);
    }
  }
  if (!current.converged || !atTarget(current))
    return retained(`Requested operating point was not reached on the refined grid. ${current.reason ?? ''}`);
  return { ...current, sourceCase: structuredClone(input), gridSequence: sequence(), automaticRefinement: sequence(),
    machContinuation: { ...current.machContinuation, actualMach: input.mach, targetMach: input.mach, reachedTarget: true } };
}
