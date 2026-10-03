// SPDX-License-Identifier: GPL-2.0-or-later
// Opt-in MPOLAR-style Mach continuation from a complete accepted Euler ISES
// checkpoint. Only Mach changes; every stage retains density, physical mass,
// inlet fractions and the last-SMOVE reference. No new initial SMOVE.
import { createStreamtubeBodySystem } from '../streamtube-body.js';
import { solveStreamtubeIses } from '../streamtube-ises-update.js';
import { initializeStreamtubeBodyFromFlow } from './streamtube-body-flow-restart.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';

export function solveStreamtubeIsesAutomatic(targetMach, { initialCheckpoint, tolerance = 1e-10,
  stageMaxIterations = 12, maxSubdivisions = 5, maxStages = 32,
  onStage, onIteration, onMesh, onCheckpoint } = {}) {
  if (!Number.isFinite(targetMach) || targetMach <= 0 || targetMach >= 1
    || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isInteger(stageMaxIterations) || stageMaxIterations < 0
    || !Number.isInteger(maxSubdivisions) || maxSubdivisions < 0 || maxSubdivisions > 20
    || !Number.isInteger(maxStages) || maxStages < 1 || maxStages > 256
    || [onStage, onIteration, onMesh, onCheckpoint].some(f => f !== undefined && typeof f !== 'function'))
    throw new Error('Invalid automatic Euler ISES continuation controls.');
  if (initialCheckpoint?.version !== 1 || !initialCheckpoint.input || !initialCheckpoint.continuation)
    throw new Error('Automatic Euler ISES requires a complete accepted checkpoint.');
  initialCheckpoint = structuredClone(initialCheckpoint);
  const c = initialCheckpoint.continuation;
  const controls = { tolerance, stagnationLimiter: c.stagnationLimiter,
    iterationGeometry: c.iterationGeometry, stepAcceptance: c.stepAcceptance };
  // Validate before recovery: a corrupt/nonconverged/folded source cannot
  // become an accepted origin merely by reducing the requested increment.
  let initial = solveStreamtubeIses(undefined, { ...controls, resume: initialCheckpoint, maxIterations: 0 });
  if (!initial.converged || !initial.finalQuality.valid)
    throw new Error('Automatic Euler ISES source must be converged on a valid convex grid.');
  const record = result => {
    const checkpoint = result.checkpoint;
    const system = createStreamtubeBodySystem(checkpoint.input);
    // ISES builds a fresh normal chart after each accepted step. Never pair
    // its final x with the stale preparation chart from an earlier stage.
    const state = system.adoptGeometry(Float64Array.from(checkpoint.initialEuler.x), checkpoint.initialEuler.nodes);
    return { result, checkpoint, system, state };
  };
  let current = record(initial);
  initial = null; initialCheckpoint = null;
  const sourceMach = current.system.conditions.mach;
  if (!(sourceMach > 0 && sourceMach < 1) || current.system.conditions.flowModel !== 'compressible')
    throw new Error('Automatic Euler ISES requires a compressible subsonic-freestream source.');
  const attempts = [];
  const formulation = () => `Research Euler ISES automatic freestream-Mach continuation; ${current.system.conditions.streamwiseMode} streamwise rows${current.system.conditions.hybrid ? ` with explicit research momentum/artificial-entropy blend epsilonP=${current.system.conditions.hybrid.epsilonP}` : ''}${current.system.conditions.upwind ? '; shared MSES speed upwinding, locally transonic interior and subsonic remote endpoints' : '; no shock dissipation'}; additive-density Newton, inlet adjustment, DEKINK and stagnation-triggered five-pair SMOVE; complete maintenance checkpoint retained, initial SMOVE not repeated; final convexity and residual required; no evolving BL or physical shock acceptance.`;
  const finish = (reachedTarget, reason) => ({ ...current.result,
    converged: reachedTarget && current.result.converged, reason,
    stateConverged: true,
    status: reachedTarget ? 'research-euler-ises-converged' : 'research-euler-ises-target-not-reached',
    checkpoint: structuredClone(current.checkpoint),
    continuation: { method: 'freestream-mach-ises', sourceMach, targetMach,
      currentMach: current.system.conditions.mach, reachedTarget, attempts, stageMaxIterations, maxSubdivisions, maxStages,
      initialSMOVERepeated: false, densityReinitialized: false,
      analyticShockInformationUsed: false, acceptedSourceReplayExact: true },
    formulation: formulation(), physicalAcceptance: false, fullSolver: false });
  if (sourceMach === targetMach) return finish(true, 'residual');

  const direction = Math.sign(targetMach - sourceMach);
  let step = Math.abs(targetMach - sourceMach), subdivisions = 0;
  while (current.system.conditions.mach !== targetMach) {
    if (attempts.length >= maxStages) return finish(false, 'Automatic Euler ISES stage limit.');
    const fromMach = current.system.conditions.mach, remaining = Math.abs(targetMach - fromMach);
    const mach = step >= remaining ? targetMach : fromMach + direction * step;
    if (mach === fromMach) return finish(false, 'Automatic Euler ISES Mach increment is not representable.');
    const increment = Math.abs(mach - fromMach), attempt = attempts.length + 1;
    const metadata = { stage: mach === targetMach ? 'target' : 'mach-continuation',
      attempt, fromMach, mach, sourceMach, targetMach,
      continuationFraction: (mach - sourceMach) / (targetMach - sourceMach) };
    const history = []; let prepared, transfer, result, observerFailed = false, observerError;
    const observe = (callback, ...values) => {
      if (!callback) return;
      try { callback(...values.map(v => structuredClone(v))); }
      catch (error) { observerFailed = true; observerError = error; throw error; }
    };
    try {
      observe(onStage, metadata);
      prepared = initializeStreamtubeBodyFromFlow({ ...current.checkpoint.input, mach }, current.system,
        { initial: current.state });
      // A changed Mach has a different h0 and residual. Recompute those
      // physical target values, while retaining the source maintenance
      // history. Editing only checkpoint.input.mach would fail exact replay.
      const resume = { version: 1, input: prepared.input,
        initialEuler: { x: Array.from(prepared.initial), nodes: prepared.initialEuler.nodes },
        residual: Array.from(prepared.flow.residual),
        continuation: structuredClone(current.checkpoint.continuation) };
      transfer = prepared.diagnostics;
      // The target resume contains detached coordinates and scalar data;
      // release its temporary full flow/system before allocating the solve.
      prepared = null;
      result = solveStreamtubeIses(undefined, { ...controls, resume, maxIterations: stageMaxIterations,
        onIteration: entry => { history.push(structuredClone(entry)); observe(onIteration, { ...entry, ...metadata }); },
        onMesh: onMesh ? snapshot => observe(onMesh, { ...metadata,
          mesh: streamtubeMeshSnapshot(snapshot), iteration: snapshot.iteration }) : undefined,
        onCheckpoint: onCheckpoint ? (checkpoint, details) => observe(onCheckpoint, checkpoint, { ...details, ...metadata, kind: 'iterate' }) : undefined });
      // Existing numerical drivers may catch an observer error during an
      // accepted update. Never turn that cancellation into a Mach retry.
      if (observerFailed) throw observerError;
      attempts.push({ ...metadata, converged: result.converged, reason: result.reason,
        iterations: result.history.length - 1, linearSolves: result.linearDiagnostics.solves,
        residual: result.diagnostics.residual, gridValid: result.finalQuality.valid,
        initialSMOVERepeated: !result.initialRedistribution.resumed,
        transfer, history });
    } catch (error) {
      if (observerFailed) throw observerError;
      attempts.push({ ...metadata, converged: false, reason: error.message,
        iterations: Math.max(0, history.length - 1), history,
        ...(transfer ? { transfer } : {}), ...(error.code ? { code: error.code } : {}) });
      result = null;
    }
    if (result?.converged && result.finalQuality.valid) {
      current = record(result); subdivisions = 0;
      observe(onCheckpoint, current.checkpoint, { ...metadata, kind: 'accepted',
        residual: result.diagnostics.residual, quality: result.finalQuality,
        stateConverged: true, reachedTarget: mach === targetMach });
      step = Math.min(1.5 * increment, Math.abs(targetMach - mach));
    } else {
      if (subdivisions >= maxSubdivisions)
        return finish(false, `Automatic Euler ISES stopped at Mach ${current.system.conditions.mach}.`);
      // The failed terminal iterate is never used as a continuation source.
      step = .5 * increment; subdivisions++;
    }
  }
  return finish(true, 'residual');
}
