// SPDX-License-Identifier: GPL-2.0-or-later
import { observableFlow, observableBL } from './streamtube-flow-preview.js';
// Bounded startup work selection. Every retained continuation source remains
// a complete root; an intermediate Ncrit never qualifies the requested case.
import { coupledCheckpointShearCoordinate, coupledShearCoordinate } from './streamtube-coupled-shear-policy.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { continueCoupledNcrit } from './streamtube-coupled-ncrit-continuation.js';
import { recoverCoupledTransition } from './streamtube-transition-recovery.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';

const maximum = f => f && ['euler', 'boundaryLayer', 'edgeMatching'].every(k => Number.isFinite(f[k]) && f[k] >= 0)
  ? Math.max(f.euler, f.boundaryLayer, f.edgeMatching) : Infinity;

export function coupledNcritStartupRecoveryPlan(input, result, { startupAttempt, maxStartupAttempts,
  maxIterations, thicknessFactor, coarseInitialization } = {}) {
  const options = result?.checkpoint?.restart?.options, gas = result?.checkpoint?.restart?.input;
  if (result?.converged || startupAttempt !== 1 || maxStartupAttempts !== 2 || !(maxIterations > 0)
    || !Number.isInteger(maxIterations) || coarseInitialization || result?.automaticRefinement
    || !result?.mesh?.quality?.valid || !Number.isFinite(maximum(result?.families))
    || options?.transitionMode !== 'automatic' || !Array.isArray(options.tripFractions)
    || !options.tripFractions.length || options.tripFractions.length !== gas?.bodies?.length
    || !options.tripFractions.every(pair => Array.isArray(pair) && pair.length === 2 && pair.every(x => x === 1))
    || input?.coupledNativeHk !== undefined && typeof input.coupledNativeHk !== 'boolean'
    || ![undefined, 'exact', 'native'].includes(options?.hkFloorLinearization)
    || !Number.isFinite(input?.ncrit) || !(input.ncrit > 4) || options.ncrit !== input.ncrit || gas?.mach !== input.mach
    || gas?.wakeGeometry !== 'independent-banks' || gas.wakeDisplacementMotion !== 'te-center'
    || !(thicknessFactor > 0 && thicknessFactor <= 1)) return null;
  return { kind: 'ncrit-transition-startup', sourceNcrit: 4, targetNcrit: input.ncrit,
    // This selects the original-XFOIL derivative continuation only for
    // the existing second startup; the governing residuals are unchanged.
    sourceHkFloorLinearization: options.hkFloorLinearization ?? 'exact',
    // Explicit numerical policy for this existing recovery only. Ordinary
    // starts stay linear, independently of their Hk derivative policy.
    sourceShearCoordinate: coupledCheckpointShearCoordinate(result.checkpoint),
    // Preserve the original thin-seed transient. A separate bounded
    // continuation can change coordinates after this linear startup stalls.
    shearCoordinate: 'linear',
    hkFloorLinearization: input.coupledNativeHk === false ? 'exact' : 'native',
    thicknessFactor: .25 * thicknessFactor, coarseAdvanceNcrit: Math.min(4 + 1.5, input.ncrit),
    sourceIterations: maxIterations, coarseAdvanceIterations: Math.min(maxIterations, 20),
    refinedIterations: Math.min(maxIterations, 20), maxNodes: 50000,
    selection: 'Bounded initial-guess heuristic: advance at most 1.5 in Ncrit before refining natural-transition windows.',
    equationsChanged: false };
}

export function coupledStartupTransitionPlan(result, { tolerance = 1e-10, maxNodes = 50000 } = {}) {
  if (!Number.isFinite(tolerance) || tolerance <= 0 || !result?.converged || maximum(result.families) > tolerance
    || !Number.isInteger(maxNodes) || maxNodes < 1 || maxNodes > 50000)
    throw new Error('Startup transition refinement requires a strict complete root and a bounded node budget.');
  const f = result.checkpoint.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const value = system.evaluate(system.initial);
  if (Object.keys(result.families).some(k => value.families[k] !== result.families[k]))
    throw new Error('Startup transition root does not replay exactly.');
  const parentNx = system.euler.layout.nx, counts = Array(parentNx).fill(1);
  const surfaces = value.layers.transitions.filter(t => t.kind === 'natural').map(t => {
    const station = system.bl.stations[t.id], body = system.euler.layout.bodies[t.body];
    const first = Math.max(body.leadingIndex, station.i - 3), last = Math.min(body.trailingIndex - 1, station.i + 1);
    for (let i = first; i <= last; i++) counts[i] = 4;
    return { body: t.body, side: t.side, transition: t.s, transitionStation: station.i,
      intervals: Array.from({ length: last - first + 1 }, (_, i) => first + i) };
  });
  if (!surfaces.length) return null;
  const refinedNx = counts.reduce((sum, n) => sum + n, 0);
  const nodeCount = (refinedNx + 1) * system.euler.layout.tubes.reduce((sum, n) => sum + n + 1, 0);
  if (nodeCount > maxNodes) throw new Error(`Startup transition refinement needs ${nodeCount} nodes, exceeding ${maxNodes}.`);
  return { reason: 'ncrit-startup-transition-resolution', source: 'converged-ncrit-root',
    surfaces, streamwiseSubdivisions: counts, normalFactor: 1, maxNodes, nodeCount, parentNx, refinedNx };
}

export function finishCoupledNcritStartup(source, { plan, maxIterations, tolerance, normalization, startupAttempt,
  onStage, onIteration, onMesh, onFlow, onIterationCheckpoint } = {}) {
  if (!plan || !Number.isInteger(maxIterations) || maxIterations < 1 || !Number.isFinite(tolerance) || tolerance <= 0
    || !source?.converged || maximum(source.families) > tolerance
    || source.checkpoint?.restart?.options?.ncrit !== plan.sourceNcrit)
    throw new Error('Ncrit startup completion requires its converged intermediate source.');
  if (plan.hkFloorLinearization !== undefined && (!['exact', 'native'].includes(plan.hkFloorLinearization)
    || plan.hkFloorLinearization !== (source.checkpoint.restart.options.hkFloorLinearization ?? 'exact')))
    throw new Error('Ncrit startup plan and source checkpoint disagree on Hk-floor linearization.');
  if (plan.shearCoordinate !== undefined
    && coupledShearCoordinate(plan.shearCoordinate) !== coupledCheckpointShearCoordinate(source.checkpoint))
    throw new Error('Ncrit startup plan and source checkpoint disagree on shear coordinate.');
  const diagnostics = { ...plan, source: { converged: true, families: { ...source.families },
    iterations: source.history.length - 1 }, reachedTarget: false, attempted: true };
  let retained = source, observerFailed = false, observerError;
  let publishedNcrit = source.checkpoint.restart.options.ncrit;
  const observe = fn => (...args) => { try { return fn?.(...args); }
    catch (error) { observerFailed = true; observerError = error; throw error; } };
  const labels = actual => ({ actualNcrit: actual, targetNcrit: plan.targetNcrit, startupRecovery: true,
    ncritContinuation: { actualNcrit: actual, targetNcrit: plan.targetNcrit, reachedTarget: actual === plan.targetNcrit } });
  const relabel = (packet, actual) => ({ ...packet, ...labels(actual),
    ncritContinuation: { ...packet.ncritContinuation, ...labels(actual).ncritContinuation } });
  const hooks = {
    onStage: observe(event => onStage?.(relabel(event, event.actualNcrit ?? retained.conditions.ncrit))),
    onIteration: observe(h => onIteration?.(relabel({ ...h, stage: h.stage ?? 'coupled' }, h.actualNcrit ?? publishedNcrit))),
    onMesh: observe((mesh, ...details) => onMesh?.({ ...mesh, ...(mesh.iteration ? {
      iteration: relabel({ ...mesh.iteration, stage: mesh.iteration.stage ?? 'coupled' }, publishedNcrit),
    } : {}) }, ...details)),
    onFlow: onFlow ? observe(frame => onFlow(relabel(frame, frame.checkpoint.restart.options.ncrit))) : undefined,
    onIterationCheckpoint: observe((checkpoint, details) => {
      // A continuation retry republishes its retained checkpoint before its
      // mesh. Keep that mesh's labels tied to the same complete state.
      publishedNcrit = checkpoint.restart.options.ncrit;
      onIterationCheckpoint?.(checkpoint, relabel(details, publishedNcrit));
    }),
  };
  const restore = () => {
    const f = retained.checkpoint.restart, system = createCoupledStreamtubeBody(f.input, { ...f.options,
      initialEuler: f.initialEuler, initialBL: f.initialBL });
    const frame = { system: system.euler, flow: retained.flow, nodes: retained.flow.nodes,
      iteration: { ...retained.history.at(-1), ...labels(f.options.ncrit) } };
    const mesh = streamtubeMeshSnapshot(frame);
    mesh.initialization = { ...mesh.initialization, ...retained.mesh.initialization };
    hooks.onIterationCheckpoint(structuredClone(retained.checkpoint), { stage: 'coupled', startupAttempt, retained: true });
    hooks.onMesh(mesh, 'solving', 'coupled');
    hooks.onFlow?.(structuredClone({ checkpoint: retained.checkpoint,
      flow: observableFlow(retained.flow),
      bl: observableBL(system.bl),
      bodies: system.euler.layout.bodies, normalization, iteration: frame.iteration,
      stage: 'coupled', startupAttempt, mach: f.input.mach }));
  };
  try {
    const advance = continueCoupledNcrit(retained, { targetNcrit: plan.coarseAdvanceNcrit, phase: 'fine',
      maxIterations: plan.coarseAdvanceIterations, tolerance, normalization, startupAttempt, ...hooks });
    retained = advance.result; diagnostics.coarseAdvance = advance.diagnostics;
    if (!advance.diagnostics.reachedTarget) throw new Error('The coarse Ncrit startup increment did not converge.');
    if (retained.conditions.ncrit < plan.targetNcrit) {
      const refinement = coupledStartupTransitionPlan(retained, { tolerance, maxNodes: plan.maxNodes });
      if (refinement) {
        const refined = recoverCoupledTransition(retained, { plan: refinement, maxIterations: plan.refinedIterations,
          tolerance, blPredictor: 'xfoil-mrchdu', normalization, startupAttempt, ...hooks });
        const refinedRestart = refined.checkpoint?.restart;
        diagnostics.refinement = { ...refinement, converged: refined.converged && !!refinedRestart, reason: refined.reason,
          families: { ...refined.families }, iterations: refined.history.length - 1,
          cells: refined.mesh.cells.length, ...(refinedRestart
            ? { actualIntervals: refinedRestart.input.outerLower.length - 1 } : { checkpointAvailable: false }) };
        if (!refined.converged || !refinedRestart) throw Object.assign(
          new Error(`The refined Ncrit startup root did not converge: ${refined.reason}`), {
            ...(refined.lastRejectedStep?.code ? { code: refined.lastRejectedStep.code } : {}),
            diagnostics: structuredClone({ reason: refined.reason, families: refined.families,
              initialRedistribution: refined.initialRedistribution, lastRejectedStep: refined.lastRejectedStep,
              meshQuality: refined.mesh.quality, checkpointAvailable: !!refinedRestart }),
          });
        retained = refined;
      }
      const continued = continueCoupledNcrit(retained, { targetNcrit: plan.targetNcrit, phase: 'fine',
        maxIterations, tolerance, normalization, startupAttempt, ...hooks });
      retained = continued.result; diagnostics.fine = continued.diagnostics;
    }
  } catch (error) {
    if (observerFailed) throw observerError;
    diagnostics.failure = { message: error?.message ?? String(error), code: error?.code, diagnostics: error?.diagnostics };
    restore();
  }
  diagnostics.actualNcrit = retained.checkpoint.restart.options.ncrit;
  diagnostics.reachedTarget = retained.converged && diagnostics.actualNcrit === plan.targetNcrit;
  diagnostics.stateConverged = retained.converged;
  const count = continuation => (continuation?.attempts ?? []).reduce((sum, attempt) => sum + (attempt.iterations ?? 0), 0);
  diagnostics.iterations = diagnostics.source.iterations + count(diagnostics.coarseAdvance)
    + (diagnostics.refinement?.iterations ?? 0) + count(diagnostics.fine);
  return { result: retained, diagnostics };
}
