// SPDX-License-Identifier: GPL-2.0-or-later
// Automatic target-Mach startup for the conservative moving-grid body.
// Follow MPOLAR's last-converged restart and failed-increment subdivision
// (MSES manual §2.5), without specifying a sonic branch or shock position.
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../streamtube-body.js';
import { initializeStreamtubeStartup } from '../streamtube-startup.js';
import { initializeStreamtubeBodyFromFlow } from './streamtube-body-flow-restart.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';

export function solveStreamtubeBodyAutomatic(input, { initialGeometry, initialEuler, initialFlow, tolerance = 1e-10,
  directMaxIterations = 12, stageMaxIterations = 20, maxSubdivisions = 5,
  initialFractionStep = .25, maxStages = 64, onIteration, onMesh, onStage } = {}) {
  if (![directMaxIterations, stageMaxIterations, maxSubdivisions].every(v => Number.isInteger(v) && v >= 0)
    || maxSubdivisions > 20 || !Number.isInteger(maxStages) || maxStages < 2
    || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isFinite(initialFractionStep) || initialFractionStep <= 0 || initialFractionStep > 1)
    throw new Error('Invalid automatic body continuation controls.');
  if (input.flowModel !== undefined && input.flowModel !== 'compressible'
    || input.streamwiseMode !== undefined && !['momentum', 'hybrid'].includes(input.streamwiseMode))
    throw new Error('Automatic transonic body flow requires compressible conservative momentum or hybrid rows.');
  if (initialFlow !== undefined && (initialGeometry !== undefined || initialEuler !== undefined))
    throw new Error('initialFlow is mutually exclusive with initialGeometry and initialEuler.');
  if (initialEuler !== undefined) throw new Error('Use initialGeometry for a geometry-only seed; complete flow restarts require initializeStreamtubeBodyFromFlow.');
  input = structuredClone({ ...input, flowModel: 'compressible', streamwiseMode: input.streamwiseMode ?? 'momentum',
    upwind: input.upwind ?? { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } });
  // This optional seed supplies geometry and global coordinates only.
  // Densities are initialized for a cold solve; do not accept a full-flow
  // option under a name that would promise preservation of shock entropy.
  initialGeometry = initialGeometry === undefined ? undefined : structuredClone(initialGeometry);
  const warm = initialFlow !== undefined;
  if (warm) {
    initialFlow = structuredClone(initialFlow);
    if (!initialFlow?.input || !initialFlow.initialEuler || !Array.isArray(initialFlow.initialEuler.nodes)
      || !(Array.isArray(initialFlow.initialEuler.x) || ArrayBuffer.isView(initialFlow.initialEuler.x)))
      throw new Error('initialFlow requires a source input and complete initialEuler x/nodes data.');
  }
  const targetSystem = createStreamtubeBodySystem(input), targetMach = targetSystem.conditions.mach;
  const targetGeometry = initialGeometry ? targetSystem.adoptGeometry(Float64Array.from(initialGeometry.x), initialGeometry.nodes)
    : targetSystem.initial.slice();
  // Reject a folded initial grid before changing the gas conditions. Mach
  // continuation is not a geometry-repair method.
  const initialQuality = warm ? null : streamtubeMeshSnapshot({ system: targetSystem, nodes: targetSystem.decode(targetGeometry).nodes }).quality;
  if (initialQuality && !initialQuality.valid) throw new Error('Automatic body startup requires a valid initial quadrilateral grid.');
  const attempts = []; let warmDiagnostics;
  const attempt = (prepare, label, mach, fraction, maxIterations) => {
    const history = []; let prepared, result, observerError;
    const observe = (callback, value) => {
      try { callback?.(value); } catch (error) { observerError = error; throw error; }
    };
    onStage?.({ label, mach, targetMach, fraction });
    try {
      prepared = prepare();
      result = solveStreamtubeBody(prepared.system, { initial: prepared.initial, maxIterations, tolerance,
        onIteration: h => {
          history.push({ ...h });
          observe(onIteration, { ...h, stage: label, mach, targetMach, continuationFraction: fraction });
        },
        onMesh: onMesh ? snapshot => observe(onMesh, {
          mesh: structuredClone(streamtubeMeshSnapshot(snapshot)), iteration: structuredClone(snapshot.iteration),
          stage: label, mach, targetMach, continuationFraction: fraction }) : undefined });
      const quality = streamtubeMeshSnapshot({ system: prepared.system, nodes: result.nodes }).quality;
      if (!quality.valid) result = { ...result, converged: false, reason: 'Final body grid fails the convexity/quality gate.' };
      attempts.push({ label, mach, fraction, converged: result.converged, reason: result.reason,
        iterations: history.length - 1, maximumResidual: result.diagnostics.residual,
        startup: prepared.diagnostics, history, gridValid: quality.valid });
      return { system: prepared.system, result };
    } catch (error) {
      // Cancellation and observer failures must reach the caller, rather
      // than silently causing more continuation stages.
      if (observerError) throw observerError;
      attempts.push({ label, mach, fraction, converged: false, reason: error.message,
        iterations: Math.max(0, history.length - 1), history, startup: prepared?.diagnostics,
        ...(error.code ? { code: error.code } : {}) });
      return null;
    }
  };
  const finish = (record, reachedTarget, reason) => ({ ...(record?.result ?? {}), system: record?.system ?? null,
    converged: reachedTarget && Boolean(record?.result.converged), ...(reason ? { reason } : {}),
    continuation: { method: 'freestream-mach', targetMach, currentMach: record?.system.conditions.mach ?? null,
      reachedTarget, attempts, densityReinitializedOnWarmRestart: false, analyticShockInformationUsed: false,
      ...(warmDiagnostics ? { initialFlow: warmDiagnostics } : {}) },
    physicalAcceptance: false, fullSolver: false });
  if (warm) {
    const source = createStreamtubeBodySystem(initialFlow.input);
    const encoded = Float64Array.from(initialFlow.initialEuler.x);
    if (encoded.length !== source.layout.n || !encoded.every(Number.isFinite))
      throw new Error('initialFlow requires the complete finite encoded source state.');
    const state = source.adoptGeometry(encoded, initialFlow.initialEuler.nodes);
    const sourceQuality = streamtubeMeshSnapshot({ system: source, nodes: source.decode(state).nodes }).quality;
    if (!sourceQuality.valid) throw new Error('Initial flow requires a valid convex quadrilateral grid.');
    const sourceFlow = source.evaluate(state);
    if (!(sourceFlow.diagnostics.residual <= tolerance))
      throw new Error('Initial flow source must already satisfy the residual tolerance.');
    // Check compatibility at the accepted source gas state, before any
    // target-state failure can trigger recovery. This rejects changed
    // topology, equations, epsilonP or upwind controls as invalid input.
    // In particular, no isentropic-to-hybrid reinterpretation is implicit.
    const sourceMach = source.conditions.mach;
    const compatible = initializeStreamtubeBodyFromFlow({ ...input, mach: sourceMach }, source, { initial: state });
    // The existing zero-update reporting path supplies the complete physical
    // result if all target attempts fail. It performs no Newton/LU work.
    const sourceResult = solveStreamtubeBody(compatible.system, { initial: compatible.initial, maxIterations: 0, tolerance });
    const compatibleQuality = streamtubeMeshSnapshot({ system: compatible.system, nodes: sourceResult.nodes }).quality;
    if (!sourceResult.converged || !compatibleQuality.valid)
      throw new Error('Compatible initial flow must remain converged on a valid convex quadrilateral grid.');
    warmDiagnostics = { sourceMach, sourceResidual: sourceResult.diagnostics.residual,
      sourceGridValid: true, unknowns: compatible.system.layout.n, completePhysicalState: true,
      coldInitializationSkipped: true, recoveryFractionStep: initialFractionStep,
      compatibility: compatible.diagnostics };
    let current = { system: compatible.system, result: sourceResult };
    const direct = attempt(() => initializeStreamtubeBodyFromFlow(input, current.system,
      { initial: current.result.x }), 'target-warm', targetMach, 1, directMaxIterations);
    if (direct?.result.converged) return finish(direct, true);
    // Retry from the last accepted complete flow. As in the cold driver's
    // continuation, initialFractionStep sets the first recovery increment;
    // failed increments halve and successful ones grow by at most 1.5.
    let fraction = 0, step = initialFractionStep, subdivisions = 0;
    while (fraction < 1) {
      if (attempts.length >= maxStages) return finish(current, false, 'Automatic warm Mach-continuation stage limit.');
      const next = Math.min(1, fraction + step), mach = next === 1 ? targetMach : sourceMach + next * (targetMach - sourceMach);
      const trial = attempt(() => initializeStreamtubeBodyFromFlow({ ...input, mach }, current.system,
        { initial: current.result.x }), 'warm-mach-continuation', mach, next, stageMaxIterations);
      if (trial?.result.converged) {
        current = trial; fraction = next; subdivisions = 0; step = Math.min(initialFractionStep, 1.5 * step);
      } else {
        if (subdivisions >= maxSubdivisions || step <= initialFractionStep / 2 ** maxSubdivisions)
          return finish(current, false, `Automatic warm Mach continuation stopped at Mach ${current.system.conditions.mach}.`);
        step *= .5; subdivisions++;
      }
    }
    return finish(current, true);
  }
  const cold = (mach, system = createStreamtubeBodySystem({ ...input, mach }), geometry) => {
    const seed = geometry ?? (initialGeometry ? system.adoptGeometry(Float64Array.from(initialGeometry.x), initialGeometry.nodes) : system.initial.slice());
    return { system, ...initializeStreamtubeStartup(system, seed) };
  };
  const direct = attempt(() => cold(targetMach, targetSystem, targetGeometry), 'target-direct', targetMach, 1, directMaxIterations);
  if (direct?.result.converged) return finish(direct, true);

  // At low Mach the fixed captured mass is compatible with a much greater
  // sonic capacity. Select an easier starting condition automatically;
  // the final target's masses, h0 and momentum equations are unchanged.
  let startMach = Math.min(.2, .5 * targetMach), current;
  for (let subdivision = 0; subdivision <= maxSubdivisions && attempts.length < maxStages; subdivision++) {
    current = attempt(() => cold(startMach), 'subcritical-start', startMach, 0, stageMaxIterations);
    if (current?.result.converged) break;
    startMach *= .5;
  }
  if (!current?.result.converged) return finish(direct ?? current, false, 'Automatic low-Mach body startup did not converge within its bounded attempts.');
  let fraction = 0, step = initialFractionStep, subdivisions = 0;
  while (fraction < 1) {
    if (attempts.length >= maxStages) return finish(current, false, 'Automatic Mach-continuation stage limit.');
    const next = Math.min(1, fraction + step), mach = next === 1 ? targetMach : startMach + next * (targetMach - startMach);
    // Only the last converged flow is used. Geometry and density jumps
    // survive transfer; target thermodynamics are recomputed, not relabeled
    // as an isentropic solution. Newton corrects the new target residual.
    const trial = attempt(() => initializeStreamtubeBodyFromFlow({ ...input, mach }, current.system,
      { initial: current.result.x }), 'mach-continuation', mach, next, stageMaxIterations);
    if (trial?.result.converged) {
      current = trial; fraction = next; subdivisions = 0; step = Math.min(initialFractionStep, 1.5 * step);
    } else {
      if (subdivisions >= maxSubdivisions || step <= initialFractionStep / 2 ** maxSubdivisions)
        return finish(current, false, `Automatic Mach continuation stopped at Mach ${current.system.conditions.mach}.`);
      step *= .5; subdivisions++;
    }
  }
  return finish(current, true);
}
