// SPDX-License-Identifier: GPL-2.0-or-later
import { observableFlow, observableBL } from './streamtube-flow-preview.js';
// Cold assembly path for the research simultaneous Euler/grid/BL/wake model.
// No saved solution or prescribed reference velocity is an input to this API.
import { planCoupledLogarithmicShearRecovery } from './streamtube-coupled-log-shear-recovery.js';
import { coupledFreshShearCoordinate, coupledResultShearCoordinate } from './streamtube-coupled-shear-policy.js';
import { coupledResolutionRecoveryCase, recoverCoupledResolution } from './streamtube-coupled-resolution-recovery.js';
import { coupledMachRecoveryCase, recoverCoupledMach } from './streamtube-coupled-mach-recovery.js';
import { retryCoupledEulerSeed } from './streamtube-coupled-seed-recovery.js';
import { solveStreamtubeAssembly, prepareStreamtubeAssembly } from './streamtube-result.js';
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { transferStreamtubeGeometry } from './streamtube-geometry.js';
import { initializeCoupledStreamtubeBody, panelCoupledEdgeGuess } from './streamtube-coupled-initializer.js';
import { solveCoupledStreamtubeIses as solveIses } from './streamtube-coupled-ises.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { initialStreamtubeDisplacement } from './streamtube-geometry.js';
import { transitionRecoveryPlan, recoverCoupledTransition, retainBestCoupledCheckpoint, stalledBoundaryLayerResolutionPlan } from './streamtube-transition-recovery.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { prepareCoupledGridSequence } from './streamtube-coupled-grid-sequence.js';
import { prepareCoupledGridLevels } from './streamtube-coupled-grid-levels.js';
import { solveCoupledGridLevel } from './streamtube-coupled-grid-levels.js';
import { continueCoupledNcrit } from './streamtube-coupled-ncrit-continuation.js';
import { coupledNcritStartupRecoveryPlan, finishCoupledNcritStartup } from './streamtube-coupled-ncrit-recovery.js';
import { streamtubeEquationControls } from './streamtube-equation-selection.js';
import { assertConvexStreamtubeGrid } from '../geometry/streamtube-convex-step.js';
import { capturePreparedStreamtubeAssembly, restorePreparedStreamtubeAssembly } from './streamtube-prepared-assembly.js';
import { selectCoupledEulerStartup, coupledStartupExtensionPlan, continueCoupledWakeCoordinates, finishSubsonicPressureCoupling } from './streamtube-coupled-startup.js';

// MSES manual §5.3 permits partial inviscid convergence before viscosity.
// This only recognizes a retained accepted endpoint at its iteration limit;
// the exact geometry/gas state and ordinary BL initializer are checked below.
function partialEulerStartup(precursor, maxIterations) {
  const flow = precursor?.flow, diagnostics = precursor?.diagnostics;
  const finiteVector = value => (Array.isArray(value) || ArrayBuffer.isView(value))
    && value.length > 0 && value.every(Number.isFinite);
  if (!(maxIterations > 0) || precursor?.status !== 'unconverged'
    || flow?.converged !== false || flow.residualConverged !== false
    || flow.reason !== 'iteration limit' || diagnostics?.reason !== 'iteration limit'
    || diagnostics.solverStopReason !== 'iteration limit' || flow.lastRejectedStep != null
    || flow.initialRedistribution?.accepted !== true
    || precursor.mesh?.quality?.valid !== true || flow.finalQuality?.valid !== true
    || !Array.isArray(precursor.mesh.quality.invalidCells) || precursor.mesh.quality.invalidCells.length
    || !Number.isInteger(diagnostics.iterations) || diagnostics.iterations < 1
    || flow.history?.at(-1)?.iteration !== diagnostics.iterations
    || !Number.isFinite(diagnostics.equationResidual) || diagnostics.equationResidual <= 0
    || !finiteVector(flow.x) || !finiteVector(flow.residual) || flow.x.length !== flow.residual.length)
    return null;
  return { startup: 'partial-inviscid', converged: false, reason: diagnostics.reason,
    iterations: diagnostics.iterations, residual: diagnostics.equationResidual };
}

// Sharp and finite-base profiles share complete-state grid sequencing.
// A fine request starts from converged coarse flow, then subdivides its fluid
// grid. Each level solves the same equations and reports its actual topology.
// A convex full-thickness initial guess can still be far from the coupled
// solution. Fine-grid sequencing must not depend on whether the geometric
// initializer happened to thin that guess.
export function coupledCoarseStartupCase(input, prepared, { maxIterations, enabled = true } = {}) {
  const bl = prepared?.system?.bl;
  if (!enabled || !(maxIterations > 0)
    || bl?.transitionMode !== 'automatic'
    || !bl.trips.every(pair => pair.every(v => v === 1))
    || !Number.isInteger(input.gridIntervals) || input.gridIntervals <= 16) return null;
  const tubes = input.gridTubes ?? 7;
  if (!Number.isInteger(tubes) || tubes < 1) return null;
  const coarse = { ...structuredClone(input), gridIntervals: 16, gridTubes: Math.min(tubes, 7) };
  // Nested refinement subdivides every streamwise interval, including the
  // inlet and wake. Explicit target counts must be coarsened too; copying
  // them unchanged multiplies a 128-interval farfield into 1024 intervals.
  const subdivisions = 2 ** Math.ceil(Math.log2(input.gridIntervals / 16));
  for (const key of ['gridInletIntervals', 'gridOutletIntervals'])
    if (Number.isInteger(input[key]) && input[key] >= 4 && input[key] <= 256)
      coarse[key] = Math.max(4, Math.ceil(input[key] / subdivisions));
  return coarse;
}

// Case trips are indexed by the user's elements. Streamtube topology sorts
// them by passage; map by the explicit element index, never by array position.
// Fractions refer to the initial contour-side parameter and become fixed
// material points. They are not chordwise x/c or a free-transition request.
// A warm checkpoint owns its Jacobian policy. An omitted old control means
// exact; an explicit request may verify that choice, never silently replace it.
export function coupledCheckpointHkPolicy(input, options) {
  if (input?.coupledNativeHk !== undefined && typeof input.coupledNativeHk !== 'boolean')
    throw new Error('Coupled Hk-floor linearization controls must be boolean.');
  if (!options || typeof options !== 'object' || Array.isArray(options)
    || ![undefined, 'exact', 'native'].includes(options.hkFloorLinearization))
    throw new Error('Invalid checkpoint Hk-floor linearization policy.');
  const policy = options.hkFloorLinearization ?? 'exact';
  if (input.coupledNativeHk !== undefined && (input.coupledNativeHk ? 'native' : 'exact') !== policy)
    throw new Error('The coupled case and checkpoint disagree on Hk-floor linearization.');
  return policy;
}

export function coupledAssemblyConditions(input, bodies, lengthScale, { checkpointOptions } = {}) {
  streamtubeEquationControls(input.eulerIsmom);
  if (input.coupledNativeHk !== undefined && typeof input.coupledNativeHk !== 'boolean')
    throw new Error('Coupled Hk-floor linearization controls must be boolean.');
  // Ordinary fresh starts and omitted legacy checkpoints retain exact
  // partials. A failed projected update may later select native sensitivities
  // at the same state, unless the caller explicitly requires exact partials.
  const nativeHk = checkpointOptions === undefined
    ? input.coupledNativeHk === true
    : coupledCheckpointHkPolicy(input, checkpointOptions) === 'native';
  const referenceChord = input.referenceChord ?? 1, referenceReynolds = input.reynolds ?? 1e6, ncrit = input.ncrit ?? 9;
  const transitionMode = input.transitionMode ?? 'fixed-trip', automatic = transitionMode === 'automatic';
  if (!['fixed-trip', 'automatic'].includes(transitionMode)) throw new Error('Unknown coupled transition mode.');
  const materialTrips = input.materialTrips ?? input.elements.map(() => automatic ? [1, 1] : [.05, .05]);
  if (![referenceChord, referenceReynolds, lengthScale, ncrit].every(Number.isFinite)
    || referenceChord <= 0 || referenceReynolds <= 0 || lengthScale <= 0 || ncrit <= 0)
    throw new Error('Coupled flow requires positive reference chord, Reynolds number, length scale and Ncrit.');
  if (!Array.isArray(materialTrips) || materialTrips.length !== input.elements.length || materialTrips.some(pair =>
    !Array.isArray(pair) || pair.length !== 2 || pair.some(v => !Number.isFinite(v) || v <= 0 || (automatic ? v > 1 : v >= 1))))
    throw new Error('Supply upper/lower material trip fractions strictly between zero and one for every element.');
  if (bodies.length !== input.elements.length || new Set(bodies.map(b => b.element)).size !== bodies.length
    || bodies.some(b => !Number.isInteger(b.element) || b.element < 0 || b.element >= input.elements.length))
    throw new Error('Coupled topology must identify each original element exactly once.');
  const reynolds = referenceReynolds * lengthScale / referenceChord;
  if (!Number.isFinite(reynolds) || reynolds <= 0) throw new Error('Invalid normalized coupled Reynolds number.');
  return { options: { reynolds, ncrit, tripFractions: bodies.map(b => materialTrips[b.element].slice()), edgeMatching: 'section-velocity',
      ...(automatic ? { transitionMode } : {}),
      ...(nativeHk ? { hkFloorLinearization: 'native' } : {}),
      ...(input.eulerIsmom === undefined ? {} : { blThermodynamics: 'historical-common-isentrope' }) },
    normalization: { referenceChord, referenceReynolds, solverLength: lengthScale, kernelReynolds: reynolds },
    materialTrips: materialTrips.map(pair => pair.slice()), elementOrder: bodies.map(b => b.element) };
}

// A refinement may fail before its first complete ISES checkpoint. Its
// retained provisional state still has the refined dimensions, so never
// split it using the earlier coarse initializer's Euler unknown count.
export function coupledAssemblyRestart(result) {
  if (result.checkpoint?.restart) return result.checkpoint.restart;
  const ne = result.x.length - 4 * result.boundaryLayer.stations.length;
  return { input: result.solverInput, options: result.coupledOptions,
    initialEuler: { x: result.x.slice(0, ne), nodes: result.flow.nodes, undisplacedNodes: result.flow.undisplacedNodes },
    initialBL: result.x.slice(ne) };
}

export function solveCoupledStreamtubeAssembly(input, { maxIterations = 20, eulerMaxIterations, tolerance = 1e-10,
  maxStartupAttempts = 2, coarseStartup = true, direct = false, resolutionRecovery = true, machRecovery = true, preparedEuler, onEulerPrepared,
  onMesh, onIteration, onStage, onCheckpoint, onIterationCheckpoint, onFlow, convergence = 'residual' } = {}) {
  if (!['residual', 'mses'].includes(convergence)) throw new Error('Invalid coupled convergence policy.');
  const solveCoupledStreamtubeIses = (input, options) => solveIses(input, { ...options, convergence });
  const solveAssembly = (input, options) => solveCoupledStreamtubeAssembly(input, { ...options, convergence });
  // Automatic untripped cases may start on a nested coarse grid at the same
  // operating point. The refined grid must pass its own convergence checks.
  // Direct solves retain one startup attempt and the requested transition controls.
  if (typeof direct !== 'boolean') throw new Error('Invalid direct solve control.');
  if (direct) maxStartupAttempts = 1;
  if (eulerMaxIterations === undefined) eulerMaxIterations = direct ? 40 : 20;
  if (!Array.isArray(input.elements) || !input.elements.length) throw new Error('Supply an airfoil assembly.');
  if (![maxIterations, eulerMaxIterations].every(n => Number.isInteger(n) && n >= 0)
    || ![1, 2].includes(maxStartupAttempts)
    || !Number.isFinite(tolerance) || tolerance <= 0 || typeof coarseStartup !== 'boolean' || typeof resolutionRecovery !== 'boolean' || typeof machRecovery !== 'boolean'
    || onEulerPrepared !== undefined && typeof onEulerPrepared !== 'function'
    || onIterationCheckpoint !== undefined && typeof onIterationCheckpoint !== 'function'
    || onFlow !== undefined && typeof onFlow !== 'function')
    throw new Error('Invalid coupled assembly iteration controls.');
  // Reject invalid physical controls before constructing a volume grid.
  const validated = coupledAssemblyConditions(input, input.elements.map((_, element) => ({ element })), input.referenceChord ?? 1);
  const earlyCoarseCase = coupledCoarseStartupCase(input, { system: { bl: {
    transitionMode: validated.options.transitionMode ?? 'fixed-trip', trips: validated.options.tripFractions,
  } } }, { maxIterations, enabled: coarseStartup });
  let stage = 'euler', startupAttempt = 0, displayedMethod;
  const changeStage = (value, details = {}) => {
    stage = value;
    onStage?.({ stage, startupAttempt,
      ...(direct ? { startupStrategy: 'direct-requested-conditions', mach: input.mach ?? .2,
        actualMach: input.mach ?? .2, targetMach: input.mach ?? .2,
        actualAlpha: input.alpha ?? 0, targetAlpha: input.alpha ?? 0, gridLevel: input.gridIntervals } : {}),
      ...(stage === 'boundary-layer-initialization' || stage.startsWith('coupled') ? displayedMethod : {}), ...details });
  };
  try {
    changeStage('euler');
    let precursor, partialEuler, solverInput, settings, initialEuler, actualMach, eulerInitialization;
    let requestedGeometry, requestedPacket = preparedEuler, fineEulerAttempted = false, fineEulerFailure;
    let prepared, result, restart, coarseInitialization, initialThicknessFactor = 1;
    let ncritInitialization, ncritRecoveryPlan, ncritRecovery, restorePreviousStartup, activeNcrit = validated.options.ncrit;
    const attempts = [];
    // The complete-state route needs target geometry/physics, not an unused
    // fine gas guess. Preserve the prepared packet for the ordinary fallback
    // and the outer Mach controller; neither capture nor restore solves gas.
    if (earlyCoarseCase) {
      requestedGeometry = preparedEuler === undefined ? prepareStreamtubeAssembly(input, {
        onMesh: (mesh, phase) => onMesh?.(mesh, phase, stage),
      }) : restorePreparedStreamtubeAssembly(preparedEuler, input);
      if (preparedEuler !== undefined) onMesh?.(requestedGeometry.mesh, 'initial', stage);
      requestedPacket = capturePreparedStreamtubeAssembly(requestedGeometry, input);
      onEulerPrepared?.(requestedPacket);
      solverInput = { ...requestedGeometry.input, wakeGeometry: 'independent-banks', wakeOutlet: 'banks',
        wakeDisplacementMotion: 'te-center' };
      settings = coupledAssemblyConditions(input, solverInput.bodies, requestedGeometry.system.conditions.lengthScale);
    }
    const ensureFineEuler = () => {
      if (fineEulerAttempted) {
        if (fineEulerFailure) throw fineEulerFailure.error;
        return;
      }
      fineEulerAttempted = true;
      try {
        if (stage !== 'euler') changeStage('euler');
        precursor = solveStreamtubeAssembly(input, { maxIterations: eulerMaxIterations, tolerance,
          // Use the same shock startup as standalone Euler. Both orders share
          // one precursor budget, on the requested grid and operating point.
          // The Euler initializer keeps smooth isentropic cases unchanged.
          firstOrderStartup: direct,
          retainBestCheckpoint: maxIterations > 0,
          ...(direct ? { stepAcceptance: 'armijo' } : {}),
          ...(requestedPacket === undefined ? {} : { preparedEuler: requestedPacket }),
          // The early path already published this same geometry handoff.
          ...(!requestedGeometry && onEulerPrepared !== undefined ? { onEulerPrepared } : {}),
          onMesh: (mesh, phase) => onMesh?.(mesh, phase, stage),
          onIteration: h => onIteration?.({ ...h, stage: 'euler' }) });
        onCheckpoint?.({ stage: 'euler', result: precursor });
        partialEuler = partialEulerStartup(precursor, maxIterations);
        if (precursor.status !== 'research-converged' && !partialEuler)
          throw Object.assign(new Error(`Coupled Euler precursor did not converge: ${precursor.diagnostics.reason}`),
            { code: 'coupled-euler-precursor', initialization: { stage, result: precursor } });
        const selected = partialEuler ? selectCoupledEulerStartup(precursor, { tolerance }) : null;
        const seed = selected ?? precursor;
        const source = createStreamtubeBodySystem(seed.solverInput);
        const sourceState = source.adoptGeometry(seed.flow.x, seed.flow.nodes);
        // Metadata alone cannot certify a partial state. Retain the exact
        // existing gas replay and convex-corner check before BL transfer.
        if (partialEuler) assertConvexStreamtubeGrid(source.evaluate(sourceState).nodes);
        solverInput = { ...seed.solverInput, wakeGeometry: 'independent-banks', wakeOutlet: 'banks',
          // A budget-limited precursor can still use temporary broadening.
          // Transfer its physical state, but retain the requested viscous law,
          // including second order if the inviscid budget ended in startup.
          ...(precursor.flow.adaptiveMcrit || precursor.flow.firstOrderStartup ? { upwind: { ...seed.solverInput.upwind,
            ...(precursor.flow.adaptiveMcrit ? { mcrit: precursor.flow.targetMcrit } : {}),
            ...(precursor.flow.firstOrderStartup ? { mucon: precursor.flow.targetMucon } : {}) } } : {}),
          wakeDisplacementMotion: 'te-center' };
        const displacement = initialStreamtubeDisplacement(source.layout, source.baseGeometry);
        const target = createStreamtubeBodySystem({ ...solverInput, displacement });
        const eulerState = transferStreamtubeGeometry(source, sourceState, target);
        settings = coupledAssemblyConditions(input, solverInput.bodies, source.conditions.lengthScale);
        // Start low-Mach, fully subsonic cases with the final pressure
        // matching equations. An intermediate section-speed solve can pinch
        // the stagnation inlet before the pressure handoff. The choice
        // depends on the converged flow, not how its grid was recovered.
        // Keep transonic runs and supplied checkpoint formulations intact.
        if (direct && !partialEuler && (input.mach ?? .2) <= .3
          && seed.flow?.diagnostics?.maxMach < 1) {
          settings.options.edgeMatching = 'pressure';
          delete settings.options.blThermodynamics;
        }
        initialEuler = { x: eulerState, ...target.decode(eulerState) };
        actualMach = precursor.mach;
        eulerInitialization = { iterations: precursor.diagnostics.iterations, residual: precursor.diagnostics.equationResidual,
          cells: precursor.diagnostics.cells, gridSmoothing: precursor.mesh.initialization.gridSmoothing,
          gasInitialization: precursor.mesh.initialization.gasInitialization,
          ...(precursor.flow.passageRecovery ? { passageRecovery: precursor.flow.passageRecovery } : {}),
          ...(partialEuler ? partialEuler : {}),
          ...(selected ? { selectedState: selected.selection } : {}) };
        if (selected) {
          const mesh = streamtubeMeshSnapshot({ system: source, nodes: seed.flow.nodes, flow: seed.flow,
            iteration: { ...seed.flow.history.at(-1), iteration: selected.selection.selectedIteration } });
          mesh.initialization = { ...precursor.mesh.initialization, ...mesh.initialization, eulerStartup: selected.selection };
          onMesh?.(mesh, 'initial', stage);
          onCheckpoint?.({ stage: 'euler-startup-selection', selection: selected.selection, result: seed });
        }
      } catch (error) {
        fineEulerFailure = { error };
        if (coarseInitialization && error && (typeof error === 'object' || typeof error === 'function') && Object.isExtensible(error))
          error.initialization = { ...error.initialization, coarseInitialization: structuredClone(coarseInitialization) };
        throw error;
      }
    };
    if (!earlyCoarseCase) ensureFineEuler();
    const ncritLabels = () => ({ ...displayedMethod,
      ...(ncritInitialization || ncritRecoveryPlan ? {
        actualNcrit: activeNcrit, targetNcrit: validated.options.ncrit,
        ncritContinuation: { actualNcrit: activeNcrit, targetNcrit: validated.options.ncrit,
          reachedTarget: activeNcrit === validated.options.ncrit } } : {}) });
    const initializeDirect = () => {
      ensureFineEuler();
      if (stage !== 'boundary-layer-initialization') changeStage('boundary-layer-initialization', {
        initialThicknessFactor, ...ncritLabels(), ...(coarseInitialization ? { coarseInitialization } : {}),
      });
      let panelGuess;
      // Panel/BL supplies only a startup guess. It also avoids a poor cold
      // BL march through the precursor's developing shock: the simultaneous
      // solve still uses the requested Mach, thermodynamics and edge equations.
      // The ordinary geometric/domain and resolved-transition checks below
      // decide whether its edge-speed and complete profile transfers are usable.
      if (settings.options.edgeMatching === 'pressure' || direct) {
        try { panelGuess = panelCoupledEdgeGuess(input, solverInput.bodies, settings.normalization.solverLength); }
        catch (error) { panelGuess = { accepted: false, reason: error.message }; }
      }
      const initializationOptions = { ...settings.options, initialEuler,
        ...(panelGuess?.accepted ? { initialEdgeVelocity: panelGuess.initialEdgeVelocity } : {}),
        ...(ncritRecoveryPlan ? { ncrit: activeNcrit,
          ...(ncritRecoveryPlan.hkFloorLinearization === 'native' ? { hkFloorLinearization: 'native' } : {}) } : {}) };
      let initialized = initializeCoupledStreamtubeBody(solverInput, initializationOptions, { initialThicknessFactor });
      let panelProfile;
      if (panelGuess?.accepted && initialized.system.bl.transitionMode === 'automatic') {
        try {
          const { profileTransfer, ...mapped } = panelGuess.initialBoundaryLayer(initialized.system,
            { thicknessFactor: initialized.initialization.thicknessFactor });
          // Rebuild from the undisplaced precursor, not the already extended
          // interiors of the edge-speed seed. Otherwise thickness is added twice.
          const candidate = initializeCoupledStreamtubeBody(solverInput,
            { ...initializationOptions, ...mapped }, { maximumBacktracks: 0, initialThicknessFactor });
          initialized = candidate;
          panelProfile = { accepted: true, fields: ['theta', 'deltaStar', 'ue', 'shear', 'transition'],
            ...(profileTransfer ? { transfer: profileTransfer } : {}) };
        } catch (error) { panelProfile = { accepted: false, reason: error.message }; }
      }
      if (panelGuess) initialized.initialization.panelEdgeGuess = { accepted: panelGuess.accepted,
        ...(panelGuess.accepted ? { residual: panelGuess.residual } : { reason: panelGuess.reason }),
        ...(panelProfile ? { profile: panelProfile } : {}),
        initialGuessOnly: true };
      return initialized;
    };
    for (startupAttempt = 1; startupAttempt <= maxStartupAttempts; startupAttempt++) {
      let activeSolverInput, gridSequence;
      displayedMethod = { shearCoordinate: ncritRecoveryPlan?.shearCoordinate ?? 'linear',
        hkFloorLinearization: ncritRecoveryPlan?.hkFloorLinearization ?? settings?.options.hkFloorLinearization ?? 'exact' };
      changeStage('boundary-layer-initialization', { initialThicknessFactor, ...ncritLabels(),
        ...(partialEuler ? { eulerPreparation: { ...partialEuler } } : {}),
        ...(attempts.length ? { retryReason: attempts.at(-1).reason } : {}) });
      // A complete-state sequence does not use a cold fine BL profile. Do
      // not require that unused guess to succeed before solving its parent.
      // Keep the ordinary initializer as the fallback if sequencing fails.
      const deferredCoarseCase = startupAttempt === 1 ? earlyCoarseCase : null;
      try { prepared = deferredCoarseCase ? null : initializeDirect(); }
      catch (error) {
        if (!ncritRecoveryPlan || !result?.checkpoint || !restorePreviousStartup) throw error;
        const failure = { message: error?.message ?? String(error), code: error?.code, diagnostics: error?.diagnostics };
        const attemptedNcrit = activeNcrit;
        activeNcrit = result.checkpoint.restart.options.ncrit;
        ncritRecovery = { ...ncritRecoveryPlan, attempted: true, actualNcrit: activeNcrit,
          reachedTarget: false, stateConverged: false, failure, failedStage: 'boundary-layer-initialization' };
        result = { ...result, actualNcrit: activeNcrit, targetNcrit: validated.options.ncrit,
          reason: `Lower-Ncrit startup initialization failed: ${failure.message}. Retained the previous admissible state at Ncrit ${activeNcrit}.`,
          ncritContinuation: { startup: ncritRecovery, actualNcrit: activeNcrit, targetNcrit: validated.options.ncrit,
            reachedTarget: activeNcrit === validated.options.ncrit, stateConverged: false } };
        attempts.push({ startupAttempt, actualNcrit: attemptedNcrit, thicknessFactor: initialThicknessFactor,
          converged: false, reason: failure.message, initializationFailed: true, failure, retainedNcrit: activeNcrit });
        restorePreviousStartup();
        break;
      }
      activeSolverInput = solverInput;
      const coarseCase = deferredCoarseCase ?? (startupAttempt === 1
        ? coupledCoarseStartupCase(input, prepared, { maxIterations, enabled: coarseStartup }) : null);
      if (coarseCase) {
        const ncritStartup = !direct && deferredCoarseCase && requestedGeometry.system.baseGeometry?.some(Boolean)
          && validated.options.ncrit > 4;
        if (ncritStartup) {
          coarseCase.ncrit = 4;
          coarseCase.coupledNativeHk = input.coupledNativeHk !== false;
        }
        const fineInitial = prepared;
        coarseInitialization = { attempted: true, accepted: false,
          requestedGrid: { intervals: input.gridIntervals, tubes: input.gridTubes ?? 7 },
          sourceGrid: { intervals: coarseCase.gridIntervals, tubes: coarseCase.gridTubes },
          ...(fineInitial ? { originalThicknessFactor: fineInitial.initialization.thicknessFactor }
            : { directFineBLInitialization: 'deferred', directFineEulerInitialization: 'deferred' }),
          equationsChanged: false, requestedGridRetained: true };
        let observerFailed = false, observerError, retainedCoarse, retainedIntervals = coarseCase.gridIntervals, levelDiagnostics;
        const observe = callback => (...args) => {
          try { return callback?.(...args); }
          catch (error) { observerFailed = true; observerError = error; throw error; }
        };
        try {
          observe(changeStage)('coupled-coarse-initialization', { coarseStage: 'euler', ...coarseInitialization,
            gridLevel: coarseCase.gridIntervals, requestedGridIntervals: input.gridIntervals });
          let coarse = solveCoupledStreamtubeAssembly(coarseCase, { maxIterations, eulerMaxIterations, tolerance, convergence,
            maxStartupAttempts, coarseStartup: false, direct,
            onStage: observe(event => changeStage('coupled-coarse-initialization', {
              coarseStage: event.stage, coarseStartupAttempt: event.startupAttempt,
              gridLevel: coarseCase.gridIntervals, requestedGridIntervals: input.gridIntervals,
              ...(ncritStartup ? { actualNcrit: event.actualNcrit ?? coarseCase.ncrit, targetNcrit: validated.options.ncrit } : {}) })),
            onIteration: observe(h => onIteration?.({ ...h, coarseStage: h.stage,
              coarseStartupAttempt: h.startupAttempt, startupAttempt, stage: 'coupled-coarse-initialization',
              ...(ncritStartup ? { actualNcrit: h.actualNcrit ?? coarseCase.ncrit, targetNcrit: validated.options.ncrit } : {}) })),
            ...(onFlow ? { onFlow: observe(frame => onFlow({ ...frame, coarseStage: frame.stage,
              coarseStartupAttempt: frame.startupAttempt, startupAttempt, stage: 'coupled-coarse-initialization',
              ...(ncritStartup ? { actualNcrit: frame.checkpoint.restart.options.ncrit, targetNcrit: validated.options.ncrit } : {}) })) } : {}),
            onMesh: observe((mesh, phase) => onMesh?.(mesh, phase, 'coupled-coarse-initialization')) });
          if (observerFailed) throw observerError;
          coarseInitialization.source = { converged: coarse.converged, reason: coarse.reason,
            families: coarse.families, cells: coarse.mesh.quality.valid ? coarse.mesh.cells.length : undefined,
            quality: coarse.mesh.quality,
            iterations: coarse.initialization.attempts.reduce((sum, attempt) => sum + attempt.iterations, 0),
            eulerIterations: coarse.initialization.euler.iterations };
          if (!coarse.converged) throw new Error(`Coarse BL startup did not converge: ${coarse.reason}`);
          if (ncritStartup) {
            const continued = continueCoupledNcrit(coarse, { targetNcrit: validated.options.ncrit, phase: 'coarse',
              maxIterations, tolerance, normalization: settings.normalization, startupAttempt,
              onStage: observe(event => changeStage(event.stage, event)), onIteration: observe(onIteration),
              onMesh: observe(onMesh), onFlow: observe(onFlow), onIterationCheckpoint: observe(onIterationCheckpoint) });
            coarse = continued.result;
            ncritInitialization = { coarse: continued.diagnostics };
            activeNcrit = coarse.checkpoint.restart.options.ncrit;
          }
          retainedCoarse = coarse;
          const source = coarse.checkpoint.restart;
          // A deferred fine precursor has not selected its edge formulation.
          // Retain the converged parent's choice on every refinement level.
          if (direct) {
            settings.options.edgeMatching = source.options.edgeMatching;
            if (source.options.geometryReplay !== undefined) settings.options.geometryReplay = source.options.geometryReplay;
            if (source.options.hkFloorLinearization !== undefined) settings.options.hkFloorLinearization = source.options.hkFloorLinearization;
            if (source.options.blThermodynamics === undefined) delete settings.options.blThermodynamics;
            else settings.options.blThermodynamics = source.options.blThermodynamics;
          }
          const sourceSystem = createCoupledStreamtubeBody(source.input, { ...source.options,
            initialEuler: source.initialEuler, initialBL: source.initialBL });
          {
            let parentUnknowns = sourceSystem.n;
            // Refinement can move natural transition across a new station.
            // Rebuild sharp-profile BL guesses with the native mixed profile
            // initializer so their thickness and N/Ctau phases agree. The
            // qualified finite-base path retains its complete interpolation.
            const blPredictor = sourceSystem.bl.hasFiniteBase ? 'interpolate' : 'xfoil-mrchdu';
            const sequenceArgs = { sourceSystem, sourceResult: coarse,
              input: solverInput, options: { ...settings.options, ncrit: activeNcrit,
                ...(source.options.hkFloorLinearization === 'native' ? { hkFloorLinearization: 'native' } : {}) },
              requestedGridIntervals: input.gridIntervals, sourceGridIntervals: coarseCase.gridIntervals };
            // Keep the qualified 16→32 path unchanged. Larger requests solve
            // intermediate levels instead of taking an eightfold cold jump.
            if (input.gridIntervals <= 2 * coarseCase.gridIntervals)
              prepared = prepareCoupledGridSequence(sequenceArgs, { tolerance, blPredictor });
            else {
              const levels = prepareCoupledGridLevels(sequenceArgs, { tolerance, blPredictor,
                solveIntermediateLevel: (seed, level) => solveCoupledGridLevel(seed, {
                  ...level, targetGridIntervals: input.gridIntervals, targetNcrit: validated.options.ncrit,
                }, { maxIterations, tolerance, normalization: settings.normalization, startupAttempt,
                  onStage: observe(event => changeStage(event.stage, event)), onIteration: observe(onIteration),
                  onMesh: observe(onMesh), onFlow: onFlow ? observe(onFlow) : undefined,
                  onIterationCheckpoint: observe(onIterationCheckpoint) }) });
              levelDiagnostics = levels.diagnostics;
              retainedCoarse = { ...coarse, ...levels.sourceResult,
                // Intermediate results describe their BL/grid transfer, not a
                // new Euler precursor. Preserve the original startup record.
                initialization: { ...coarse.initialization, ...levels.sourceResult?.initialization } };
              retainedIntervals = levels.sourceGridIntervals;
              parentUnknowns = levels.sourceSystem.n;
              prepared = levels.prepared;
              if (!prepared) throw new Error(levels.diagnostics.failure?.reason ?? 'An intermediate grid level did not converge.');
            }
            activeSolverInput = prepared.input;
            actualMach = coarse.mach;
            eulerInitialization = { ...structuredClone(coarse.initialization.euler), source: 'coarse-grid',
              sourceGridIntervals: coarseCase.gridIntervals, requestedFineEulerPerformed: false };
            coarseInitialization.directFineEulerInitialization = 'skipped';
            gridSequence = { kind: 'coarse-to-fine', parentUnknowns, unknowns: prepared.system.n,
              requestedGridIntervals: input.gridIntervals, actualGridIntervals: prepared.transfer.nominalGridIntervals,
              initialization: prepared.transfer, ...(levelDiagnostics ? { levels: levelDiagnostics } : {}) };
            prepared.mesh.initialization.gridRefinement = gridSequence;
            coarseInitialization.requestedGridRetained = false;
            coarseInitialization.actualGrid = { cells: prepared.mesh.cells.length,
              intervals: prepared.system.euler.layout.nx, tubes: prepared.system.euler.layout.tubes.slice() };
          }
          coarseInitialization.accepted = true;
          coarseInitialization.mapping = prepared.transfer ?? prepared.initialization?.transfer;
        } catch (error) {
          if (observerFailed) throw observerError;
          coarseInitialization.reason = error.message;
          if (retainedCoarse) {
            // An accepted parent is useful evidence and a restart. Never
            // discard it to retry the unused fine cold guess after a failed
            // refinement, or present its coefficients as a fine-grid root.
            const retained = retainedCoarse, f = retained.checkpoint.restart;
            const actualNcrit = f.options.ncrit;
            const sequence = { kind: 'coarse-to-fine', reachedTarget: false,
              requestedGridIntervals: input.gridIntervals, actualGridIntervals: retainedIntervals,
              ...(levelDiagnostics ? { levels: levelDiagnostics } : {}), failure: { reason: error.message } };
            const reason = `Requested grid (${input.gridIntervals} surface intervals) was not reached; retained converged ${retainedIntervals}-interval level. ${error.message}`;
            observe(changeStage)('coupled-grid-refinement', { gridLevel: retainedIntervals,
              requestedGridIntervals: input.gridIntervals, retained: true, reason });
            observe(onMesh)(retained.mesh, 'solving', stage);
            observe(onCheckpoint)({ stage, startupAttempt, families: retained.families, restart: f });
            if (onFlow) {
              const system = createCoupledStreamtubeBody(f.input, { ...f.options,
                initialEuler: f.initialEuler, initialBL: f.initialBL });
              const value = system.evaluate(system.initial);
              observe(onFlow)(structuredClone({ checkpoint: retained.checkpoint,
                flow: observableFlow(value.outer),
                bl: observableBL(system.bl),
                bodies: system.euler.layout.bodies, normalization: settings.normalization,
                iteration: retained.history.at(-1), mach: f.input.mach, actualNcrit, targetNcrit: validated.options.ncrit,
                stage, startupAttempt, gridLevel: retainedIntervals, requestedGridIntervals: input.gridIntervals, retained: true }));
            }
            return { ...retained, converged: false, stateConverged: true, status: 'unconverged', reason,
              restart: f, mach: f.input.mach, actualNcrit, targetNcrit: validated.options.ncrit,
              ...(actualNcrit !== validated.options.ncrit ? { ncritContinuation: {
                ...ncritInitialization, actualNcrit, targetNcrit: validated.options.ncrit,
                reachedTarget: false, stateConverged: true } } : {}),
              gridSequence: sequence, automaticRefinement: sequence,
              initialization: { ...retained.initialization, coarseInitialization: {
                ...coarseInitialization, requestedGridRetained: false,
                actualGrid: { cells: retained.mesh.cells.length, intervals: f.input.outerLower?.length - 1,
                  tubes: f.input.weights?.map(row => row.length) } } } };
          }
          prepared = fineInitial ?? initializeDirect();
          if (!fineInitial) {
            coarseInitialization.directFineBLInitialization = 'fallback';
            coarseInitialization.directFineEulerInitialization = 'fallback';
            coarseInitialization.originalThicknessFactor = prepared.initialization.thicknessFactor;
          }
          activeSolverInput = solverInput; gridSequence = undefined;
          activeNcrit = validated.options.ncrit; ncritInitialization = undefined;
          coarseInitialization.requestedGridRetained = true;
          delete coarseInitialization.actualGrid;
        }
        changeStage('boundary-layer-initialization', { coarseInitialization });
        // Restore the requested topology before any subsequent observer can
        // cancel. A failed coarse attempt must not leave its preview active.
        onMesh?.(prepared.mesh, 'initial', stage);
      }
      const s = prepared.system, seed = s.evaluate(s.initial);
      const stateOptions = { ...settings.options, ncrit: activeNcrit,
        ...(s.conditions.hkFloorLinearization === 'native' ? { hkFloorLinearization: 'native' } : {}),
        ...(s.bl.transitionMode === 'automatic' ? { transitionState: s.bl.snapshotActive() } : {}) };
      const checkpoint = { input: activeSolverInput, options: stateOptions,
        initialEuler: { x: s.initial.slice(0, s.ne), nodes: seed.outer.nodes, undisplacedNodes: seed.outer.undisplacedNodes },
        initialBL: s.initial.slice(s.ne) };
      displayedMethod = { shearCoordinate: coupledFreshShearCoordinate(prepared, ncritRecoveryPlan),
        hkFloorLinearization: stateOptions.hkFloorLinearization ?? 'exact' };
      changeStage('boundary-layer-initialization', { ...ncritLabels(), boundaryLayerInitialization: {
        method: prepared.initialization?.method, wakeInitialization: prepared.initialization?.wakeInitialization,
        ...(prepared.initialization?.panelEdgeGuess?.profile
          ? { panelProfile: prepared.initialization.panelEdgeGuess.profile } : {}),
        thicknessFactor: prepared.initialization?.thicknessFactor, accepted: true } });
      onCheckpoint?.({ stage, startupAttempt, ...ncritLabels(), initialization: prepared.initialization, families: seed.families, restart: checkpoint });
      if (!coarseCase) onMesh?.(prepared.mesh, 'initial', stage);
      changeStage('coupled', ncritLabels());
      // Check for a repeated transition cycle after one ordinary startup
      // chunk. If none occurs, use the caller's remaining iteration budget
      // on the same state. A single local refinement replaces the thinner
      // retry when it succeeds; it never changes the requested flow data.
      const automatic = s.bl.transitionMode === 'automatic';
      const subsonicStartup = direct && automatic && maxIterations >= 2 && (input.mach ?? 0) <= .3
        && stateOptions.edgeMatching === 'section-velocity' && seed.outer.diagnostics.maxMach < 1
        && prepared.initialization?.panelEdgeGuess?.accepted === true;
      // Reserve most of the budget for the final pressure formulation. A
      // prolonged section-speed startup can pinch the stagnation inlet cell
      // before that formulation has a chance to regularize the edge state.
      const firstBudget = subsonicStartup ? Math.min(12, Math.floor(maxIterations / 2))
        : automatic && !direct ? Math.min(maxIterations, 20) : maxIterations;
      const iterationCallback = h => {
        if (h.hkFloorLinearization) displayedMethod = { ...displayedMethod, hkFloorLinearization: h.hkFloorLinearization };
        onIteration?.({ ...h, stage: 'coupled', startupAttempt, ...ncritLabels() });
      };
      // ISES publishes its complete checkpoint immediately before the same
      // accepted flow's mesh callback. Never borrow a prior startup's state.
      let flowCheckpoint, bestState, checkpointIterationOffset = 0;
      const checkpointObserver = automatic || onIterationCheckpoint || onFlow ? { onCheckpoint: (value, details) => {
        if (automatic && startupAttempt === 1) bestState = retainBestCoupledCheckpoint(bestState, value, {
          iteration: (details.history?.at(-1)?.iteration ?? -1) + checkpointIterationOffset,
          initialRedistribution: details.initialRedistribution,
        });
        if (onFlow) flowCheckpoint = value;
        onIterationCheckpoint?.(structuredClone(value), structuredClone({ ...details, stage: 'coupled', startupAttempt, ...ncritLabels(),
          ...(checkpointIterationOffset ? { iterationOffset: checkpointIterationOffset } : {}) }));
      } } : {};
      const meshCallback = state => {
        const mesh = streamtubeMeshSnapshot({ ...state, iteration: { ...state.iteration, startupAttempt } });
        if (gridSequence) mesh.initialization.gridRefinement = gridSequence;
        onMesh?.(mesh,
          state.iteration.iteration === 0 ? 'initial' : 'solving', stage);
        if (onFlow && flowCheckpoint) onFlow(structuredClone({ checkpoint: flowCheckpoint,
          flow: observableFlow(state.flow),
          // Only invariant station topology is needed; phase lives in the
          // matching checkpoint, not the initializer's mutable metadata.
          bl: observableBL(s.bl), bodies: state.system.layout.bodies,
          normalization: settings.normalization, iteration: state.iteration,
          mach: flowCheckpoint.restart.input.mach, stage, startupAttempt, ...ncritLabels() }));
      };
      restorePreviousStartup = () => {
        changeStage('coupled', { ...ncritLabels(), ncritContinuation: result.ncritContinuation, retained: true });
        flowCheckpoint = result.checkpoint;
        restart = result.checkpoint.restart;
        onIterationCheckpoint?.(structuredClone(result.checkpoint), { stage, startupAttempt,
          ...ncritLabels(), ncritContinuation: result.ncritContinuation, retained: true });
        const f = result.checkpoint.restart;
        const retainedSystem = createCoupledStreamtubeBody(f.input, { ...f.options,
          initialEuler: f.initialEuler, initialBL: f.initialBL });
        iterationCallback({ ...result.history.at(-1), retained: true });
        meshCallback({ system: retainedSystem.euler, flow: result.flow, nodes: result.flow.nodes,
          iteration: { ...result.history.at(-1), retained: true } });
        onCheckpoint?.({ stage, startupAttempt, ...ncritLabels(), ncritContinuation: result.ncritContinuation,
          families: result.families, restart, retained: true });
      };
      result = solveCoupledStreamtubeIses(activeSolverInput, { ...stateOptions,
        // Preserve the section-speed startup's accepted-step sequence. Its
        // event boundaries can amplify roundoff from a different LU order.
        ...(direct && stateOptions.edgeMatching === 'pressure' ? { linearOrdering: 'aligned-auto' } : {}),
        initialEuler: checkpoint.initialEuler, initialBL: checkpoint.initialBL, maxIterations: firstBudget, tolerance,
        // Only this cold inviscid-to-viscous handoff requests a new wake seed.
        // Direct library solves and checkpoint resumes preserve their chart.
        wakeGridInitialization: direct && automatic && maxIterations > 0,
        // Direct Newton already solves the interior grid with BL displacement.
        // An extra interior extension after DSLIM can turn its admissible
        // trial into a negative-pressure state. Keep the corrected wall and
        // Newton interior coordinates; all final domain/merit gates still run.
        // Legacy continuation paths and resumed chunks retain their policy.
        iterationGeometry: 'ises-sampled', stepAcceptance: direct ? 'event-armijo' : 'admissible',
        ...(automatic ? { blUpdate: 'xfoil', projectionGeometry: direct
          && !(solverInput.stagnationMotion === 'interpolated' && settings.options.edgeMatching === 'section-velocity')
          ? 'fixed' : 'boundary-increment' } : {}),
        ...(input.coupledNativeHk === false ? { hkProjectionRecovery: false } : {}),
        shearCoordinate: coupledFreshShearCoordinate(prepared, ncritRecoveryPlan),
        ...checkpointObserver, onIteration: iterationCallback, onMesh: meshCallback });
      if (subsonicStartup) {
        result = finishSubsonicPressureCoupling(result, { maxIterations, tolerance, linearOrdering: 'aligned-auto',
          ...checkpointObserver, onIteration: iterationCallback, onMesh: meshCallback,
          onCorrection: couplingStartup => changeStage('coupled', { ...ncritLabels(), couplingStartup }) });
        if (result.couplingStartup) {
          settings.options.edgeMatching = 'pressure';
          delete settings.options.blThermodynamics;
        }
      }
      if (direct) result = continueCoupledWakeCoordinates(result, { maxIterations, tolerance,
        ...checkpointObserver, onIteration: iterationCallback, onMesh: meshCallback,
        onCorrection: wakeCorrection => changeStage('coupled', { ...ncritLabels(), wakeCorrection }) });
      let coarseHistory, transitionRecoveryFailure;
      const recoveryPlan = !direct && maxIterations > 0 && startupAttempt === 1 ? transitionRecoveryPlan(result) : null;
      if (recoveryPlan) {
        coarseHistory = result.history;
        result.mesh.initialization.gridSmoothing = eulerInitialization.gridSmoothing;
        let recovered, observerFailed = false;
        const observe = callback => (...args) => {
          try { return callback?.(...args); }
          catch (error) { observerFailed = true; throw error; }
        };
        try {
          recovered = recoverCoupledTransition(result, { plan: recoveryPlan, maxIterations, tolerance,
            onIteration: observe(h => iterationCallback({ ...h, iteration: h.iteration + offset, transitionRecovery: true })),
            onMesh: observe(onMesh), onStage: observe(detail => changeStage(detail.stage, detail)),
            onIterationCheckpoint: observe((cp, details) => onIterationCheckpoint?.(cp, { ...details, iterationOffset: offset })),
            onFlow: onFlow ? observe(flow => onFlow({ ...flow, iteration: { ...flow.iteration, iteration: flow.iteration.iteration + offset } })) : undefined,
            normalization: settings.normalization, startupAttempt });
        } catch (error) {
          if (observerFailed) throw error;
          // Retain the admissible coarse state if a finer initial guess is
          // impossible. This is a failed recovery, not a certified result.
          transitionRecoveryFailure = error.message;
          // A preview may already have published the finer topology. Restore
          // the retained mesh too, so a later failure cannot overlay that
          // abandoned flow snapshot on the coarse returned state.
          changeStage('coupled', { transitionRecoveryFailure });
          onMesh?.(result.mesh, 'solving', stage);
        }
        if (recovered) result = recovered;
      } else if (!result.converged && (result.solverStopReason ?? result.reason) === 'iteration limit'
        && result.checkpoint && !result.couplingStartup && firstBudget < maxIterations) {
        const before = result;
        checkpointIterationOffset = firstBudget;
        result = solveCoupledStreamtubeIses(undefined, { resume: before.checkpoint,
          maxIterations: maxIterations - firstBudget, tolerance, iterationGeometry: 'ises-sampled', stepAcceptance: direct ? 'event-armijo' : 'admissible',
          ...checkpointObserver,
          onIteration: h => { if (h.iteration) iterationCallback({ ...h, iteration: h.iteration + firstBudget }); },
          onMesh: state => { if (state.iteration.iteration) meshCallback({ ...state,
            iteration: { ...state.iteration, iteration: state.iteration.iteration + firstBudget } }); } });
        result.history = [...before.history, ...result.history.slice(1).map(h => ({ ...h, iteration: h.iteration + firstBudget }))];
        const a = before.linearDiagnostics, b = result.linearDiagnostics;
        result.linearDiagnostics = { solves: a.solves + b.solves, refinements: a.refinements + b.refinements,
          pivotRecoveries: a.pivotRecoveries + b.pivotRecoveries,
          maxRelativeResidual: Math.max(a.maxRelativeResidual, b.maxRelativeResidual),
          iterations: [...a.iterations, ...b.iterations.map(h => ({ ...h, iteration: h.iteration + firstBudget }))] };
      }
      // Continue a productive retained state before discarding it for a new
      // startup. Every chunk must earn its budget again; a finite reserve
      // caps the total extra work even when transition changes recur.
      let extensionCount = 0;
      for (;;) {
        const extensionPlan = !recoveryPlan && !ncritInitialization && !ncritRecoveryPlan
          ? coupledStartupExtensionPlan(result, { startupAttempt, maxIterations,
            extensionCount, tolerance, coarseInitialization,
            requestedConditions: { mach: input.mach ?? 0, reynolds: settings.options.reynolds,
              ncrit: validated.options.ncrit, transitionMode: settings.options.transitionMode,
              tripFractions: settings.options.tripFractions, ismom: input.eulerIsmom } }) : null;
        if (!extensionPlan) break;
        const before = result, offset = extensionPlan.originalIterations, controls = before.checkpoint.continuation;
        checkpointIterationOffset = offset;
        changeStage('coupled', { startupExtension: { kind: extensionPlan.kind,
          originalIterations: offset, additionalIterations: extensionPlan.additionalIterations } });
        result = solveCoupledStreamtubeIses(undefined, { resume: extensionPlan.resume,
          maxIterations: extensionPlan.additionalIterations, tolerance,
          iterationGeometry: controls.iterationGeometry, stepAcceptance: controls.stepAcceptance,
          stagnationLimiter: controls.stagnationLimiter, ...checkpointObserver,
          onIteration: h => { if (h.iteration) iterationCallback({ ...h, iteration: h.iteration + offset }); },
          onMesh: state => { if (state.iteration.iteration) meshCallback({ ...state,
            iteration: { ...state.iteration, iteration: state.iteration.iteration + offset } }); } });
        const extensionHistory = result.history, a = before.linearDiagnostics, b = result.linearDiagnostics;
        result.history = [...before.history, ...extensionHistory.slice(1).map(h => ({ ...h, iteration: h.iteration + offset }))];
        result.linearDiagnostics = { solves: a.solves + b.solves, refinements: a.refinements + b.refinements,
          pivotRecoveries: a.pivotRecoveries + b.pivotRecoveries,
          maxRelativeResidual: Math.max(a.maxRelativeResidual, b.maxRelativeResidual),
          iterations: [...a.iterations, ...b.iterations.map(h => ({ ...h, iteration: h.iteration + offset }))] };
        const { resume: unusedResume, ...evidence } = extensionPlan;
        result.startupExtension = { ...evidence,
          original: before.startupExtension?.original ?? evidence.original,
          chunks: [...(before.startupExtension?.chunks ?? []), { offset, budget: extensionPlan.additionalIterations,
            iterations: extensionHistory.length - 1, selection: extensionPlan.selection, progress: extensionPlan.progress }],
          iterations: (before.startupExtension?.iterations ?? 0) + extensionHistory.length - 1,
          converged: result.converged, reason: result.reason, families: { ...result.families } };
        extensionCount++;
      }
      // A coarse transition bubble can stall without reversing its active
      // interval. After the ordinary budget, inspect the actual BL limiters
      // and refine once from the best complete accepted state. The selection
      // heuristic changes work and resolution; all physical gates remain.
      const stallPlan = !direct && !recoveryPlan && maxIterations > 0 && startupAttempt === 1
        && activeNcrit === validated.options.ncrit && actualMach === (input.mach ?? 0)
        ? transitionRecoveryPlan(result, { bestState, tolerance }) : null;
      if (stallPlan?.reason === 'boundary-layer-transition-stall') {
        const terminal = result;
        coarseHistory = terminal.history;
        const offset = coarseHistory.length - 1;
        let observerFailed = false, observerError;
        const observe = callback => (...args) => {
          try { return callback?.(...args); }
          catch (error) { observerFailed = true; observerError = error; throw error; }
        };
        // This exact zero-update replay restores the chart, phase and full
        // geometry. It does not perform a factorization or initial SMOVE.
        const controls = bestState.checkpoint.continuation;
        const best = solveCoupledStreamtubeIses(undefined, { resume: bestState.checkpoint, maxIterations: 0, tolerance,
          iterationGeometry: controls.iterationGeometry, stepAcceptance: controls.stepAcceptance,
          stagnationLimiter: controls.stagnationLimiter });
        best.mesh.initialization.gridSmoothing = eulerInitialization.gridSmoothing;
        result = best;
        try {
          const recovered = recoverCoupledTransition(best, { plan: stallPlan,
            maxIterations: Math.min(maxIterations, 20), tolerance,
            onIteration: observe(h => iterationCallback({ ...h, iteration: h.iteration + offset, transitionRecovery: true })),
            onMesh: observe(onMesh), onStage: observe(detail => changeStage(detail.stage, detail)),
            onIterationCheckpoint: observe((cp, details) => onIterationCheckpoint?.(cp, { ...details, iterationOffset: offset })),
            onFlow: onFlow ? observe(flow => onFlow({ ...flow, iteration: { ...flow.iteration, iteration: flow.iteration.iteration + offset } })) : undefined,
            normalization: settings.normalization, startupAttempt });
          if (!recovered?.converged) {
            transitionRecoveryFailure = recovered?.reason ?? 'Local BL-stall refinement did not return a result.';
            result.recoveryAttempt = { converged: false, reason: transitionRecoveryFailure,
              families: recovered?.families, iterations: recovered?.history?.length - 1,
              plan: stallPlan, retained: 'best accepted unrefined state' };
          } else result = recovered;
        } catch (error) {
          if (observerFailed) throw observerError;
          transitionRecoveryFailure = error.message;
          result.recoveryAttempt = { converged: false, reason: transitionRecoveryFailure,
            plan: stallPlan, retained: 'best accepted unrefined state' };
        }
        if (result === best) {
          result.reason = `Local BL-stall recovery failed: ${transitionRecoveryFailure}. Retained the best accepted state.`;
          changeStage('coupled', { transitionRecoveryFailure, retainedBestIteration: bestState.iteration });
          // Restore both pressure and mesh after any abandoned finer preview.
          flowCheckpoint = best.checkpoint;
          const retainedSystem = createCoupledStreamtubeBody(best.solverInput, { ...best.checkpoint.restart.options,
            initialEuler: best.checkpoint.restart.initialEuler, initialBL: best.checkpoint.restart.initialBL });
          meshCallback({ system: retainedSystem.euler, flow: best.flow, nodes: best.flow.nodes,
            iteration: { ...best.history.at(-1), iteration: bestState.iteration } });
          onIterationCheckpoint?.(structuredClone(best.checkpoint), {
            stage: 'coupled', startupAttempt, retainedBestIteration: bestState.iteration, transitionRecoveryFailure });
        }
      }
      // Retain the original thin Ncrit startup. Only its complete stalled
      // endpoint earns one explicitly recorded change of shear coordinate.
      const shearPlan = direct ? null : planCoupledLogarithmicShearRecovery(result, {
        startupPlan: ncritRecoveryPlan, startupAttempt, maxIterations, tolerance });
      if (shearPlan) {
        const before = result, offset = shearPlan.originalIterations, controls = shearPlan.resume.continuation;
        const { resume: unusedResume, ...selection } = shearPlan;
        displayedMethod = { shearCoordinate: 'logarithmic',
          hkFloorLinearization: shearPlan.resume.restart.options.hkFloorLinearization ?? 'exact' };
        checkpointIterationOffset = offset;
        changeStage('coupled', { ...ncritLabels(), shearRecovery: selection });
        const continued = solveCoupledStreamtubeIses(undefined, { resume: shearPlan.resume,
          maxIterations: shearPlan.additionalIterations, tolerance,
          iterationGeometry: controls.iterationGeometry, stepAcceptance: controls.stepAcceptance,
          stagnationLimiter: controls.stagnationLimiter, ...checkpointObserver,
          onIteration: h => { if (h.iteration) iterationCallback({ ...h, iteration: h.iteration + offset,
            shearRecovery: { kind: shearPlan.kind, shearCoordinate: 'logarithmic' } }); },
          onMesh: state => { if (state.iteration.iteration) meshCallback({ ...state,
            iteration: { ...state.iteration, iteration: state.iteration.iteration + offset } }); } });
        const additionalIterations = continued.history.length - 1;
        const evidence = { ...selection, iterations: additionalIterations,
          totalIterations: offset + additionalIterations, converged: continued.converged,
          reason: continued.reason, families: { ...continued.families }, checkpointAvailable: !!continued.checkpoint };
        if (continued.checkpoint) {
          result = continued;
          result.history = [...before.history, ...continued.history.slice(1).map(h => ({ ...h, iteration: h.iteration + offset }))];
          const a = before.linearDiagnostics, b = continued.linearDiagnostics;
          result.linearDiagnostics = { solves: a.solves + b.solves, refinements: a.refinements + b.refinements,
            pivotRecoveries: a.pivotRecoveries + b.pivotRecoveries,
            maxRelativeResidual: Math.max(a.maxRelativeResidual, b.maxRelativeResidual),
            iterations: [...a.iterations, ...b.iterations.map(h => ({ ...h, iteration: h.iteration + offset }))] };
        } else {
          result = before;
          displayedMethod = { shearCoordinate: coupledResultShearCoordinate(before),
            hkFloorLinearization: before.checkpoint.restart.options.hkFloorLinearization ?? 'exact' };
          restorePreviousStartup();
        }
        result.shearRecovery = evidence;
        // Completion accepts the actual retained coordinate; it does not
        // reinterpret the original linear checkpoint as logarithmic.
        ncritRecoveryPlan = { ...ncritRecoveryPlan, shearCoordinate: coupledResultShearCoordinate(result),
          shearRecovery: evidence };
      }
      const reportGridSolve = () => {
        if (!gridSequence) return;
        // The nested levels log ends at preparation of the last seed. Keep
        // its history intact and report the actual retained solve separately.
        // Replace the record so prior continuation frames keep their Ncrit.
        gridSequence = { ...gridSequence, reachedTarget: result.converged,
          finalSolve: { converged: result.converged, reason: result.reason,
            actualNcrit: result.checkpoint?.restart?.options.ncrit ?? activeNcrit,
            families: { ...result.families }, cells: result.mesh.cells.length,
            iterations: result.history.length - 1 } };
        result.gridSequence = gridSequence;
        if (!result.automaticRefinement || result.automaticRefinement.kind === 'coarse-to-fine')
          result.automaticRefinement = gridSequence;
        if (!result.mesh.initialization.gridRefinement || result.mesh.initialization.gridRefinement.kind === 'coarse-to-fine')
          result.mesh.initialization.gridRefinement = gridSequence;
      };
      reportGridSolve();
      if (ncritRecoveryPlan && result.converged) {
        result.mesh.initialization.gridSmoothing = eulerInitialization.gridSmoothing;
        const completed = finishCoupledNcritStartup(result, { plan: ncritRecoveryPlan, maxIterations, tolerance,
          normalization: settings.normalization, startupAttempt,
          onStage: event => changeStage(event.stage, event), onIteration, onMesh, onFlow, onIterationCheckpoint });
        result = completed.result; ncritRecovery = { ...completed.diagnostics,
          actualThicknessFactor: prepared.initialization.thicknessFactor };
        activeNcrit = result.checkpoint.restart.options.ncrit;
      }
      if (ncritInitialization && result.converged && activeNcrit < validated.options.ncrit) {
        const continued = continueCoupledNcrit(result, { targetNcrit: validated.options.ncrit, phase: 'fine',
          maxIterations, tolerance, normalization: settings.normalization, startupAttempt,
          onStage: event => changeStage(event.stage, event), onIteration, onMesh, onFlow, onIterationCheckpoint });
        result = continued.result; ncritInitialization.fine = continued.diagnostics;
        activeNcrit = result.checkpoint.restart.options.ncrit;
        // Grid completion describes the retained root on this resolution;
        // the separate Ncrit status below may still mark the request unfinished.
        reportGridSolve();
      }
      if (ncritInitialization || ncritRecoveryPlan) {
        result.actualNcrit = activeNcrit; result.targetNcrit = validated.options.ncrit;
        result.ncritContinuation = { ...ncritInitialization, ...(ncritRecoveryPlan ? { startup: ncritRecovery ?? {
          ...ncritRecoveryPlan, attempted: true, reachedTarget: false, actualNcrit: activeNcrit,
          stateConverged: result.converged, reason: result.reason } } : {}),
          actualNcrit: activeNcrit, targetNcrit: validated.options.ncrit,
          reachedTarget: activeNcrit === validated.options.ncrit, stateConverged: result.converged };
        if (activeNcrit !== validated.options.ncrit) {
          result.stateConverged = result.converged; result.converged = false; result.status = 'unconverged';
          result.reason = `Requested Ncrit ${validated.options.ncrit} was not reached; retained Ncrit ${activeNcrit}.`;
        }
      }
      const recoveryCase = coupledResolutionRecoveryCase(input, result,
        { enabled: direct && resolutionRecovery, maxIterations });
      if (recoveryCase) {
        const recovered = recoverCoupledResolution(result, input, { sourceInput: recoveryCase,
          solve: solveAssembly, maxIterations, eulerMaxIterations, tolerance,
          onIteration, onMesh, onStage: event => changeStage(event.stage, event),
          onCheckpoint, onIterationCheckpoint, onFlow });
        if (recovered.result) {
          result = recovered.result;
          prepared = { ...prepared, initialization: recovered.initialization.boundaryLayer };
        } else {
          result.resolutionRecovery = recovered.diagnostics;
          restorePreviousStartup();
        }
      }
      if (direct && !result.converged && maxIterations > 0) {
        const before = result, offset = result.history.length - 1;
        checkpointIterationOffset = offset;
        const retried = retryCoupledEulerSeed({ result, input: solverInput, options: settings.options,
          initialEuler, initialization: prepared.initialization, maxIterations, tolerance,
          ...checkpointObserver,
          onIteration: h => iterationCallback({ ...h, iteration: h.iteration + offset, seedRecovery: 'euler-edge' }),
          onMesh: state => meshCallback({ ...state, iteration: { ...state.iteration,
            iteration: state.iteration.iteration + offset, seedRecovery: 'euler-edge' } }) });
        if (retried) {
          if (retried.result) {
            result = retried.result; prepared = retried.prepared;
            result.history = [...before.history, ...result.history.slice(1).map(h => ({ ...h, iteration: h.iteration + offset }))];
            const a = before.linearDiagnostics, b = result.linearDiagnostics;
            result.linearDiagnostics = { ...b, solves: a.solves + b.solves,
              refinements: a.refinements + b.refinements, pivotRecoveries: a.pivotRecoveries + b.pivotRecoveries,
              maxRelativeResidual: Math.max(a.maxRelativeResidual, b.maxRelativeResidual),
              iterations: [...a.iterations, ...b.iterations.map(h => ({ ...h, iteration: h.iteration + offset }))] };
          }
          result.seedRecovery = retried.diagnostics;
          if (!retried.result) { checkpointIterationOffset = 0; restorePreviousStartup(); }
        } else checkpointIterationOffset = 0;
      }
      const resolutionPlan = direct && maxIterations > 0 ? stalledBoundaryLayerResolutionPlan(result) : null;
      if (resolutionPlan) {
        const before = result;
        const discardedSeedHistory = before.seedRecovery?.accepted === false ? before.seedRecovery.history ?? [] : [];
        const seedOffset = before.history.length - 1;
        const precedingHistory = [...before.history, ...discardedSeedHistory.slice(1).map(h => ({
          ...h, iteration: h.iteration + seedOffset, seedRecovery: 'euler-edge', retained: false }))];
        const offset = precedingHistory.length - 1;
        before.mesh.initialization ??= {};
        before.mesh.initialization.gridSmoothing ??= eulerInitialization.gridSmoothing;
        let refined, observerFailure = false, failure;
        const observe = callback => (...args) => {
          try { return callback?.(...args); }
          catch (error) { observerFailure = true; throw error; }
        };
        try {
          refined = recoverCoupledTransition(result, { plan: resolutionPlan, maxIterations, tolerance,
            onIteration: observe(h => iterationCallback({ ...h, iteration: h.iteration + offset, transitionRecovery: true })),
            onMesh: observe(onMesh), onStage: observe(detail => changeStage(detail.stage, detail)),
            onIterationCheckpoint: observe((cp, details) => onIterationCheckpoint?.(cp, { ...details, iterationOffset: offset })),
            onFlow: onFlow ? observe(flow => onFlow({ ...flow, iteration: { ...flow.iteration, iteration: flow.iteration.iteration + offset } })) : undefined,
            normalization: settings.normalization, startupAttempt });
        } catch (error) { if (observerFailure) throw error; failure = error.message; }
        const diagnostics = { plan: resolutionPlan, accepted: refined?.converged === true,
          originalFamilies: before.families, originalIterations: seedOffset, precedingIterations: offset,
          iterations: refined ? refined.history.length - 1 : 0,
          reason: failure ?? refined?.reason, families: refined?.families };
        if (refined?.converged) {
          result = refined; result.seedRecovery = before.seedRecovery;
          result.history = [...precedingHistory, ...refined.history.slice(1).map(h => ({
            ...h, iteration: h.iteration + offset, transitionRecovery: true }))];
          const phases = [before.linearDiagnostics, before.seedRecovery?.accepted === false
            ? before.seedRecovery.linearDiagnostics : null, refined.linearDiagnostics].filter(Boolean);
          result.linearDiagnostics = { ...refined.linearDiagnostics,
            solves: phases.reduce((sum, d) => sum + d.solves, 0),
            refinements: phases.reduce((sum, d) => sum + d.refinements, 0),
            pivotRecoveries: phases.reduce((sum, d) => sum + d.pivotRecoveries, 0),
            maxRelativeResidual: Math.max(...phases.map(d => d.maxRelativeResidual)),
            iterations: phases.flatMap((d, i) => d.iterations.map(h => ({ ...h,
              iteration: h.iteration + (i === 0 ? 0 : i === phases.length - 1 ? offset : seedOffset) }))) };
        } else { result = before; restorePreviousStartup(); }
        result.transitionResolutionRecovery = diagnostics;
      }
      const machCase = coupledMachRecoveryCase(input, result, { enabled: direct && machRecovery, maxIterations });
      if (machCase) {
        const recovered = recoverCoupledMach(result, input, { sourceInput: machCase,
          solve: solveAssembly, maxIterations, eulerMaxIterations, tolerance,
          normalization: settings.normalization, onStage: event => changeStage(event.stage, event),
          onIteration, onMesh, onFlow, onIterationCheckpoint });
        if (recovered.result) result = recovered.result;
        else { result.machRecovery = recovered.diagnostics; restorePreviousStartup(); }
      }
      restart = coupledAssemblyRestart(result);
      attempts.push({ startupAttempt, thicknessFactor: prepared.initialization.thicknessFactor,
        converged: result.converged, reason: result.reason, families: result.families, quality: result.mesh.quality,
        iterations: (ncritRecovery?.iterations ?? result.history.length - 1) + (coarseHistory && (result.automaticRefinement || result.recoveryAttempt)
          ? coarseHistory.length - 1 + (result.recoveryAttempt?.iterations ?? 0) : 0),
        history: result.history, lastRejectedStep: result.lastRejectedStep,
        ...(result.startupExtension ? { startupExtension: result.startupExtension } : {}),
        ...(coarseHistory ? { coarseHistory, coarseIterations: coarseHistory.length - 1,
          automaticRefinement: result.automaticRefinement, transitionRecoveryFailure } : {}) });
      onCheckpoint?.({ stage, startupAttempt, ...ncritLabels(), families: result.families, restart });
      if (result.converged || maxIterations === 0 || result.automaticRefinement) break;
      // Keep a successful ordinary cold solve unchanged. Only its failed,
      // admitted complete state can select the second startup at lower Ncrit;
      // initialization exceptions and invalid grids cannot enter this route.
      // The selected equations, geometry and non-Ncrit conditions are reused.
      if (startupAttempt === 1 && maxStartupAttempts === 2 && activeNcrit > 4 && !coarseInitialization) {
        ncritRecoveryPlan = coupledNcritStartupRecoveryPlan({ ...input, ncrit: validated.options.ncrit,
          mach: input.mach ?? 0 }, result, { startupAttempt, maxStartupAttempts, maxIterations,
          thicknessFactor: prepared.initialization.thicknessFactor, coarseInitialization });
        if (ncritRecoveryPlan) {
          activeNcrit = ncritRecoveryPlan.sourceNcrit;
          initialThicknessFactor = ncritRecoveryPlan.thicknessFactor;
          continue;
        }
      }
      initialThicknessFactor = .25 * prepared.initialization.thicknessFactor;
    }
    return { ...result, model: 'research-streamtube-euler-bl', physicalAcceptance: false,
      alpha: input.alpha ?? 0, mach: actualMach, ...settings.normalization,
      materialTrips: settings.materialTrips, elementOrder: settings.elementOrder,
      bodies: result.solverInput.bodies, restart,
      initialization: { euler: eulerInitialization,
        boundaryLayer: prepared.initialization, attempts,
        ...(ncritRecoveryPlan ? { ncritStartupRecovery: result.ncritContinuation.startup } : {}),
        ...(coarseInitialization ? { coarseInitialization } : {}) },
      solverSettings: { maxIterations, eulerMaxIterations, maxStartupAttempts, tolerance, convergence, wakeGeometry: 'independent-banks', wakeOutlet: 'banks',
        shearCoordinate: coupledResultShearCoordinate(result),
        hkFloorLinearization: result.checkpoint?.restart.options.hkFloorLinearization
          ?? result.conditions?.hkFloorLinearization ?? 'exact',
        stepMethod: 'ises-density-newton', ...(direct ? { stepAcceptance: result.stepAcceptance } : {}), blUpdate: result.blUpdate ?? 'giles', automaticTransitionRefinement: !direct || !!result.transitionResolutionRecovery, automaticResolutionRecovery: resolutionRecovery, automaticMachRecovery: machRecovery, ...(direct ? { direct: true } : {}),
        // Before the first checkpoint, the driver still reports its actual
        // configured policy. A complete legacy checkpoint remains fixed.
        ...((result.checkpoint ? result.checkpoint.continuation?.projectionGeometry : result.projectionGeometry) === 'boundary-increment'
          ? { projectionGeometry: 'boundary-increment' } : {}),
        streamwiseMode: result.solverInput.streamwiseMode, edgeMatching: settings.options.edgeMatching,
        ...(input.eulerIsmom === undefined ? {} : { eulerIsmom: input.eulerIsmom,
          hybrid: structuredClone(result.solverInput.hybrid), upwind: structuredClone(result.solverInput.upwind),
          blThermodynamics: settings.options.blThermodynamics }),
        ...(coarseInitialization ? { coarseStartup: { maximumAttempts: 1, attempted: true,
          accepted: coarseInitialization.accepted } } : {}) },
      coefficientStatus: 'unavailable', forceStatus: 'Solid-wall forces and viscous drag are not validated for this model.' };
  } catch (error) {
    if (error && (typeof error === 'object' || typeof error === 'function') && Object.isExtensible(error)) error.stage ??= stage;
    throw error;
  }
}
