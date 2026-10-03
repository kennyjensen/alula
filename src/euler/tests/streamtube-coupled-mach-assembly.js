import { solveCoupledStreamtubeAlpha } from './streamtube-coupled-alpha.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// Public normalization/cache adapter for the explicit research hybrid route.
// The .3 cold-routing threshold, .2 baseline and hybrid controls below are
// workbench policies, not published MSES constants or a physical Mach limit.
import { targetMachGridPlan, solveCoupledTargetGrid } from './streamtube-coupled-target-grid.js';
import { coupledResultShearCoordinate, coupledCheckpointShearCoordinate } from '../streamtube-coupled-shear-policy.js';
import { coupledCheckpointHkPolicy, coupledAssemblyConditions, solveCoupledStreamtubeAssembly } from '../streamtube-coupled-assembly.js';
import { initializeCoupledStreamtubeFromFlow } from './streamtube-coupled-flow-restart.js';
import { certifyCoupledStreamtubeHybrid } from './streamtube-coupled-hybrid-certification.js';
import { solveCoupledStreamtubeAutomatic } from './streamtube-coupled-automatic.js';
import { solveCoupledStreamtubeIses } from '../streamtube-coupled-ises.js';
import { checkpointDataEqual } from '../streamtube-nested-checkpoint.js';
import { streamtubeEquationControls } from '../streamtube-equation-selection.js';
import { colderStreamtubeStartupMach } from './streamtube-cold-startup.js';
import { retargetPreparedStreamtubeAssembly } from '../streamtube-prepared-assembly.js';
import { normalizeStreamtubeEquationRegions, validateStreamtubeEquationRegionTopology } from '../streamtube-equation-selection.js';

const require = (condition, message) => { if (!condition) throw new Error(message); };
const copy = value => structuredClone(value);
const families = ['euler', 'boundaryLayer', 'edgeMatching'];
const withoutMach = value => { const result = { ...value }; delete result.mach; return result; };
const defaultUpwind = () => ({ mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } });
const normalization = result => Object.fromEntries(['referenceChord', 'referenceReynolds', 'solverLength', 'kernelReynolds']
  .map(key => [key, result[key]]));

// No numerical system, mesh, residual evaluation or initialization is created
// here. A supplied stale parent is an error; only omission selects cold start.
export function coupledMachPlan(caseData, parentResult, {
  tolerance = parentResult?.solverSettings?.tolerance ?? 1e-10,
  epsilonP = 1e-5, upwind = defaultUpwind(), entropyRegions,
} = {}) {
  require(caseData?.flowModel === 'streamtube-grid' && caseData.quadBoundaryLayers === true
    && Array.isArray(caseData.elements) && caseData.elements.length > 0,
  'Coupled Mach continuation requires the coupled quadrilateral flow model.');
  const targetMach = caseData.mach ?? .2;
  const selectedEquations = streamtubeEquationControls(caseData.eulerIsmom);
  require(Number.isFinite(targetMach) && targetMach > 0 && targetMach < 1
    && Number.isFinite(tolerance) && tolerance > 0, 'Invalid coupled Mach target or tolerance.');
  require(Number.isFinite(epsilonP) && epsilonP > 0 && upwind && typeof upwind === 'object'
    && Object.keys(upwind).every(k => ['mucon', 'mcrit', 'boundary'].includes(k))
    && Number.isFinite(upwind.mucon) && Number.isFinite(upwind.mcrit) && upwind.mcrit >= 0 && upwind.mcrit <= 1
    && upwind.boundary?.kind === 'unfiltered-first-two' && Object.keys(upwind.boundary).length === 1,
  'Invalid explicit research hybrid/upwind controls.');
  if (caseData.eulerIsmom !== undefined) require(epsilonP === selectedEquations.hybrid.epsilonP
    && checkpointDataEqual(upwind, selectedEquations.upwind),
  'Explicit ISMOM uses the public dissipation controls: epsilonP=1e-5, mucon=1 and Mcrit=.99.');
  coupledAssemblyConditions(caseData, caseData.elements.map((_, element) => ({ element })), caseData.referenceChord ?? 1);
  // The caller names the exact inherited equation region it intends to
  // retain. Omission still rejects regionalized parents; no cold route or
  // default logical region is silently converted into this research mode.
  let expectedRegions;
  if (entropyRegions !== undefined) {
    require(parentResult !== undefined && caseData.eulerIsmom === 3,
      'Explicit entropy regions require a warm parent and explicit ISMOM3.');
    expectedRegions = normalizeStreamtubeEquationRegions(entropyRegions);
  }
  if (parentResult === undefined) return { route: targetMach > .3 ? 'cold-baseline' : 'cold', targetMach,
    sourceMach: targetMach > .3 ? Math.min(.2, targetMach / 2) : targetMach,
    tolerance, epsilonP, upwind: copy(upwind) };
  require(parentResult?.model === 'research-streamtube-euler-bl' && parentResult.converged === true
    && parentResult.mesh?.quality?.valid === true, 'Mach continuation requires a converged coupled parent with a valid grid.');
  require(parentResult.sourceCase && checkpointDataEqual(withoutMach(parentResult.sourceCase), withoutMach(caseData)),
    'The coupled parent is stale: every case control except Mach must remain unchanged.');
  const checkpoint = parentResult.checkpoint, f = checkpoint?.restart, input = f?.input;
  require(checkpoint?.version === 1 && input && f.options && f.initialEuler?.x?.length
    && Array.isArray(f.initialEuler.nodes) && Array.isArray(f.initialEuler.undisplacedNodes)
    && f.initialBL?.length && checkpoint.continuation && checkpoint.families,
  'The coupled parent has no complete continuation checkpoint.');
  coupledCheckpointHkPolicy(caseData, f.options);
  require(checkpointDataEqual(checkpoint.families, parentResult.families), 'The coupled parent and checkpoint residuals differ.');
  require(Number.isFinite(parentResult.solverSettings?.tolerance) && parentResult.solverSettings.tolerance > 0
    && Object.keys(checkpoint.families).length === 3 && families.every(k => Number.isFinite(checkpoint.families[k])
      && checkpoint.families[k] >= 0 && checkpoint.families[k] <= Math.min(tolerance, parentResult.solverSettings.tolerance)),
  'The coupled parent must satisfy both source and requested convergence tolerances.');
  const sourceMach = parentResult.sourceCase.mach ?? .2;
  require(input.mach === sourceMach && (parentResult.mach === undefined || parentResult.mach === sourceMach)
    && (input.alpha ?? 0) === (caseData.alpha ?? 0), 'The coupled parent has inconsistent incidence or actual Mach metadata.');
  const settings = coupledAssemblyConditions(caseData, input.bodies, parentResult.solverLength, { checkpointOptions: f.options });
  require(checkpointDataEqual(settings.normalization, normalization(parentResult)),
    'The coupled parent changed the public physical normalization.');
  for (const key of ['reynolds', 'ncrit', 'edgeMatching', 'tripFractions'])
    require(checkpointDataEqual(settings.options[key], f.options[key]), `The coupled case and checkpoint disagree on ${key}.`);
  require((settings.options.transitionMode ?? 'fixed-trip') === (f.options.transitionMode ?? 'fixed-trip'),
    'The coupled case and checkpoint disagree on transition mode.');
  if (f.options.transitionMode === 'automatic') require(Array.isArray(f.options.transitionState)
    && f.options.transitionState.length === 2 * input.bodies.length && f.options.transitionState.every(Number.isInteger),
  'The coupled checkpoint has no complete automatic transition map.');
  if (expectedRegions !== undefined) {
    require(input.hybrid?.ismom === 3 && Array.isArray(input.outerLower) && Array.isArray(input.weights)
      && input.weights.every(row => Array.isArray(row)), 'Explicit entropy regions require an ISMOM3 checkpoint with complete grid topology.');
    validateStreamtubeEquationRegionTopology(expectedRegions, { nx: input.outerLower.length - 1,
      tubes: input.weights.map(row => row.length), bodies: input.bodies });
  }
  const baseline = caseData.eulerIsmom === undefined && input.streamwiseMode === 'isentropic' && input.upwind === undefined && input.hybrid === undefined
    && f.options.blThermodynamics === undefined;
  const hybrid = input.streamwiseMode === 'hybrid' && f.options.blThermodynamics === 'historical-common-isentrope'
    && checkpointDataEqual(input.hybrid, { epsilonP, ...(caseData.eulerIsmom === undefined ? {} : { ismom: caseData.eulerIsmom }),
      ...(expectedRegions === undefined ? {} : { entropyRegions: expectedRegions }) })
    && checkpointDataEqual(input.upwind, upwind);
  require((input.flowModel === undefined || input.flowModel === 'compressible') && (baseline || hybrid),
  'The coupled parent must use certifiable isentropic equations or the unchanged explicit historical hybrid controls.');
  for (const key of ['streamwiseMode', 'wakeGeometry', 'wakeOutlet']) if (parentResult.solverSettings[key] !== undefined)
    require(parentResult.solverSettings[key] === (input[key] ?? (key === 'streamwiseMode' ? undefined : 'centerline')),
      `The coupled parent has stale ${key} solver metadata.`);
  return { route: baseline ? sourceMach === targetMach ? 'warm-isentropic' : 'warm-certified' : 'warm-hybrid', sourceMach, targetMach,
    tolerance, epsilonP, upwind: copy(upwind), settings,
    ...(expectedRegions === undefined ? {} : { entropyRegions: copy(expectedRegions) }) };
}

export function solveCoupledStreamtubeMach(caseData, { parentResult, preparedEuler,
  tolerance = parentResult?.solverSettings?.tolerance ?? 1e-10,
  epsilonP = 1e-5, upwind = defaultUpwind(), entropyRegions, stageMaxIterations = 20, maxMachStep = .05, maxStages = 48, maxSubdivisions = 5,
  maxFirstOrderRecoveries = 2, maxWakeRecoveries = 4, maxBacktracks = 12, blPredictor, maxIterations = caseData?.maxIterations ?? (caseData?.transitionMode === 'automatic' ? 40 : 20),
  eulerMaxIterations = caseData?.eulerMaxIterations ?? 20, maxStartupAttempts = 2,
  onStage, onIteration, onMesh, onCheckpoint, onPrepared, onFlow, stopAtShock = false,
  alphaContinuation = true, maxAlphaStep = .1, targetInviscidStartup = true, targetGridSequencing = true, dissipationEnhancement = true, iterationRecovery = true,
  maxProgressExtraIterations = caseData?.maxIterations === undefined ? 8 : 0,
} = {}) {
  const targetGridPlan = targetMachGridPlan(caseData, { parentResult, preparedEuler, enabled: targetGridSequencing });
  if (typeof alphaContinuation !== 'boolean' || !Number.isFinite(maxAlphaStep) || maxAlphaStep <= 0)
    throw new Error('Invalid alpha continuation policy.');
  if (typeof stopAtShock !== 'boolean') throw new Error('Invalid shock-refinement handoff policy.');
  if (typeof targetInviscidStartup !== 'boolean') throw new Error('Invalid target inviscid startup policy.');
  if (typeof dissipationEnhancement !== 'boolean') throw new Error('Invalid dissipation enhancement policy.');
  if (typeof iterationRecovery !== 'boolean' || !Number.isInteger(maxProgressExtraIterations)
    || maxProgressExtraIterations < 0 || maxProgressExtraIterations > 20) throw new Error('Invalid progress recovery policy.');
  const resumeOperatingPoint = (seed, observers) => {
    const result = solveCoupledStreamtubeAlpha(caseData.alpha ?? 0, {
      initialCheckpoint: seed.checkpoint, targetMach: caseData.mach,
      maxMachStep, maxAlphaStep, tolerance, stageMaxIterations, maxStages, maxSubdivisions,
      normalization: normalization(seed), stopAtShock, preferAlphaApproach: parentResult === undefined, ...observers,
    });
    return { ...seed, ...result, sourceCase: { ...copy(seed.sourceCase),
      mach: result.mach, alpha: result.alpha }, requestedCase: copy(caseData),
      restart: result.checkpoint.restart,
      machContinuation: { ...seed.machContinuation, actualMach: result.mach, targetMach: caseData.mach,
        reachedTarget: result.mach === caseData.mach },
      initialization: seed.initialization,
      solverSettings: seed.solverSettings,
      ...(result.converged ? {} : { status: 'research-coupled-target-not-reached' }) };
  };
  if (targetGridPlan) return solveCoupledTargetGrid(caseData, { plan: targetGridPlan, parentResult, maxIterations, tolerance,
    dissipationEnhancement, iterationRecovery, onStage, onIteration, onMesh, onCheckpoint, onPrepared, onFlow,
    resumeOperatingPoint,
    solveCoarse: (coarseCase, observers) => solveCoupledStreamtubeMach(coarseCase, {
      tolerance, epsilonP, upwind, entropyRegions, stageMaxIterations, maxMachStep, maxStages, maxSubdivisions,
      maxFirstOrderRecoveries, maxWakeRecoveries, maxBacktracks, blPredictor, maxIterations, eulerMaxIterations, maxStartupAttempts,
      alphaContinuation, maxAlphaStep, targetInviscidStartup, targetGridSequencing: false, dissipationEnhancement, iterationRecovery, maxProgressExtraIterations, ...observers }) });
  let alphaSource = parentResult?.sourceCase?.alpha ?? Math.sign(caseData.alpha ?? 0) * Math.max(0, Math.abs(caseData.alpha ?? 0) - 1);
  if (alphaContinuation && parentResult !== undefined && alphaSource !== caseData.alpha
    && parentResult.sourceCase?.mach !== caseData.mach) {
    coupledMachPlan({ ...caseData, alpha: alphaSource }, parentResult, { tolerance, epsilonP, upwind, entropyRegions });
    return resumeOperatingPoint(parentResult, { onStage, onIteration, onCheckpoint, onMesh, onFlow });
  }
  const warmAlpha = parentResult !== undefined && parentResult.sourceCase?.mach === caseData?.mach
    && alphaSource !== caseData?.alpha;
  if (alphaContinuation && preparedEuler === undefined && Number.isFinite(caseData?.alpha)
    && (warmAlpha || parentResult === undefined && caseData.mach >= .6 && caseData.alpha !== 0
      && maxIterations > 0 && eulerMaxIterations > 0)) {
    const targetAlpha = caseData.alpha;
    const seedAtIncidence = () => solveCoupledStreamtubeMach({ ...copy(caseData), alpha: alphaSource }, {
      parentResult, alphaContinuation: false, stopAtShock, targetInviscidStartup, targetGridSequencing: false,
      tolerance, epsilonP, upwind, entropyRegions, stageMaxIterations, maxMachStep, maxStages, maxSubdivisions,
      maxFirstOrderRecoveries, maxWakeRecoveries, maxBacktracks, blPredictor, maxIterations, eulerMaxIterations,
      maxStartupAttempts, dissipationEnhancement, iterationRecovery, maxProgressExtraIterations,
      onStage: e => onStage?.({ ...e, alpha: alphaSource, actualAlpha: alphaSource, targetAlpha, alphaSeed: true }),
      onIteration: e => onIteration?.({ ...e, alpha: alphaSource, actualAlpha: alphaSource, targetAlpha, alphaSeed: true }),
      onMesh, onFlow,
      // Intermediate incidence must never be cached under requested alpha.
      onCheckpoint: (cp, e) => onCheckpoint?.(cp, { ...e, alpha: alphaSource, actualAlpha: alphaSource, targetAlpha,
        kind: e.kind === 'accepted' ? 'alpha-seed' : e.kind, reachedTarget: false }),
    });
    let seed = seedAtIncidence();
    if (!warmAlpha && !seed.converged && !seed.stateConverged && !seed.refinementRequested && alphaSource !== 0) {
      alphaSource *= .5;
      onStage?.({ stage: 'coupled-cold-recovery', actualAlpha: alphaSource, targetAlpha,
        mach: caseData.mach, targetMach: caseData.mach,
        reason: 'Retry a lower incidence seed for alpha continuation.' });
      seed = seedAtIncidence();
    }
    if (seed.refinementRequested) return { ...seed, requestedCase: copy(caseData),
      alphaContinuation: { sourceAlpha: alphaSource, actualAlpha: alphaSource, targetAlpha, reachedTarget: false } };
    if ((seed.converged || seed.stateConverged) && seed.checkpoint?.restart.input.streamwiseMode === 'hybrid') {
      if (!warmAlpha) return resumeOperatingPoint(seed, { onStage, onIteration, onCheckpoint, onMesh, onFlow });
      const result = solveCoupledStreamtubeAlpha(targetAlpha, { initialCheckpoint: seed.checkpoint,
        tolerance, maxAlphaStep, stageMaxIterations, onStage, onIteration, onCheckpoint, onMesh, onFlow,
        normalization: normalization(seed) });
      if (!result.converged && result.stateConverged)
        return resumeOperatingPoint({ ...seed, ...result,
          sourceCase: { ...copy(caseData), alpha: result.alpha } },
        { onStage, onIteration, onCheckpoint, onMesh, onFlow });
      return { ...seed, ...result, mach: caseData.mach, actualMach: caseData.mach,
        sourceCase: { ...copy(caseData), alpha: result.alpha }, requestedCase: copy(caseData),
        initialization: seed.initialization, solverSettings: seed.solverSettings,
        ...(result.converged ? {} : { status: 'research-coupled-target-not-reached' }) };
    }
    onStage?.({ stage: 'coupled-cold-recovery', alpha: alphaSource, actualAlpha: alphaSource,
      targetAlpha, mach: seed.mach ?? caseData.mach, targetMach: caseData.mach,
      reason: 'Lower-incidence seed did not converge; preserving the alpha-startup failure instead of switching to final-alpha Mach continuation.' });
    return { ...seed, converged: false, requestedCase: copy(caseData),
      alphaContinuation: { sourceAlpha: alphaSource, actualAlpha: alphaSource, targetAlpha, reachedTarget: false },
      reason: `Lower-incidence startup at ${alphaSource}° failed: ${seed.reason}` };
  }
  const requestedCase = copy(caseData), plan = coupledMachPlan(requestedCase, parentResult, { tolerance, epsilonP, upwind, entropyRegions });
  // A cold recovery inside a warm request inherits the accepted source's
  // Jacobian policy; it is not an unrelated fresh-default selection.
  const inheritedHkPolicy = parentResult === undefined ? undefined
    : coupledCheckpointHkPolicy(requestedCase, parentResult.checkpoint.restart.options);
  const coldSolverCase = value => inheritedHkPolicy === undefined ? value
    : { ...value, coupledNativeHk: inheritedHkPolicy === 'native' };
  require(preparedEuler === undefined || parentResult === undefined, 'A warm coupled restart cannot also supply a cold prepared mesh.');
  require([stageMaxIterations, maxIterations, eulerMaxIterations].every(n => Number.isInteger(n) && n >= 0)
    && Number.isFinite(maxMachStep) && maxMachStep > 0 && maxMachStep <= 1
    && Number.isInteger(maxStages) && maxStages >= 1 && maxStages <= 256
    && Number.isInteger(maxSubdivisions) && maxSubdivisions >= 0 && maxSubdivisions <= 20
    && Number.isInteger(maxFirstOrderRecoveries) && maxFirstOrderRecoveries >= 0 && maxFirstOrderRecoveries <= 4
    && Number.isInteger(maxWakeRecoveries) && maxWakeRecoveries >= 0 && maxWakeRecoveries <= 4
    && Number.isInteger(maxBacktracks) && maxBacktracks >= 0 && maxBacktracks <= 20
    && [1, 2].includes(maxStartupAttempts)
    && (blPredictor === undefined || ['preserve', 'xfoil-mrchdu'].includes(blPredictor))
    && [onStage, onIteration, onMesh, onCheckpoint, onPrepared, onFlow].every(fn => fn === undefined || typeof fn === 'function'),
  'Invalid public coupled Mach-continuation controls.');
  let observerFailed = false, observerError, latestMesh, targetStartup;
  const observe = (fn, ...args) => {
    if (!fn) return;
    try { fn(...args.map(copy)); }
    catch (error) { observerFailed = true; observerError = error; throw error; }
  };
  const info = (mach, details = {}) => ({ ...details, mach, actualMach: mach, targetMach: plan.targetMach });
  // This new opt-in reports the checkpoint's actual linear controls, not
  // potentially stale labels inherited from an earlier public result.
  const regionalLinearSettings = checkpoint => plan.entropyRegions === undefined ? {} : {
    linearOrdering: checkpoint.continuation.linearOrdering ?? 'auto',
    pivotTolerance: checkpoint.continuation.pivotTolerance,
  };
  const projectionSettings = (checkpoint, previous, reportedPolicy) => {
    // Only a result returned before checkpoint creation needs this fallback.
    // An existing checkpoint without the optional control means fixed.
    const policy = checkpoint ? checkpoint.continuation?.projectionGeometry ?? 'fixed' : reportedPolicy ?? 'fixed';
    return policy === 'boundary-increment' || previous?.projectionGeometry !== undefined
      ? { projectionGeometry: policy } : {};
  };
  const attach = (result, actualCase, source, details = {}) => {
    const mach = result.conditions?.mach ?? result.mach ?? actualCase.mach;
    require(mach === actualCase.mach, 'Returned coupled flow and actual case Mach disagree.');
    const actualNcrit = result.checkpoint?.restart?.options?.ncrit ?? result.conditions?.ncrit ?? actualCase.ncrit ?? 9;
    const targetNcrit = requestedCase.ncrit ?? 9, ncritContinuation = result.ncritContinuation ?? source.ncritContinuation;
    const ncritTracked = ncritContinuation !== undefined || result.actualNcrit !== undefined || actualNcrit !== targetNcrit;
    const ncritReached = actualNcrit === targetNcrit && ncritContinuation?.reachedTarget !== false;
    if (ncritTracked) actualCase = { ...actualCase, ncrit: actualNcrit };
    const reachedTarget = mach === plan.targetMach && ncritReached && result.converged === true && result.mesh?.quality?.valid === true;
    const settings = result.solverInput
      ? coupledAssemblyConditions(actualCase, result.solverInput.bodies, source.solverLength,
        { checkpointOptions: result.checkpoint?.restart.options ?? source.checkpoint?.restart.options }) : null;
    return { ...result, model: 'research-streamtube-euler-bl', physicalAcceptance: false, fullSolverComplete: false,
      converged: reachedTarget, stateConverged: result.stateConverged ?? (result.converged === true && result.mesh?.quality?.valid === true),
      ...(reachedTarget ? {} : { status: 'research-coupled-target-not-reached' }),
      mach, actualMach: mach, targetMach: plan.targetMach, alpha: actualCase.alpha ?? 0,
      ...(ncritTracked ? { actualNcrit, targetNcrit,
        ...(ncritContinuation ? { ncritContinuation: copy(ncritContinuation) } : {}) } : {}),
      ...(settings ? { ...settings.normalization, materialTrips: settings.materialTrips, elementOrder: settings.elementOrder } : {}),
      ...(result.solverInput ? { bodies: result.solverInput.bodies } : {}),
      ...(result.checkpoint ? { restart: copy(result.checkpoint.restart) } : {}),
      sourceCase: copy(actualCase), ...(reachedTarget ? {} : { requestedCase: copy(requestedCase) }),
      initialization: copy(source.initialization), ...(source.refinement ? { refinement: copy(source.refinement) } : {}),
      solverSettings: { ...copy(source.solverSettings), tolerance, ...details.solverSettings,
        shearCoordinate: coupledResultShearCoordinate(result.checkpoint ? result : source),
        hkFloorLinearization: coupledCheckpointHkPolicy(actualCase, result.checkpoint?.restart.options
          ?? source.checkpoint?.restart.options ?? result.coupledOptions ?? result.conditions
          ?? source.coupledOptions ?? source.conditions),
        ...regionalLinearSettings(result.checkpoint ?? source.checkpoint),
        ...projectionSettings(result.checkpoint ?? source.checkpoint, { ...source.solverSettings, ...details.solverSettings },
          result.projectionGeometry ?? source.projectionGeometry) },
      machContinuation: { route: plan.route, sourceMach: plan.sourceMach, actualMach: mach, targetMach: plan.targetMach,
        reachedTarget: !ncritReached ? mach === plan.targetMach && (result.stateConverged ?? result.converged) === true
          && result.mesh?.quality?.valid === true : reachedTarget,
        coldBaselineUsed: plan.route.endsWith('cold-baseline'),
        ...(targetStartup ? { targetInviscidStartup: copy(targetStartup) } : {}), ...details.diagnostics },
      coefficientStatus: 'unavailable', forceStatus: 'Solid-wall forces and total drag require the explicit research force adapter.' };
  };
  const coldSourceAttempt = (mach, route, preparedEuler, onEulerPrepared) => {
    const coldCase = { ...copy(requestedCase), mach };
    const targetAttempt = route === 'cold-target-inviscid';
    const startupInfo = event => info(mach, { ...event, route,
      ...(targetAttempt ? { startupStrategy: 'target-inviscid-then-viscous' } : {}) });
    try {
      const source = solveCoupledStreamtubeAssembly(coldSolverCase(coldCase), { maxIterations, eulerMaxIterations, tolerance, maxStartupAttempts,
        preparedEuler, onEulerPrepared,
        ...(targetAttempt ? { coarseStartup: false, maxStartupAttempts: 1 } : {}),
        onStage: event => observe(onStage, startupInfo(event)),
        onIteration: event => observe(onIteration, startupInfo(event)),
        ...(onFlow ? { onFlow: frame => observe(onFlow, { ...frame, ...info(mach) }) } : {}),
        onMesh: (mesh, phase, stage) => {
          latestMesh = copy(mesh);
          observe(onMesh, { ...mesh, ...info(mach) }, phase, info(mach, { stage }));
        } });
      if (observerFailed) throw observerError;
      const actualNcrit = source.checkpoint?.restart?.options?.ncrit ?? source.conditions?.ncrit ?? source.actualNcrit;
      const ncritTracked = source.ncritContinuation !== undefined || source.actualNcrit !== undefined
        || Number.isFinite(actualNcrit) && actualNcrit !== (coldCase.ncrit ?? 9);
      source.sourceCase = ncritTracked ? { ...coldCase, ncrit: actualNcrit } : coldCase;
      if (ncritTracked && (actualNcrit !== (coldCase.ncrit ?? 9) || source.ncritContinuation?.reachedTarget === false)) {
        source.stateConverged ??= source.converged === true && source.mesh?.quality?.valid === true;
        source.converged = false;
        source.requestedCase = copy(coldCase);
      }
      return { source };
    } catch (error) {
      if (observerFailed) throw observerError;
      return { error, failure: { model: 'research-streamtube-euler-bl', status: 'research-coupled-initialization-failed', converged: false,
        stateConverged: false, availableFlow: false, physicalAcceptance: false, fullSolverComplete: false,
        ...info(mach), sourceCase: coldCase, requestedCase: copy(requestedCase), mesh: latestMesh,
        ...(error?.actualNcrit === undefined ? {} : { actualNcrit: error.actualNcrit, targetNcrit: requestedCase.ncrit ?? 9 }),
        reason: error?.message ?? String(error), failure: { stage: error?.stage, code: error?.code,
          ...(error?.actualNcrit === undefined ? {} : { actualNcrit: error.actualNcrit, targetNcrit: requestedCase.ncrit ?? 9 }),
          ...(error?.diagnostics === undefined ? {} : { diagnostics: copy(error.diagnostics) }),
          ...(error?.cause === undefined ? {} : { cause: { message: error.cause?.message ?? String(error.cause),
            code: error.cause?.code, stage: error.cause?.stage, diagnostics: copy(error.cause?.diagnostics) } }) },
        initialization: copy(error?.initialization), solverSettings: { maxIterations, eulerMaxIterations, maxStartupAttempts, tolerance },
        machContinuation: { route, sourceMach: mach, actualMach: mach,
          targetMach: plan.targetMach, reachedTarget: false, coldBaselineUsed: route.endsWith('cold-baseline') } } };
    }
  };
  const coldSource = (initialMach, route, initialPrepared) => {
    let mach = initialMach, prepared = initialPrepared, originalPrepared = initialPrepared;
    const attempts = [];
    for (;;) {
      const result = coldSourceAttempt(mach, route, prepared, packet => {
        prepared = packet; originalPrepared ??= packet;
      });
      if (attempts.length) {
        const output = result.failure ?? result.source;
        output.initialization = { ...output.initialization, coldMachStartup: {
          initialMach, actualMach: mach, targetMach: plan.targetMach, meshReused: true, attempts: copy(attempts) } };
      }
      // A lower baseline is useful only when the public controller will
      // subsequently continue back to the requested operating condition.
      const nextMach = route.endsWith('cold-baseline') && prepared && result.error
        ? colderStreamtubeStartupMach(result.error, mach) : null;
      if (nextMach === null) return { ...result, mach };
      attempts.push({ mach, code: result.error.code, reason: result.error.message,
        diagnostics: copy(result.error.diagnostics) });
      const nextCase = coldSolverCase({ ...copy(requestedCase), mach: nextMach });
      try { prepared = retargetPreparedStreamtubeAssembly(originalPrepared, nextCase); }
      catch (error) {
        // A failed geometry handoff is terminal. Retain the original gas
        // rejection and actual Mach; no replacement state has been solved.
        result.failure.reason += ` Automatic colder initialization could not retain the prepared mesh: ${error.message}`;
        result.failure.initialization = { ...result.failure.initialization, coldMachStartup: {
          initialMach, actualMach: mach, targetMach: plan.targetMach, meshReused: false, attempts: copy(attempts),
          meshReuseFailure: { message: error.message, ...(error.code ? { code: error.code } : {}) } } };
        return { ...result, mach };
      }
      observe(onStage, info(nextMach, { stage: 'euler-cold-startup', previousMach: mach,
        startupAttempt: attempts.length + 1, reason: result.error.message, meshReused: true, route }));
      mach = nextMach;
    }
  };
  let source = parentResult === undefined ? undefined : copy(parentResult);
  if (!source) {
    // MSES §5.3: first move the inviscid flow at the requested operating
    // condition, then introduce the requested viscosity in one step.
    // The .6 trigger is workbench policy, not a published sonic threshold.
    // Fine requests perform this on the coarse target-Mach level above.
    let fallbackPrepared = preparedEuler;
    if (targetInviscidStartup && plan.targetMach >= .6 && maxIterations > 0 && eulerMaxIterations > 0) {
      let targetPrepared;
      const target = coldSourceAttempt(plan.targetMach, 'cold-target-inviscid', preparedEuler,
        packet => { targetPrepared = packet; });
      targetStartup = { attempted: true, mach: plan.targetMach, alpha: requestedCase.alpha,
        reynolds: requestedCase.reynolds, converged: target.source?.converged === true,
        reason: target.error?.message ?? target.source?.reason,
        initialization: copy(target.source?.initialization ?? target.failure?.initialization) };
      if (target.source?.converged && target.source.mesh?.quality?.valid) {
        plan.route = 'cold-target-inviscid'; plan.sourceMach = plan.targetMach;
        const result = attach(target.source, target.source.sourceCase, target.source);
        if (result.stateConverged && result.checkpoint) observe(onCheckpoint, result.checkpoint,
          info(plan.targetMach, { kind: 'accepted', stage: 'coupled', reachedTarget: result.converged }));
        return result;
      }
      if (!fallbackPrepared && targetPrepared) {
        try { fallbackPrepared = retargetPreparedStreamtubeAssembly(targetPrepared,
          coldSolverCase({ ...requestedCase, mach: plan.sourceMach })); }
        catch (error) { targetStartup.meshReuseFailure = error.message; }
      }
      targetStartup.meshReusedForFallback = fallbackPrepared !== undefined;
      observe(onStage, info(plan.sourceMach, { stage: 'coupled-cold-recovery',
        reason: targetStartup.reason, startupStrategy: 'low-mach-fallback' }));
    }
    const cold = coldSource(plan.sourceMach, plan.route, fallbackPrepared);
    if (cold.failure) {
      if (targetStartup) cold.failure.machContinuation.targetInviscidStartup = copy(targetStartup);
      return cold.failure;
    }
    plan.sourceMach = cold.mach;
    source = cold.source;
    if (plan.route === 'cold' || !source.converged || !source.mesh.quality.valid) {
      const result = attach(source, source.sourceCase, source);
      if (result.stateConverged && result.checkpoint) observe(onCheckpoint, result.checkpoint,
        info(plan.sourceMach, { kind: 'accepted', stage: 'coupled', reachedTarget: result.converged }));
      return result;
    }
    // Reuse the same guards for the public cold result. No extra cold attempt
    // is introduced when certification or continuation rejects a state.
    coupledMachPlan(requestedCase, source, { tolerance, epsilonP: plan.epsilonP, upwind: plan.upwind,
      ...(plan.entropyRegions === undefined ? {} : { entropyRegions: plan.entropyRegions }) });
  }
  const replayIsentropic = source => {
    // Repeating an accepted operating point must retain its equations and
    // original force adapter. The unchanged ISES resume validates complete
    // state/history and skips initial SMOVE; it performs no Newton update.
    const checkpoint = copy(source.checkpoint);
    const { iterationGeometry, stepAcceptance, stagnationLimiter } = checkpoint.continuation;
    const result = solveCoupledStreamtubeIses(undefined, { resume: checkpoint, iterationGeometry,
      stepAcceptance, stagnationLimiter, maxIterations: 0, tolerance, maxBacktracks });
    require(result.converged && result.mesh.quality.valid && result.conditions.mach === checkpoint.restart.input.mach
      && checkpointDataEqual(result.families, checkpoint.families)
      && result.residual.every(r => Number.isFinite(r) && Math.abs(r) <= tolerance),
    'The same-Mach isentropic source must replay exactly as an admissible converged convex state.');
    const savedPoints = checkpoint.restart.initialEuler.nodes, actualPoints = result.flow.nodes;
    require(savedPoints.length === actualPoints.length && savedPoints.every((group, g) =>
      group.length === actualPoints[g].length && group.every((row, i) => row.length === actualPoints[g][i].length
        && row.every((point, j) => point.x === actualPoints[g][i][j].x && point.y === actualPoints[g][i][j].y))),
    'The same-Mach isentropic source does not replay its actual physical nodes exactly.');
    const settings = coupledAssemblyConditions(source.sourceCase, result.solverInput.bodies, result.conditions.referenceChord,
      { checkpointOptions: checkpoint.restart.options });
    require(checkpointDataEqual(settings.normalization, normalization(source)),
      'Same-Mach isentropic replay changed the public physical normalization.');
    return result;
  };
  if (plan.route === 'warm-isentropic') {
    const result = replayIsentropic(source);
    observe(onCheckpoint, result.checkpoint, info(plan.sourceMach, { kind: 'accepted', stage: 'coupled-source', reachedTarget: true }));
    return attach(result, { ...source.sourceCase, mach: plan.sourceMach }, source, { diagnostics: { sameMachReplay: true,
      modelUnchanged: true, certificationPerformed: false, initialRedistributionSkippedOnWarmRestart: true } });
  }
  let checkpoint = copy(source.checkpoint), certification, coldRecovery;
  // At most two certifications: supplied source, then one fresh lower-Mach
  // baseline. The baseline's own existing startup policy is unchanged.
  for (let certificationAttempt = 0; checkpoint.restart.input.streamwiseMode === 'isentropic' && certificationAttempt < 2; certificationAttempt++) {
    observe(onStage, info(plan.sourceMach, { stage: 'hybrid-certification', startupAttempt: 1 }));
    try {
      const certified = certifyCoupledStreamtubeHybrid(checkpoint, { tolerance, epsilonP: plan.epsilonP, upwind: plan.upwind });
      checkpoint = certified.checkpoint; certification = certified.diagnostics;
      break;
    } catch (error) {
      if (plan.route.endsWith('cold-baseline')) return attach({ ...source,
        reason: `Cold baseline could not be certified for hybrid continuation: ${error.message}` },
      source.sourceCase, source, { diagnostics: { coldRecovery, failureStage: 'hybrid-certification', reason: error.message } });
      if (error?.code !== 'coupled-hybrid-certification-speed-bias' || coldRecovery) throw error;
      // A typed model incompatibility is not permission to hide corrupt
      // history, residuals, phase or physical geometry. Validate the complete
      // ORIGINAL isentropic source before entering the single recovery route.
      const retained = replayIsentropic(source), retainedSource = source;
      const mach = Math.min(.2, plan.targetMach / 2);
      coldRecovery = { reason: error.message, code: error.code, bias: copy(error.diagnostics),
        originalSourceMach: source.sourceCase.mach, baselineMach: mach,
        originalSourceReplayed: true, originalSourceConverged: true, originalSourceUnchanged: true, coldAttempts: 1 };
      observe(onStage, info(mach, { stage: 'coupled-cold-recovery', startupAttempt: 1, recovery: coldRecovery }));
      const cold = coldSource(mach, 'incompatible-source-cold-baseline');
      if (cold.failure) return attach({ ...retained, reason: `Lower-Mach recovery failed: ${cold.failure.reason}` },
        { ...retainedSource.sourceCase, mach: plan.sourceMach }, retainedSource, { diagnostics: {
          route: 'incompatible-source-cold-baseline', coldBaselineUsed: true, coldRecovery, failureStage: 'cold-initialization' } });
      source = cold.source;
      plan.route = 'incompatible-source-cold-baseline'; plan.sourceMach = cold.mach;
      if (cold.mach !== mach) Object.assign(coldRecovery, { initialBaselineMach: mach, baselineMach: cold.mach });
      if (!source.converged || !source.mesh.quality.valid) return attach(source, source.sourceCase, source,
        { diagnostics: { coldRecovery, failureStage: 'cold-coupled' } });
      coupledMachPlan(requestedCase, source, { tolerance, epsilonP: plan.epsilonP, upwind: plan.upwind,
        ...(plan.entropyRegions === undefined ? {} : { entropyRegions: plan.entropyRegions }) });
      checkpoint = copy(source.checkpoint);
    }
  }
  // Native mixed-profile preparation applies only to the verified automatic
  // terminal-trip XFOIL path. Earlier physical trips keep their existing
  // material mapping; the geometry-free MRCHDU adapter cannot represent them.
  const predictor = blPredictor ?? (checkpoint.continuation.blUpdate === 'xfoil'
    && checkpoint.restart.options.transitionMode === 'automatic'
    && checkpoint.restart.options.tripFractions?.every(pair => pair.every(v => v === 1))
    ? 'xfoil-mrchdu' : 'preserve');
  const solverSettings = { ...copy(source.solverSettings), tolerance, streamwiseMode: 'hybrid',
    shearCoordinate: coupledCheckpointShearCoordinate(checkpoint),
    hkFloorLinearization: coupledCheckpointHkPolicy(source.sourceCase, checkpoint.restart.options),
    hybrid: copy(checkpoint.restart.input.hybrid), upwind: copy(plan.upwind), blThermodynamics: 'historical-common-isentrope',
    edgeMatching: 'section-velocity', stepMethod: 'ises-density-newton',
    maxIterations: stageMaxIterations, stageMaxIterations, maxMachStep, maxStages, maxSubdivisions, maxFirstOrderRecoveries, maxWakeRecoveries, maxBacktracks,
    ...(predictor === 'preserve' ? {} : { blPredictor: predictor }),
    ...regionalLinearSettings(checkpoint), ...projectionSettings(checkpoint, source.solverSettings) };
  // Reconstruct only the supplied complete source. This verifies the physical
  // length normalization and provides worker BL topology, never a cold seed.
  {
    const f = checkpoint.restart;
    const { system } = initializeCoupledStreamtubeFromFlow(plan.sourceMach, checkpoint, { tolerance });
    const settings = coupledAssemblyConditions(source.sourceCase, f.input.bodies, system.euler.conditions.lengthScale, { checkpointOptions: f.options });
    require(checkpointDataEqual(settings.normalization, normalization(source)), 'Mach continuation changed the public physical normalization.');
    require((settings.options.transitionMode ?? 'fixed-trip') === system.bl.transitionMode,
      'Mach continuation changed the source transition mode.');
    // The topology system is temporary and is never used by the solver below.
    // A callback may retain it for coefficient extraction without owning any
    // subsequent Newton state or continuation checkpoint.
    const parentMetadata = { model: 'research-streamtube-euler-bl', converged: true,
      mesh: { quality: { valid: true } }, checkpoint: copy(checkpoint), families: copy(checkpoint.families),
      sourceCase: copy(source.sourceCase), ...settings.normalization, solverSettings: copy(solverSettings),
      initialization: copy(source.initialization), ...(source.refinement ? { refinement: copy(source.refinement) } : {}) };
    onPrepared?.({ system, settings: copy(settings), checkpoint: copy(checkpoint), parentResult: parentMetadata,
      ...info(plan.sourceMach), certification: copy(certification) });
  }
  observe(onCheckpoint, checkpoint, info(plan.sourceMach, { kind: 'accepted', stage: 'coupled-source',
    reachedTarget: plan.sourceMach === plan.targetMach }));
  const result = solveCoupledStreamtubeAutomatic(plan.targetMach, { initialCheckpoint: checkpoint, tolerance, stopAtShock,
    stageMaxIterations, maxMachStep, maxStages, maxSubdivisions, maxFirstOrderRecoveries, maxWakeRecoveries, maxBacktracks, blPredictor: predictor, dissipationEnhancement,
    iterationRecovery, maxProgressExtraIterations, includeFlowState: Boolean(onMesh),
    onStage: event => observe(onStage, info(event.mach, { ...event, stage: event.stage ?? 'coupled-mach', startupAttempt: 1 })),
    onIteration: event => observe(onIteration, info(event.mach, { ...event, stage: event.stage ?? 'coupled', startupAttempt: 1 })),
    onMesh: onMesh ? (event, state) => {
      const mesh = { ...event.mesh, ...info(event.mach) };
      mesh.initialization.gridSmoothing = source.initialization?.euler?.gridSmoothing;
      observe(onMesh, mesh, event.iteration.iteration === 0 ? 'initial' : 'solving',
        { ...state, ...info(event.mach), stage: event.stage ?? 'coupled', startupAttempt: 1 });
    } : undefined,
    onCheckpoint: (value, event) => observe(onCheckpoint, value, info(event.mach, { ...event, stage: event.stage ?? 'coupled' })),
  });
  if (observerFailed) throw observerError;
  const actualMach = result.conditions.mach, actualCase = { ...copy(requestedCase), mach: actualMach };
  return attach(result, actualCase, source, {
    solverSettings,
    diagnostics: { certification: copy(certification), coldRecovery: copy(coldRecovery), densityReinitializedOnWarmRestart: false,
      boundaryLayerReinitializedOnWarmRestart: result.continuation.boundaryLayerReinitializedOnWarmRestart,
      initialRedistributionSkippedOnWarmRestart: true } });
}
