// SPDX-License-Identifier: GPL-2.0-or-later
// Simultaneous Euler/grid/all-surface BL/wake Newton with the reconstructed
// ISES density update and grid-maintenance sequence. A single scalar applies
// to the complete direction. Optional retries test physical admissibility
// and, with explicit Armijo selection, residual decrease after maintenance.
// Final convexity remains mandatory.
// The opt-in native event-profile experiment additionally prepares physical
// BL profiles after that scalar step; its changes are reported separately.
import { auditStreamtubeShocks, msesTemporaryMcrit, auditDissipationDamping } from './streamtube-shock-audit.js';
import { createIterationProgress, physicalIterationVector, residualMerit } from './streamtube-iteration-progress.js';
import { createCoupledStreamtubeBody, coupledStreamtubeTripEvents, coupledStreamtubeResult } from './streamtube-coupled.js';
import { proposeCoupledDensityNewton } from './streamtube-density-newton.js';
import { proposeCoupledXfoilBLUpdate, proposeLogarithmicShearUpdate } from './streamtube-coupled-xfoil-update.js';
import { respondToCoupledProjectionGeometry } from './streamtube-coupled-projection-geometry.js';
import { prepareCoupledTransitionProfileTrial } from './streamtube-transition-profile-trial.js';
import { solveSparseDirect } from '../numerics/klu.js';
import { createStreamtubeStationOrdering } from './streamtube-station-ordering.js';
import { solveCoupledLinearSystem } from './streamtube-coupled-linear-solve.js';
import { solveCoupledAlignedSystem } from './streamtube-coupled-aligned-solve.js';
import { redistributeStreamtubeTangentially } from '../geometry/streamtube-tangential-redistribution.js';
import { captureStreamtubeInletFractions, adjustStreamtubeInlets, reparameterizeStreamtubeInlets, dekinkStreamtubeInteriors } from '../geometry/streamtube-grid-maintenance.js';
import { assertConvexStreamtubeGrid } from '../geometry/streamtube-convex-step.js';
import { prepareConvexWakeGridUpdate } from './streamtube-grid-update.js';
import { requireConvexGridUpdate } from './streamtube-grid-update.js';
import { interpolateStreamtubeGridNodes } from '../geometry/streamtube-convex-step.js';
import { requireCoupledResidualDecrease } from './streamtube-iteration-progress.js';
import { relaxCoupledWakeSeed } from './streamtube-coupled-initializer.js';
import { coupledRegularizedDirection, coupledNewtonRayStalled, coupledStepStagnation, coupledConstrainedDirection } from './streamtube-coupled-direction-recovery.js';
import { coupledChangeSnapshot, acceptedCoupledChanges, newtonCoupledChanges, coupledResidualChanges,
  coupledChangeDecision, smallCoupledChanges } from './streamtube-coupled-convergence.js';

// Fine-grid residuals exceed mobile engines' function argument limits.
// Reduce in place: no argument spreading or full-size temporary array.
const maximum = residual => residual.reduce((peak, value) => Math.max(peak, Math.abs(value)), -Infinity);

export function solveCoupledStreamtubeIses(input, { initialEuler, initialBL, reynolds = 1e6, ncrit = 9, tripFractions,
  transitionMode = 'fixed-trip', transitionState, blThermodynamics, hkFloorLinearization, geometryReplay,
  edgeMatching = 'pressure', maxIterations = 12, tolerance = 1e-8, stagnationLimiter = 'listing',
  iterationGeometry = 'convex', stepAcceptance = 'listing', maxBacktracks = 12, blUpdate, projectionGeometry, linearOrdering, eventProfile, shearCoordinate,
  onIteration, onMesh, onCheckpoint, resume, dissipationEnhancement, iterationRecovery, hkProjectionRecovery, directionRecovery, eliminateAmplification,
  wakeGridInitialization = false, maxProgressExtraIterations = 0, convergence } = {}) {
  convergence ??= resume?.continuation?.convergence ?? 'residual';
  if (!['residual', 'mses'].includes(convergence)) throw new Error('Invalid coupled convergence policy.');
  // Restarting or changing conditions requires fresh change evidence.
  let convergenceDecision;
  if (typeof wakeGridInitialization !== 'boolean' || resume && wakeGridInitialization)
    throw new Error('Wake grid initialization requires a fresh coupled solve.');
  if (resume !== undefined) {
    if (input !== undefined || initialEuler !== undefined || initialBL !== undefined || resume?.version !== 1
      || !resume.restart || !resume.continuation) throw new Error('Resume requires one complete ISES checkpoint and no separate initial state.');
    if (hkFloorLinearization !== undefined && hkFloorLinearization !== (resume.restart.options?.hkFloorLinearization ?? 'exact'))
      throw new Error('ISES checkpoint Hk-floor linearization controls do not match.');
    if (geometryReplay !== undefined && geometryReplay !== (resume.restart.options?.geometryReplay ?? 'legacy'))
      throw new Error('ISES checkpoint geometry replay controls do not match.');
    ({ input, options: { reynolds, ncrit, tripFractions, edgeMatching, transitionMode = 'fixed-trip', transitionState, blThermodynamics,
      hkFloorLinearization, geometryReplay }, initialEuler, initialBL } = resume.restart);
    const c = resume.continuation;
    if (eliminateAmplification !== undefined && eliminateAmplification !== (c.eliminateAmplification ?? false))
      throw new Error('ISES checkpoint amplification-update controls do not match.');
    eliminateAmplification = c.eliminateAmplification ?? false;
    const resumedShearCoordinate = c.shearCoordinate === undefined ? 'linear' : c.shearCoordinate;
    if (shearCoordinate !== undefined && shearCoordinate !== resumedShearCoordinate)
      throw new Error('ISES checkpoint shear-coordinate controls do not match.');
    shearCoordinate = resumedShearCoordinate;
    if (linearOrdering !== undefined && linearOrdering !== (c.linearOrdering ?? 'auto'))
      throw new Error('ISES checkpoint linear ordering controls do not match.');
    linearOrdering = c.linearOrdering ?? 'auto';
    if (blUpdate !== undefined && blUpdate !== (c.blUpdate ?? 'giles')) throw new Error('ISES checkpoint BL update controls do not match.');
    blUpdate = c.blUpdate ?? 'giles';
    if (projectionGeometry !== undefined && projectionGeometry !== (c.projectionGeometry ?? 'fixed'))
      throw new Error('ISES checkpoint projection geometry controls do not match.');
    projectionGeometry = c.projectionGeometry ?? 'fixed';
    const resumedEventProfile = c.eventProfile === undefined ? 'none' : c.eventProfile;
    if (eventProfile !== undefined && eventProfile !== resumedEventProfile)
      throw new Error('ISES checkpoint event profile controls do not match.');
    eventProfile = resumedEventProfile;
    if (iterationGeometry !== c.iterationGeometry || stepAcceptance !== c.stepAcceptance || stagnationLimiter !== c.stagnationLimiter)
      throw new Error('ISES checkpoint update controls do not match.');
  }
  blUpdate ??= 'giles';
  if (shearCoordinate === undefined) shearCoordinate = 'linear';
  projectionGeometry ??= 'fixed';
  if (eventProfile === undefined) eventProfile = 'none';
  eliminateAmplification ??= blUpdate === 'xfoil' && transitionMode === 'automatic' && eventProfile === 'none';
  if (typeof eliminateAmplification !== 'boolean' || eliminateAmplification
    && (blUpdate !== 'xfoil' || transitionMode !== 'automatic' || eventProfile !== 'none'))
    throw new Error('Invalid coupled amplification-elimination controls.');
  // Archived checkpoints keep their original update policy. New monotone
  // XFOIL updates may recover a proven surface-projection descent failure.
  hkProjectionRecovery ??= resume ? resume.continuation.hkProjectionRecovery ?? false
    : blUpdate === 'xfoil' && eventProfile === 'none' && ['armijo', 'event-armijo'].includes(stepAcceptance);
  directionRecovery ??= resume ? resume.continuation.directionRecovery ?? false
    : blUpdate === 'xfoil' && eventProfile === 'none' && ['armijo', 'event-armijo'].includes(stepAcceptance);
  if (typeof directionRecovery !== 'boolean' || directionRecovery
    && (blUpdate !== 'xfoil' || eventProfile !== 'none' || !['armijo', 'event-armijo'].includes(stepAcceptance)))
    throw new Error('Invalid coupled direction recovery controls.');
  if (typeof hkProjectionRecovery !== 'boolean' || hkProjectionRecovery
    && (blUpdate !== 'xfoil' || eventProfile !== 'none' || !['armijo', 'event-armijo'].includes(stepAcceptance)))
    throw new Error('Invalid coupled Hk projection recovery controls.');
  // Old checkpoints resolve omission to their original automatic ordering
  // above. Only a new solve selects guarded physical-station ordering.
  linearOrdering ??= 'station-auto';
  iterationRecovery ??= resume ? resume.continuation.iterationRecovery === true : stepAcceptance === 'admissible';
  if (typeof iterationRecovery !== 'boolean' || !Number.isInteger(maxProgressExtraIterations)
    || maxProgressExtraIterations < 0 || maxProgressExtraIterations > 20)
    throw new Error('Invalid coupled progress recovery controls.');
  let wakeCoordinateRecovery = resume?.continuation.wakeCoordinateRecovery ?? false;
  let stationFallback = resume?.continuation.stationFallback === undefined ? false : resume.continuation.stationFallback;
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0
    || !['listing', 'prose'].includes(stagnationLimiter) || !['convex', 'ises-sampled'].includes(iterationGeometry)
    || !['listing', 'admissible', 'armijo', 'event-armijo'].includes(stepAcceptance) || !['giles', 'xfoil'].includes(blUpdate)
    || !['fixed', 'boundary-increment'].includes(projectionGeometry) || !['auto', 'station', 'station-auto', 'aligned-auto'].includes(linearOrdering)
    || typeof stationFallback !== 'boolean' || stationFallback && linearOrdering !== 'station-auto'
    || typeof wakeCoordinateRecovery !== 'boolean'
    || projectionGeometry === 'boundary-increment' && blUpdate !== 'xfoil'
    || !['linear', 'logarithmic'].includes(shearCoordinate) || shearCoordinate === 'logarithmic' && blUpdate !== 'xfoil'
    || !['none', 'xfoil-mrchdu'].includes(eventProfile)
    || eventProfile === 'xfoil-mrchdu' && (transitionMode !== 'automatic' || blUpdate !== 'xfoil' || projectionGeometry !== 'boundary-increment')
    || !Number.isInteger(maxBacktracks) || maxBacktracks < 0 || maxBacktracks > 20)
    throw new Error('Invalid coupled ISES iteration controls.');
  const savedDissipation = resume?.continuation.dissipationEnhancement;
  dissipationEnhancement ??= savedDissipation !== undefined;
  if (typeof dissipationEnhancement !== 'boolean' || savedDissipation && !dissipationEnhancement)
    throw new Error('Invalid or conflicting dissipation enhancement policy.');
  const targetMcrit = savedDissipation?.targetMcrit ?? input.upwind?.mcrit;
  let previousDensityChange = savedDissipation?.previousDensityChange ?? 0;
  if (dissipationEnhancement) {
    if (!input.upwind || !['hybrid', 'momentum'].includes(input.streamwiseMode)
      || savedDissipation && savedDissipation.method !== 'mses')
      throw new Error('MSES dissipation enhancement requires explicit momentum or hybrid upwinding.');
    msesTemporaryMcrit(targetMcrit, previousDensityChange);
  }
  const geometryDomain = iterationGeometry === 'convex' ? 'convex' : 'positive-simple';
  const preserveConvexity = stepAcceptance !== 'listing';
  const requireConvex = preserveConvexity || iterationGeometry === 'convex';
  const propose = blUpdate === 'xfoil' ? (shearCoordinate === 'logarithmic' ? proposeLogarithmicShearUpdate : proposeCoupledXfoilBLUpdate) : proposeCoupledDensityNewton;
  if (projectionGeometry === 'boundary-increment' && input.wakeGeometry !== 'independent-banks')
    throw new Error('Projected boundary response currently requires independent wake banks.');
  if (input.stagnationMotion !== undefined && !['walls-only', 'interpolated'].includes(input.stagnationMotion)) throw new Error('Unknown coupled ISES stagnation motion.');
  if (input.geometryDomain !== undefined && input.geometryDomain !== geometryDomain) throw new Error('Conflicting coupled ISES geometry domains.');
  input = { ...input, stagnationMotion: input.stagnationMotion ?? 'walls-only', normalStencil: input.normalStencil ?? 'body-stations', geometryDomain };
  const coupledOptions = { reynolds, ncrit, tripFractions, edgeMatching, ...(transitionMode === 'automatic' ? { transitionMode } : {}),
    ...(blThermodynamics === undefined ? {} : { blThermodynamics }),
    ...(hkFloorLinearization === undefined || hkFloorLinearization === 'exact' ? {} : { hkFloorLinearization }),
    ...(geometryReplay === undefined || geometryReplay === 'legacy' ? {} : { geometryReplay }) };
  let system = createCoupledStreamtubeBody(input, { ...coupledOptions, initialEuler, initialBL, transitionState });
  // Pressure matching uses dynamic pressure units. Use the same units for
  // Euler pressure rows in its merit; their pInf normalization otherwise
  // hides the pressure balance at low Mach. Root tolerances are unchanged.
  const pressureRows = new Set(['streamwise', 'internalPressure', 'farfieldPressure', 'cutPressure', 'trailingKutta', 'leadingKutta']);
  const meritWeights = edgeMatching === 'pressure' ? Float64Array.from(system.initial, (_, i) =>
    i < system.ne && pressureRows.has(system.euler.layout.rows[i].kind) ? system.euler.conditions.pressureScale : 1) : undefined;
  const meritControls = step => ({ step, tolerance, ...(meritWeights ? { weights: meritWeights } : {}) });
  if (eventProfile === 'xfoil-mrchdu' && (!system.bl.trips.every(pair => pair.every(x => x === 1))
    || system.bl.snapshotActive().some(j => j <= 0) || system.euler.layout.wakeDisplacementMotion !== 'te-center'))
    throw new Error('Native event-profile trial requires terminal trips, resolved transitions and the TE-center wake chart.');
  const optionsAtState = () => ({ ...coupledOptions, ...(transitionMode === 'automatic' ? { transitionState: system.bl.snapshotActive() } : {}) });
  if (!system.euler.layout.densityCount) throw new Error('Coupled ISES density updates require compressible equations.');
  let state = system.initial.slice();
  let value = system.admissibleValue(state, { requireConvex });
  if (!value) throw new Error('Inadmissible coupled ISES initial state.');
  let reason = 'iteration limit', failed = false, lastRejectedStep = null;
  if (resume && Object.keys(value.families).some(k => value.families[k] !== resume.families?.[k]))
    throw new Error('ISES checkpoint residual does not replay exactly.');
  const fractions = resume ? structuredClone(resume.continuation.fractions) : captureStreamtubeInletFractions(value.outer.nodes, system.euler.layout.bodies);
  if (!Array.isArray(fractions) || fractions.length !== input.bodies.length || fractions.some((f, body) =>
    !Array.isArray(f) || f.length !== input.bodies[body].leadingIndex + 1 || f[0] !== 0 || f.at(-1) !== 1
    || !f.every(Number.isFinite) || f.some((v, i) => i && v <= f[i - 1]))) throw new Error('Invalid ISES checkpoint inlet fractions.');
  const adopt = (x, nodes, undisplacedNodes) => {
    // Build an independent chart, then validate the exact displaced walls
    // and normal wake gaps through adoptGeometry's reversible restore.
    // No trial may change the accepted chart or silently move a material trip.
    if (preserveConvexity) assertConvexStreamtubeGrid(nodes);
    const candidate = createCoupledStreamtubeBody(input, { ...optionsAtState(), initialBL: x.slice(system.ne),
      initialEuler: { x: x.slice(0, system.ne), undisplacedNodes } });
    const state = candidate.initial.slice();
    candidate.euler.setDisplacement(candidate.bl.thicknesses(state.subarray(candidate.ne)));
    state.set(candidate.euler.adoptGeometry(state.subarray(0, candidate.ne), nodes));
    // Amplification transport is causal: solve N for the trial thickness,
    // velocity and FINAL station geometry. The simultaneous Newton system
    // supplies its tangent. Re-integrating before grid maintenance would
    // leave the accepted N profile inconsistent with its station lengths.
    // Keep all governing rows and physical BL unknowns; this correction
    // still has to pass the complete maintained-residual acceptance test.
    let amplification;
    if (eliminateAmplification) {
      const phase = candidate.bl.snapshotActive();
      const update = candidate.bl.updateActive(state.subarray(candidate.ne), state.subarray(0, candidate.ne),
        { reconcileAmplificationSurfaces: candidate.bl.surfaces.map((_, k) => k) });
      amplification = { phaseChanged: candidate.bl.snapshotActive().some((v, k) => v !== phase[k]),
        changes: update.changes };
    }
    let admissibility;
    const value = candidate.admissibleValue(state, { requireConvex, onFailure: failure => { admissibility = failure; } });
    if (!value) {
      const detail = admissibility?.kind === 'bl-domain'
        ? ` BL station ${admissibility.station.id} (body ${admissibility.station.body}, ${admissibility.station.side ?? admissibility.station.kind}, i=${admissibility.station.i}): ${admissibility.failedConstraints.join(', ')}; shape=${admissibility.shape}, enthalpy=${admissibility.enthalpy}.`
        : admissibility?.kind === 'grid-convexity'
          ? ` Cell ${admissibility.firstCell.cell} (group ${admissibility.firstCell.group}, interval ${admissibility.firstCell.interval}, tube ${admissibility.firstCell.tube}) fails final convexity.`
          : admissibility ? ` ${admissibility.stage}: ${admissibility.message}` : '';
      const error = new Error(`Maintained coupled state fails physical/grid admissibility.${detail}`);
      if (admissibility) error.coupledAdmissibility = admissibility;
      if (admissibility?.code !== undefined) error.code = admissibility.code;
      if (admissibility?.diagnostics !== undefined) error.diagnostics = structuredClone(admissibility.diagnostics);
      throw error;
    }
    return { system: candidate, state, value, amplification };
  };
  // Reconstruct the same physical state when changing ONLY dissipation.
  // Candidate ownership keeps a failed parameter change transactional.
  const stateAtMcrit = mcrit => {
    if (input.upwind.mcrit === mcrit) return { input, system, state, value };
    const nextInput = { ...input, upwind: { ...input.upwind, mcrit } };
    const candidate = createCoupledStreamtubeBody(nextInput, { ...optionsAtState(),
      initialEuler: { x: state.slice(0, system.ne), nodes: value.outer.nodes,
        undisplacedNodes: value.outer.undisplacedNodes }, initialBL: state.slice(system.ne) });
    const nextState = candidate.initial.slice(), nextValue = candidate.admissibleValue(nextState, { requireConvex });
    if (!nextValue) throw new Error('Dissipation change fails physical/grid admissibility.');
    if (nextState.length !== state.length || nextState.some((v, i) => v !== state[i]))
      throw new Error('Dissipation change altered the packed physical state.');
    return { input: nextInput, system: candidate, state: nextState, value: nextValue };
  };
  const changeMcrit = mcrit => { ({ input, system, state, value } = stateAtMcrit(mcrit)); };
  const redistribute = (nodes, correctionScale = 1) => {
    const passages = [], moved = nodes.map((group, g) => {
      const r = redistributeStreamtubeTangentially(group, { referenceBank: g === 0 ? group[0].length - 1 : 0,
        fixedBanks: [g !== 0, g !== nodes.length - 1], correctionScale, quadratureDomain: iterationGeometry === 'convex' ? 'convex' : 'sampled-positive' });
      passages.push({ group: g, referenceBank: r.referenceBank, fixedBanks: r.fixedBanks, pairs: r.solution.pairs,
        coordinateRelativeResidual: r.solution.relativeResidual, maxDisplacement: r.maxDisplacement,
        ...(correctionScale === 1 ? {} : { correctionScale }) });
      return r.nodes;
    });
    return { nodes: moved, passages };
  };
  const initialRedistribution = { beforeFamilies: value.families, accepted: false };
  try {
    if (resume) Object.assign(initialRedistribution, { accepted: true, resumed: true, afterFamilies: value.families, passages: [] });
    else {
      // As in the inviscid driver, preserve the listed full SMOVE first.
      // An admissible BL seed can fail its finite coordinate correction;
      // retry that correction from the unchanged full coupled state.
      let correctionScale = 1;
      const rejections = [];
      for (let trial = 0; ; trial++) {
        let stage = 'SMOVE';
        try {
          const r = redistribute(value.outer.nodes, correctionScale);
          stage = 'admissibility';
          const candidate = adopt(state, r.nodes, value.outer.undisplacedNodes);
          ({ system, state, value } = candidate);
          Object.assign(initialRedistribution, { accepted: true, afterFamilies: value.families, passages: r.passages,
            ...(rejections.length ? { correctionScale, backtracks: rejections.length, rejections } : {}) });
          break;
        } catch (error) {
          if (stepAcceptance === 'listing') throw error;
          rejections.push({ correctionScale, stage, message: error.message,
            ...(error.coupledAdmissibility ? { admissibility: error.coupledAdmissibility } : {}) });
          Object.assign(initialRedistribution, { correctionScale, backtracks: rejections.length - 1, rejections });
          if (trial >= maxBacktracks) throw error;
          correctionScale *= .5;
        }
      }
    }
  } catch (error) {
    reason = `Coupled ISES initial redistribution rejected: ${error.message}`; failed = true; initialRedistribution.rejection = error.message;
    if (error.coupledAdmissibility) initialRedistribution.admissibility = error.coupledAdmissibility;
  }
  if (wakeGridInitialization && !failed && blUpdate === 'xfoil' && stepAcceptance === 'event-armijo' && transitionMode === 'automatic') {
    const relaxed = relaxCoupledWakeSeed(input, optionsAtState(), system, state, value);
    if (relaxed) {
      initialRedistribution.wakeSmoothing = relaxed.report;
      if (relaxed.report.accepted) ({ system, state, value } = relaxed);
    }
  }
  let lastRedistributedStagnation = resume ? resume.continuation.lastRedistributedStagnation.slice() : value.outer.stagnation.slice();
  let preferredOrdering = resume?.continuation.preferredOrdering ?? 'amd', pivotTolerance = resume?.continuation.pivotTolerance ?? .001;
  if (lastRedistributedStagnation.length !== input.bodies.length || !lastRedistributedStagnation.every(Number.isFinite)
    || !['amd', 'colamd'].includes(preferredOrdering) || ![.001, 1].includes(pivotTolerance)) throw new Error('Invalid ISES checkpoint continuation state.');
  const history = [], linearDiagnostics = { solves: 0, maxRelativeResidual: 0, refinements: 0, pivotRecoveries: 0, iterations: [] };
  const progress = createIterationProgress();
  const phaseKey = () => JSON.stringify(system.bl.snapshotActive());
  const seedProgress = () => progress.seed(physicalIterationVector(state, value.outer.nodes), phaseKey(), input.upwind?.mcrit);
  if (iterationRecovery) seedProgress();
  let recovery = null, recoveries = 0, progressStop, plateauRecoveryDisabled = false;
  const checkpoint = () => ({ version: 1, families: { ...value.families },
    ...(convergenceDecision?.converged ? { convergence: structuredClone(convergenceDecision) } : {}),
    restart: { input, options: optionsAtState(), initialEuler: { x: state.slice(0, system.ne), nodes: value.outer.nodes,
      undisplacedNodes: value.outer.undisplacedNodes }, initialBL: state.slice(system.ne) },
    continuation: { ...(iterationRecovery ? { iterationRecovery: true } : {}),
      ...(convergence === 'mses' ? { convergence } : {}),
      ...(eliminateAmplification ? { eliminateAmplification: true } : {}),
      ...(wakeCoordinateRecovery ? { wakeCoordinateRecovery: true } : {}),
      ...(hkProjectionRecovery ? { hkProjectionRecovery: true } : {}),
      ...(directionRecovery ? { directionRecovery: true } : {}),
      ...(dissipationEnhancement ? { dissipationEnhancement: { method: 'mses', targetMcrit, previousDensityChange } } : {}),
      fractions, lastRedistributedStagnation: lastRedistributedStagnation.slice(), preferredOrdering, pivotTolerance,
      iterationGeometry, stepAcceptance, stagnationLimiter, ...(blUpdate === 'xfoil' ? { blUpdate } : {}),
      ...(linearOrdering !== 'auto' ? { linearOrdering } : {}),
      ...(stationFallback ? { stationFallback: true } : {}),
      ...(eventProfile !== 'none' ? { eventProfile } : {}),
      ...(shearCoordinate !== 'linear' ? { shearCoordinate } : {}),
      ...(projectionGeometry === 'boundary-increment' ? { projectionGeometry } : {}) } });
  const report = (iteration, details) => {
    const entry = { ...details, iteration, residual: maximum(value.residual), ...value.families,
      shearCoordinate, hkFloorLinearization: coupledOptions.hkFloorLinearization ?? 'exact' };
    if (value.outer.diagnostics.maxMach >= .9) entry.shockAudit = auditStreamtubeShocks(value.outer, system.euler.conditions);
    if (iterationRecovery) {
      entry.residualContext = { mcrit: input.upwind?.mcrit, transition: system.bl.snapshotActive(),
        equationsChanged: history.length > 0 && history.at(-1).residualContext?.mcrit !== input.upwind?.mcrit,
        ...(iteration === 0 && resume ? { monitorRestarted: true } : {}) };
      // Observe prescribed equations periodically, without committing their chart.
      if (!dissipationEnhancement || input.upwind.mcrit === targetMcrit)
        entry.prescribedResidual = { iteration, residual: entry.residual, ...value.families };
      else if (iteration % 5 === 0 || details.progress?.cause) {
        try { const prescribed = stateAtMcrit(targetMcrit).value;
          entry.prescribedResidual = { iteration, residual: maximum(prescribed.residual), ...prescribed.families }; }
        catch (error) { entry.prescribedResidual = { iteration, unavailable: error.message }; }
      }
    }
    history.push(entry); onIteration?.(entry);
    if (initialRedistribution.accepted && onCheckpoint) onCheckpoint(checkpoint(), { history, linearDiagnostics, initialRedistribution });
    onMesh?.({ system: system.euler, nodes: value.outer.nodes, flow: value.outer, coupledFamilies: value.families, iteration: { ...entry } });
  };
  report(0, { step: 0 });
  for (let iteration = 1; !failed && !convergenceDecision?.converged && iteration <= maxIterations + maxProgressExtraIterations && (maximum(value.residual) > tolerance || dissipationEnhancement && input.upwind.mcrit !== targetMcrit); iteration++) {
    if (iteration > maxIterations && (!iterationRecovery || recovery || !progress.canExtend())) break;
    let proposal, stage = 'linear solve', accepted, rejections = [], recoveryOutcome, transitionUnchanged = false;
    let iterationDissipation;
    if (dissipationEnhancement) {
      try {
        const mcrit = maximum(value.residual) <= tolerance ? targetMcrit : recovery?.mcrit ?? msesTemporaryMcrit(targetMcrit, previousDensityChange);
        changeMcrit(mcrit);
        iterationDissipation = { method: 'mses', mcrit, targetMcrit, previousDensityChange };
        if (mcrit === targetMcrit && maximum(value.residual) <= tolerance) break;
      } catch (error) {
        reason = error.message; failed = true;
        lastRejectedStep = { stage: 'dissipation enhancement', message: error.message, rejections: [] };
        break;
      }
    }
    const events = coupledStreamtubeTripEvents(system, { blUpdate }), phase = events.snapshot();
    const beforeChange = convergence === 'mses' ? coupledChangeSnapshot(system, state, value) : null;
    let newtonChanges;
    const beforeMerit = residualMerit(value.residual);
    let projectionDescentFailure;
    const diagnoseProjectionDescent = () => {
      if (!hkProjectionRecovery || coupledOptions.hkFloorLinearization === 'native' || !transitionUnchanged
        || !proposal?.projection?.displacementChanges.some(p => p.kind === 'surface' && p.correction > 0)) return;
      const raw = proposal.x.slice();
      for (const p of proposal.projection.displacementChanges) raw[system.ne + 4 * p.id + 2] = p.beforeDeltaStar;
      for (const p of proposal.projection.auxiliaryChanges) raw[system.ne + 4 * p.id] = p.before;
      try {
        const rawValue = system.admissibleValue(raw, { requireConvex });
        if (rawValue && residualMerit(rawValue.residual) < beforeMerit) {
          const projectedValue = system.admissibleValue(proposal.x, { requireConvex });
          if (projectedValue && residualMerit(projectedValue.residual) >= beforeMerit)
            return { rawValue, projectedValue };
        }
      } finally { system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne))); }
    };
    try {
      const matrixStarted = ['station-auto', 'aligned-auto'].includes(linearOrdering) ? performance.now() : undefined;
      const matrix = system.jacobian(state);
      const jacobianMilliseconds = matrixStarted === undefined ? undefined : performance.now() - matrixStarted;
      let ordering, linear, stationPolicy, rejectedLinearSolve, initialDirection;
      try {
        if (linearOrdering === 'aligned-auto') {
          ({ ordering, linear, stationPolicy } = solveCoupledAlignedSystem(matrix, value.residual.map(v => -v), {
            layout: system.euler.layout, stations: system.bl.stations, preferredOrdering, pivotTolerance,
          }));
        } else if (linearOrdering === 'station-auto') {
          ({ ordering, linear, stationPolicy } = solveCoupledLinearSystem(matrix, value.residual.map(v => -v), {
            layout: system.euler.layout, stations: system.bl.stations, preferredOrdering, pivotTolerance,
            mode: linearOrdering, stationEnabled: !stationFallback,
          }));
          if (stationPolicy.recommendation === 'auto') stationFallback = true;
        } else {
          // Preserve explicit research ordering and archived automatic solves.
          ordering = linearOrdering === 'station' ? createStreamtubeStationOrdering({ matrix,
            layout: system.euler.layout, stations: system.bl.stations }) : null;
          linear = solveSparseDirect(matrix, value.residual.map(v => -v), ordering
            ? { ordering: 'given', rowPermutation: ordering.Puser, columnPermutation: ordering.Quser,
              btf: false, pivotTolerance, pivotFallback: false }
            : { preferredOrdering, pivotTolerance });
        }
      } catch (error) {
        if (!directionRecovery || error.code !== 'KLU_RESIDUAL_LIMIT') throw error;
        // An uncertified Newton solve is never used. A constrained residual
        // gradient needs no inverse Jacobian and must pass the same nonlinear
        // gates; this can move the next linearization away from singularity.
        try { initialDirection = coupledConstrainedDirection(system, state, matrix, value.residual, { weights: meritWeights }); }
        catch { throw error; }
        if (!initialDirection) throw error;
        rejectedLinearSolve = { code: error.code, diagnostics: error.diagnostics, attempts: error.attempts };
        linear = { x: initialDirection.direction };
      }
      if (rejectedLinearSolve) {
        linearDiagnostics.iterations.push({ iteration, rejectedLinearSolve,
          directionRecovery: initialDirection.diagnostics });
      } else {
        if (!ordering) preferredOrdering = linear.ordering ?? preferredOrdering;
        pivotTolerance = linear.pivotTolerance ?? pivotTolerance;
        linearDiagnostics.solves++; linearDiagnostics.maxRelativeResidual = Math.max(linearDiagnostics.maxRelativeResidual, linear.relativeResidual);
        linearDiagnostics.refinements += linear.refinements;
        if (linear.pivotTolerance === 1 && linear.attempts.some(a => a.pivotTolerance < 1)) linearDiagnostics.pivotRecoveries++;
        linearDiagnostics.iterations.push({ iteration, ordering: linear.ordering, pivotTolerance: linear.pivotTolerance, attempts: linear.attempts,
          ...(stationPolicy ? { stationPolicy, jacobianMilliseconds } : {}),
          ...(ordering ? { stationOrdering: { version: 1, matchedNonzeroDiagonals: ordering.diagnostics.n,
            borderPairs: ordering.diagnostics.borderPairCount, maximumLocalStationMismatch: ordering.diagnostics.maximumLocalStationMismatch,
            btf: false } } : {}) });
      }
      if (beforeChange && !rejectedLinearSolve && smallCoupledChanges(coupledResidualChanges(system, value.residual)))
        newtonChanges = newtonCoupledChanges(system, state, linear.x, value);
      let maximumStep = recovery ? .5 : 1, redistributionBacktracks = 0, geometryRedistribution = false;
      let automaticFallback, transitionFallback, initialProposalStep;
      let directionIndex = initialDirection ? 3 : 0,
        directionDetails = initialDirection ? { ...initialDirection.diagnostics, cause: 'uncertified-newton-solve' } : undefined,
        stagnationFallback, constrainedDirectionTried = Boolean(initialDirection);
      const nextDirection = cause => {
        let alternative;
        const strengths = [.5, .1, 2];
        while (!alternative && directionIndex < strengths.length) {
          try { alternative = coupledRegularizedDirection(matrix, value.residual, system.ne,
            strengths[directionIndex++], { weights: meritWeights, preferredOrdering, pivotTolerance }); }
          catch { /* A failed factorization cannot supply a search direction. */ }
        }
        if (!alternative && !constrainedDirectionTried) {
          constrainedDirectionTried = true;
          try {
            events.restore(phase);
            alternative = coupledConstrainedDirection(system, state, matrix, value.residual, { weights: meritWeights });
          } catch { /* A failed constraint projection cannot supply a direction. */ }
        }
        if (!alternative) return false;
        directionDetails = { ...alternative.diagnostics, cause };
        if (alternative.linear) {
          linear = alternative.linear;
          linearDiagnostics.solves++;
          linearDiagnostics.refinements += linear.refinements;
          linearDiagnostics.maxRelativeResidual = Math.max(linearDiagnostics.maxRelativeResidual, linear.relativeResidual);
          linearDiagnostics.iterations.push({ iteration, ordering: linear.ordering,
            pivotTolerance: linear.pivotTolerance, attempts: linear.attempts, directionRecovery: directionDetails });
        } else linear = { ...linear, x: alternative.direction };
        maximumStep = 1; initialProposalStep = undefined; geometryRedistribution = false;
        redistributionBacktracks = 0; automaticFallback = undefined; transitionFallback = undefined;
        return true;
      };
      for (let trial = 0; ; trial++) {
        let redistributionAttempted = false, repairOnly = false, profile;
        try {
          proposal = undefined;
          transitionUnchanged = false;
          events.restore(phase);
          stage = 'coupled density proposal'; proposal = propose(system, state, linear.x, { stagnationLimiter, maximumStep });
          if (directionDetails) proposal.directionRecovery = directionDetails;
          initialProposalStep ??= proposal.step;
          system.euler.setDisplacement(system.bl.thicknesses(proposal.x.subarray(system.ne)));
          let event, decoded, response;
          if (eventProfile === 'xfoil-mrchdu') {
            stage = 'decode displaced grid'; decoded = system.euler.decode(proposal.x.subarray(0, system.ne));
            stage = 'projection grid response'; response = respondToCoupledProjectionGeometry(system, proposal.x, proposal.projection, decoded);
            stage = 'native event profile'; profile = prepareCoupledTransitionProfileTrial(system, proposal.x, { decoded, baseNodes: response.nodes });
            if (profile.diagnostics.active) {
              proposal.x = profile.x; decoded = profile.decoded;
              system.euler.setDisplacement(system.bl.thicknesses(proposal.x.subarray(system.ne)));
            }
            stage = 'material-trip transfer'; event = events.prepare(proposal.x, state);
          } else {
            stage = 'material-trip transfer'; event = events.prepare(proposal.x, state);
            stage = 'decode displaced grid'; decoded = system.euler.decode(proposal.x.subarray(0, system.ne));
            // A thickness projection changes the displaced boundary after the
            // common Newton increment. Respond on trial coordinates only; the
            // independent candidate below owns geometry adoption and domain checks.
            stage = 'projection grid response';
            response = projectionGeometry === 'boundary-increment'
              ? respondToCoupledProjectionGeometry(system, proposal.x, proposal.projection, decoded) : null;
          }
          let proposedNodes = profile?.diagnostics.active ? profile.nodes : response?.nodes ?? decoded.nodes;
          transitionUnchanged = !event.changed;
          stage = 'Newton grid step';
          if (preserveConvexity) {
            // Compare every trial with the same initial Newton-step scale.
            // As backtracking shrinks the proposal, the fraction at the same
            // physical grid contact grows. It must not switch off an ongoing
            // wake recovery merely because its denominator became smaller.
            // Once the restricted chart needs pairing, retain that coordinate
            // policy across Newton updates and checkpoint resumes. Returning
            // immediately to the raw bound can trap the repaired wake again.
            const minimumGridFraction = wakeCoordinateRecovery ? 1
              : Math.min(1, initialProposalStep * 2 ** -maxBacktracks / proposal.step);
            let gridTrial;
            try {
              gridTrial = prepareConvexWakeGridUpdate(value.outer.nodes, proposedNodes,
                system.euler.layout, decoded.allocation, minimumGridFraction);
            } catch (error) {
              const cell = error.diagnostics?.cell;
              const inletLimited = error.code === 'streamtube-grid-step' && cell
                && system.euler.layout.bodies.some((b, body) => cell.i < b.leadingIndex
                  && (cell.group === body && cell.tube === system.euler.layout.tubes[body] - 1
                    || cell.group === body + 1 && cell.tube === 0));
              if (!inletLimited) throw error;
              // As in the inviscid update, inlet re-spacing belongs to the
              // complete trial. A raw normal-coordinate corner may fold even
              // when that complete movement remains convex and decreases merit.
              const repaired = reparameterizeStreamtubeInlets(proposedNodes, system.euler.layout.bodies, fractions);
              try { requireConvexGridUpdate(value.outer.nodes, repaired.nodes); }
              catch (repairError) {
                if (repairError.code === 'streamtube-grid-step') repairError.code = 'streamtube-inlet-grid-step';
                throw repairError;
              }
              gridTrial = { nodes: repaired.nodes };
            }
            proposedNodes = gridTrial.nodes;
            if (gridTrial.correction) {
              proposal.wakeCoordinateRepair = gridTrial.correction;
              wakeCoordinateRecovery = true;
            }
          }
          stage = 'inlet adjustment'; const inlet = adjustStreamtubeInlets(proposedNodes, system.euler.layout.bodies, fractions, { preserveConvexity });
          stage = 'DEKINK'; const dekink = dekinkStreamtubeInteriors(inlet.nodes, { preserveConvexity });
          const triggeredBodies = system.euler.layout.bodies.flatMap((b, body) => {
            const lower = dekink.nodes[body][b.leadingIndex + 1].at(-1), upper = dekink.nodes[body + 1][b.leadingIndex + 1][0];
            const spacing = .5 * Math.hypot(upper.x - lower.x, upper.y - lower.y);
            return Math.abs(decoded.stagnation[body] - lastRedistributedStagnation[body]) > .5 * spacing ? [body] : [];
          });
          // Raw Newton/projection rejections do not reject a coordinate
          // correction that has not run. Preserve full SMOVE on its first
          // actual attempt, then reduce only rejected redistribution trials.
          // A convexity-limited direction also requests the existing SMOVE
          // coordinate correction. Thin wall tubes can lose tangential
          // alignment before the stagnation-motion threshold is reached.
          // Only the reduced, convex raw trial reaches this operation; the
          // complete maintained Euler/BL state still passes every final gate.
          stage = 'SMOVE'; redistributionAttempted = triggeredBodies.length > 0 || geometryRedistribution;
          repairOnly = geometryRedistribution && triggeredBodies.length === 0;
          let r = redistributionAttempted ? redistribute(dekink.nodes, 2 ** -redistributionBacktracks) : { nodes: dekink.nodes, passages: [] };
          stage = 'admissibility'; let candidate;
          if (['armijo', 'event-armijo'].includes(stepAcceptance)
            && !(stepAcceptance === 'event-armijo' && event.changed)) {
            // DEKINK/SMOVE is a finite coordinate correction, separate from
            // Newton. Damping Newton cannot eliminate its zero-step jump.
            // First prove descent of the raw physical trial, then backtrack
            // maintenance independently. Every candidate is fully admitted.
            try {
              candidate = adopt(proposal.x, r.nodes, decoded.undisplacedNodes);
              stage = 'residual decrease';
              requireCoupledResidualDecrease(value.residual, candidate.value.residual, meritControls(proposal.step));
            } catch {
              stage = 'admissibility';
              candidate = adopt(proposal.x, proposedNodes, decoded.undisplacedNodes);
              stage = 'residual decrease';
              requireCoupledResidualDecrease(value.residual, candidate.value.residual, meritControls(proposal.step));
              let correctionScale = 0;
              for (let k = 1; k <= maxBacktracks; k++) {
                try {
                  const moved = adopt(proposal.x, interpolateStreamtubeGridNodes(proposedNodes, r.nodes, 2 ** -k), decoded.undisplacedNodes);
                  requireCoupledResidualDecrease(value.residual, moved.value.residual, meritControls(proposal.step));
                  candidate = moved; correctionScale = 2 ** -k; break;
                } catch { /* Retain the fully admitted raw trial. */ }
              }
              r = { ...r, nodes: candidate.value.outer.nodes, correctionScale };
            }
          } else candidate = adopt(proposal.x, r.nodes, decoded.undisplacedNodes);
          if (candidate.amplification?.phaseChanged) {
            event.changed = true;
            event.changes = [...(event.changes ?? []), ...candidate.amplification.changes];
          }
          const admitted = { candidate, event, inlet, dekink, triggeredBodies, r, redistributionAttempted,
            geometryRedistribution, response: response?.diagnostics, profile: profile?.diagnostics };
          // Optional recovery must not turn the ordinary admissible policy
          // into a mandatory monotone solve. Preserve its first fully checked
          // candidate before testing the extra merit condition.
          if (recovery && stepAcceptance === 'admissible' && !automaticFallback)
            automaticFallback = { accepted: admitted, proposal };
          // First seek descent, including across a transition event. If the
          // search cannot decrease merit, a fully admitted interval change
          // can cross the discrete N-to-shear equation boundary. Only then
          // restart descent on its new topology; no physical gate is waived.
          // Retain the smallest admitted event step encountered by the search.
          // Reusing its first, larger trial discards the closer boundary
          // crossing and can restart Newton far from the converged branch.
          if (stepAcceptance === 'event-armijo' && event.changed
            && (!transitionFallback || proposal.step < transitionFallback.proposal.step))
            transitionFallback = { accepted: { ...admitted, transitionEventAcceptance: true }, proposal };
          let residualDecrease;
          if (stepAcceptance === 'armijo' || stepAcceptance === 'event-armijo' || recovery && !event.changed) {
            // Compare the fully maintained, independently admitted candidate.
            // Reuse the same global direction on failure. A grid-repair
            // request from a larger trial is not a stagnation-motion event.
            stage = 'residual decrease';
            residualDecrease = requireCoupledResidualDecrease(value.residual, candidate.value.residual,
              meritControls(proposal.step));
          }
          accepted = { candidate, event, inlet, dekink, triggeredBodies, r, redistributionAttempted,
            geometryRedistribution, response: response?.diagnostics, profile: profile?.diagnostics, residualDecrease };
          if (directionRecovery && !directionDetails && coupledNewtonRayStalled(history, proposal, event, residualDecrease)) {
            stagnationFallback = { accepted, proposal };
            if (nextDirection('stalled-newton-ray')) { accepted = undefined; trial = -1; continue; }
          }
          // An alternate ray that makes the same negligible progress has
          // not resolved the stall. Keep the best admitted endpoint while
          // trying the remaining bounded directions, including the feasible
          // residual gradient. Do not repeatedly accept the first tiny ray.
          if (stagnationFallback && directionDetails && !event.changed && proposal.step < 1e-3
            && residualDecrease.afterSquaredNorm > .99 * residualDecrease.beforeSquaredNorm) {
            if (residualDecrease.afterSquaredNorm < stagnationFallback.accepted.residualDecrease.afterSquaredNorm)
              stagnationFallback = { accepted, proposal };
            if (nextDirection('stalled-alternative-ray')) { accepted = undefined; trial = -1; continue; }
          }
          if (stagnationFallback && residualDecrease.afterSquaredNorm >= stagnationFallback.accepted.residualDecrease.afterSquaredNorm)
            ({ accepted, proposal } = stagnationFallback);
          break;
        } catch (error) {
          // A later grid/material gate may reject a successfully prepared
          // profile. Keep its local warnings and changes with that trial too.
          if (profile && error && typeof error === 'object' && !error.eventProfile)
            error.eventProfile = profile.diagnostics;
          const recoverableProposal = stage === 'coupled density proposal'
            && error.code === 'COUPLED_XFOIL_CANDIDATE_DOMAIN' && error.recoverable === true
            && Number.isFinite(error.step) && error.step > 0 && error.step <= maximumStep;
          if (error.diagnostics?.coordinateRepairBacktrack) wakeCoordinateRecovery = true;
          const step = proposal?.step ?? (recoverableProposal ? error.step : undefined);
          rejections.push({ step, stage, message: error.message,
            ...(error.code === undefined ? {} : { code: error.code }),
            ...(error.diagnostics === undefined ? {} : { diagnostics: structuredClone(error.diagnostics) }),
            ...(error.transitionRootFailure ? { transitionRootFailure: error.transitionRootFailure } : {}),
            ...(error.eventProfile ? { eventProfile: error.eventProfile } : {}),
            ...(error.coupledAdmissibility ? { admissibility: error.coupledAdmissibility } : {}) });
          if (trial >= maxBacktracks && error.code === 'COUPLED_RESIDUAL_DECREASE' && automaticFallback) {
            ({ accepted, proposal } = automaticFallback);
            recoveryOutcome = { kind: 'ordinary-policy-restored', reason: 'Automatic merit search exhausted',
              rejectedTrials: rejections.length, acceptedStep: proposal.step };
            recovery = null; plateauRecoveryDisabled = true; seedProgress();
            break;
          }
          if (trial >= maxBacktracks && stage === 'residual decrease' && transitionFallback) {
            ({ accepted, proposal } = transitionFallback);
            break;
          }
          if (trial >= maxBacktracks && error.code === 'COUPLED_RESIDUAL_DECREASE') {
            // Diagnose the failed Newton/projection pair before a gradient
            // fallback can mask it with indefinitely slow accepted steps.
            projectionDescentFailure = diagnoseProjectionDescent();
            if (projectionDescentFailure) throw error;
          }
          if (directionRecovery && trial >= maxBacktracks && error.code === 'COUPLED_RESIDUAL_DECREASE') {
            if (nextDirection('rejected-newton-ray')) { trial = -1; continue; }
          }
          if (stagnationFallback && trial >= maxBacktracks) {
            ({ accepted, proposal } = stagnationFallback);
            break;
          }
          if (stepAcceptance === 'listing' || trial >= maxBacktracks || stage === 'coupled density proposal' && !recoverableProposal) throw error;
          if (stage === 'Newton grid step' && error.code === 'streamtube-grid-step') geometryRedistribution = true;
          // A full trial-only SMOVE repair can spoil descent, yet removing
          // it entirely leaves the next direction pinned against the same
          // wake corner. Try one half correction at this valid Newton step.
          // If that also fails, shorten Newton and discard the stale repair
          // request. Do not exhaust the search damping coordinates when the
          // Newton step itself needs shortening. A genuine stagnation-motion
          // trigger retains the ordinary full-correction merit search.
          if (stage === 'residual decrease') {
            if (repairOnly && redistributionBacktracks === 0) redistributionBacktracks++;
            else { geometryRedistribution = false; repairOnly = false; }
          }
          if (redistributionAttempted && (stage === 'SMOVE' || stage === 'admissibility')) redistributionBacktracks++;
          maximumStep = stage === 'residual decrease' && repairOnly ? step
            : (error.code === 'streamtube-grid-step' ? error.stepFraction : .5) * step;
        }
      }
    } catch (error) {
      // A rejected station factor followed by a failed legacy factor still
      // establishes the fallback policy for this retained grid/checkpoint.
      if (linearOrdering === 'station-auto' && error.diagnostics?.stationPolicy?.recommendation === 'auto')
        stationFallback = true;
      const { x, ...details } = proposal ?? {};
      lastRejectedStep = { ...details, stage, message: error.message, rejections,
        ...(error.code === undefined ? {} : { code: error.code }),
        ...(error.diagnostics === undefined ? {} : { diagnostics: structuredClone(error.diagnostics) }),
        ...(error.transitionRootFailure ? { transitionRootFailure: error.transitionRootFailure } : {}),
        ...(error.eventProfile ? { eventProfile: error.eventProfile } : {}),
        ...(error.coupledAdmissibility ? { admissibility: error.coupledAdmissibility } : {}),
        ...(error.attempts ? { linearAttempts: error.attempts } : {}) };
      reason = `Coupled ISES update rejected during ${stage}: ${error.message}`; failed = true;
    } finally {
      // Restore even on an accepted trial: it commits a different, complete
      // system below. A rejected/aborted update leaves this one restartable.
      events.restore(phase); system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne)));
    }
    if (!accepted && hkProjectionRecovery && coupledOptions.hkFloorLinearization !== 'native'
      && lastRejectedStep?.code === 'COUPLED_RESIDUAL_DECREASE' && transitionUnchanged
      && proposal?.projection?.displacementChanges.some(p => p.kind === 'surface' && p.correction > 0)) {
      // A correct fixed-branch Jacobian can lose descent after DSLIM changes
      // the surface thickness. Diagnose the unprojected trial only: never
      // accept it or bypass the original physical/grid/residual gates.
      projectionDescentFailure ??= diagnoseProjectionDescent();
      if (projectionDescentFailure) {
        const { rawValue, projectedValue } = projectionDescentFailure;
        // Native floor sensitivities provide a quasi-Newton direction for
        // the same residual. Change them once, at the retained accepted
        // state; do not remarch the BL or alter the requested conditions.
        // Reuse the accepted undisplaced chart exactly. Rebasing subtracts
        // and reapplies displacement, which can round a physical node even
        // though changing Hk sensitivities must leave every residual intact.
        const next = createCoupledStreamtubeBody(input, { ...optionsAtState(), hkFloorLinearization: 'native',
          geometryReplay: 'preserve-undisplaced',
          initialEuler: { x: state.slice(0, system.ne), nodes: value.outer.nodes, undisplacedNodes: value.outer.undisplacedNodes },
          initialBL: state.slice(system.ne) });
        const nextState = next.initial.slice(), nextValue = next.admissibleValue(nextState, { requireConvex });
        if (!nextValue || nextState.some((v, k) => v !== state[k])
          || nextValue.residual.some((v, k) => v !== value.residual[k])
          || JSON.stringify(nextValue.outer.nodes) !== JSON.stringify(value.outer.nodes))
          throw new Error('Hk projection recovery changed the retained state, grid or residual.');
        const diagnostics = { from: 'exact', to: 'native', reason: lastRejectedStep.code,
          rawMerit: residualMerit(rawValue.residual), projectedMerit: residualMerit(projectedValue.residual),
          beforeMerit, step: proposal.step,
          rejectedTrials: rejections.length, equationsChanged: false };
        coupledOptions.hkFloorLinearization = 'native';
        coupledOptions.geometryReplay = 'preserve-undisplaced';
        system = next; state = nextState; value = nextValue;
        failed = false; reason = 'iteration limit'; lastRejectedStep = null;
        recovery = null;
        if (iterationRecovery) seedProgress();
        // This consumes a linearization, not an accepted Newton update.
        // Publish the same state with its new method and restart controls.
        report(iteration, { step: 0, accepted: false, hkProjectionRecovery: diagnostics });
        continue;
      }
    }
    if (!accepted) break;
    const { candidate, event, inlet, dekink, triggeredBodies, r, redistributionAttempted,
      geometryRedistribution, response, profile, residualDecrease } = accepted;
    ({ system, state, value } = candidate);
    if (beforeChange) {
      const stable = !event.changed && beforeChange.phase === JSON.stringify(system.bl.snapshotActive())
        && beforeChange.regimes.every((r, i) => r === system.bl.stations[i].regime)
        && beforeChange.mcrit === system.euler.conditions.upwind?.mcrit
        && (!dissipationEnhancement || input.upwind.mcrit === targetMcrit);
      convergenceDecision = coupledChangeDecision({ accepted: stable ? acceptedCoupledChanges(system, beforeChange, state, value) : undefined,
        newton: newtonChanges, residual: coupledResidualChanges(system, value.residual), families: value.families, stable });
    }
    // A requested-law root no longer needs the recovery policy. Subsequent
    // operating points start with ordinary backtracking again.
    if (maximum(value.residual) <= tolerance && (!dissipationEnhancement || input.upwind.mcrit === targetMcrit))
      wakeCoordinateRecovery = false;
    if (redistributionAttempted) lastRedistributedStagnation = value.outer.stagnation.slice();
    const { x, ...step } = proposal;
    if (accepted.transitionEventAcceptance)
      step.transitionEventAcceptance = { method: 'admissible-transition-event', residualComparison: 'different-equation-topology' };
    if (dissipationEnhancement) {
      const damping = auditDissipationDamping(targetMcrit, proposal.undampedUpdate.relativeDensity.maximum, proposal.step);
      previousDensityChange = damping.acceptedDensityChange;
      iterationDissipation.damping = damping;
    }
    let progressDetails;
    if (iterationRecovery) {
      progressDetails = progress.observe({ vector: physicalIterationVector(state, value.outer.nodes),
        before: beforeMerit, after: residualMerit(value.residual), residual: maximum(value.residual), tolerance,
        phase: phaseKey(), phaseChanged: event.changed, mcrit: input.upwind?.mcrit,
        maintenance: redistributionAttempted || dekink.repairs.length > 0 });
      if (recovery) {
        if (event.changed) recovery.comparable = false;
        recovery.remaining--;
        if (recovery.remaining === 0) {
          if (recovery.comparable && residualMerit(value.residual) >= .9 * recovery.merit && maximum(value.residual) > tolerance) {
            plateauRecoveryDisabled = true;
            recoveryOutcome = { kind: 'ordinary-policy-restored', reason: 'Damping did not accelerate residual reduction' };
          }
          recovery = null; seedProgress();
        }
      }
      if (progressDetails.cause && maximum(value.residual) > tolerance) {
        if (progressDetails.cause === 'negligible-state-change' || progressDetails.cause === 'two-state-cycle' && recoveries >= 2)
          progressStop = `Coupled progress stalled: ${progressDetails.cause}.`;
        else if (!recovery && !plateauRecoveryDisabled && recoveries < 2) {
          recovery = { remaining: 3, comparable: true, merit: residualMerit(value.residual), mcrit: input.upwind?.mcrit };
          recoveries++; seedProgress();
        }
      }
      progressDetails = { ...progressDetails, action: progressStop ? 'stop-stage' : recovery ? 'damped-recovery' : 'ordinary-newton',
        recoveries, ...(recoveryOutcome ? { recoveryOutcome } : {}), ...(recovery ? { frozenMcrit: recovery.mcrit, remaining: recovery.remaining } : {}),
        ...(iteration > maxIterations ? { extendedBudget: true } : {}) };
    }
    report(iteration, { ...step, ...(iterationDissipation ? { dissipation: { ...iterationDissipation, acceptedDensityChange: previousDensityChange } } : {}), backtracks: rejections.length, rejections,
      ...(convergenceDecision ? { convergence: structuredClone(convergenceDecision) } : {}),
      ...(progressDetails ? { progress: progressDetails } : {}),
      ...(residualDecrease ? { residualDecrease } : {}),
      ...(profile ? { eventProfile: profile } : {}),
      ...(event.changed ? { activeChange: true, changes: event.changes, meritComparable: false } : {}),
      maintenance: { inlet: { maxDisplacement: inlet.maxDisplacement, maxArcErrorBefore: inlet.maxArcErrorBefore, maxArcErrorAfter: inlet.maxArcErrorAfter,
          ...(inlet.reparameterization ? { reparameterization: inlet.reparameterization } : {}),
          ...(inlet.convexity ? { convexity: inlet.convexity } : {}) },
        dekinkRepairs: dekink.repairs, ...(dekink.convexity ? { dekinkConvexity: dekink.convexity } : {}),
        triggeredBodies, ...(geometryRedistribution ? { geometryRedistribution: true } : {}), passages: r.passages,
        ...(r.correctionScale === undefined ? {} : { correctionScale: r.correctionScale }),
        ...(response ? { projectionGeometry: response } : {}) } });
    if (convergenceDecision?.converged) { reason = 'solution changes'; break; }
    const stalledStep = directionRecovery && !dissipationEnhancement ? coupledStepStagnation(history, tolerance) : null;
    if (stalledStep) {
      reason = 'Coupled updates stagnated at roundoff-sized steps after guarded direction recovery.';
      lastRejectedStep = { stage: 'coupled progress', code: 'COUPLED_STEP_STAGNATION', message: reason, diagnostics: stalledStep };
      break;
    }
    if (progressStop) { reason = progressStop; break; }
  }
  let finalDissipationFailure;
  if (dissipationEnhancement) {
    try { changeMcrit(targetMcrit); }
    catch (error) { finalDissipationFailure = error.message; reason = `Final prescribed dissipation check failed: ${error.message}`; }
  }
  const result = coupledStreamtubeResult(system, state, { tolerance, reason, history, linearDiagnostics, initialRedistribution, lastRejectedStep,
    ...(convergenceDecision ? { convergence: convergenceDecision } : {}),
    ...(iterationRecovery ? { progressControl: { recoveries, stopReason: progressStop ?? null, monitorRestarted: !!resume } } : {}),
    ...(initialRedistribution.accepted ? { checkpoint: checkpoint() } : {}),
    solverInput: input, coupledOptions: optionsAtState(), iterationGeometry, stepAcceptance, stagnationLimiter,
    ...(preserveConvexity ? { gridAcceptance: 'convex' } : {}),
    ...(blUpdate === 'xfoil' ? { blUpdate } : {}),
    ...(projectionGeometry === 'boundary-increment' ? { projectionGeometry } : {}),
    ...(linearOrdering !== 'auto' ? { linearOrdering } : {}),
    ...(eventProfile !== 'none' ? { eventProfile } : {}),
      ...(shearCoordinate !== 'linear' ? { shearCoordinate } : {}),
    stepMethod: 'ises-density-newton', transitionEvents: true });
  if (!dissipationEnhancement) return result;
  const finalEquations = !finalDissipationFailure && input.upwind.mcrit === targetMcrit;
  if (!finalEquations) result.mesh.initialization.flowSolved = false;
  return { ...result, converged: result.converged && finalEquations,
    ...(!finalEquations ? { status: 'unconverged' } : {}),
    dissipationEnhancement: { method: 'mses', targetMcrit, previousDensityChange,
      finalMcrit: input.upwind.mcrit, finalEquations, ...(finalDissipationFailure ? { failure: finalDissipationFailure } : {}) } };

}
