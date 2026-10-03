// SPDX-License-Identifier: GPL-2.0-or-later
import { observableFlow, observableBL } from './streamtube-flow-preview.js';
// Bounded startup recovery for repeatedly cycling natural transition.
// Refinement constructs a new guess. All Euler/BL/wake equations, physical
// controls and final acceptance tests still have to converge on that grid.
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from './streamtube-coupled-refinement.js';
import { solveCoupledStreamtubeIses } from './streamtube-coupled-ises.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { prepareRefinedCoupledGridProfile } from './streamtube-coupled-grid-profile.js';

const familyMaximum = families => families && ['euler', 'boundaryLayer', 'edgeMatching']
  .every(k => Number.isFinite(families[k]) && families[k] >= 0)
  ? Math.max(families.euler, families.boundaryLayer, families.edgeMatching) : Infinity;
const finiteVector = row => (Array.isArray(row) || ArrayBuffer.isView(row))
  && row.length > 0 && row.every(Number.isFinite);
function completeCheckpoint(checkpoint) {
  const f = checkpoint?.restart;
  return checkpoint?.version === 1 && Number.isFinite(familyMaximum(checkpoint.families))
    && f?.input && f.options?.transitionMode === 'automatic' && checkpoint.continuation
    && finiteVector(f.initialEuler?.x) && finiteVector(f.initialBL)
    && Array.isArray(f.initialEuler?.nodes) && Array.isArray(f.initialEuler?.undisplacedNodes)
    && Array.isArray(f.options.transitionState) && f.options.transitionState.every(Number.isInteger);
}

// One resolution retry after a constrained solve has exhausted its budget.
// Locate the actual dominant BL row; a geometry or farfield failure alone
// cannot select surface refinement. This is a work-selection heuristic,
// not a relaxed convergence or physical-accuracy criterion.
export function stalledBoundaryLayerResolutionPlan(result, { maxNodes = 50000 } = {}) {
  if (result?.converged || result?.automaticRefinement || !result?.mesh?.quality?.valid
    || !completeCheckpoint(result.checkpoint) || !Number.isInteger(maxNodes) || maxNodes < 1
    || maxNodes > 50000 || !['armijo', 'event-armijo'].includes(result.checkpoint.continuation.stepAcceptance)) return null;
  const f = result.families;
  if (!(f.boundaryLayer > 2 * Math.max(f.euler, f.edgeMatching))) return null;
  const recent = result.history?.filter(h => h.iteration > 0).slice(-8) ?? [];
  if (recent.length < 4 || recent.filter(h => h.step > 0 && h.step < .02).length < 3) return null;
  const cp = result.checkpoint, ne = cp.restart.initialEuler.x.length, residual = result.residual;
  if (!residual || residual.length !== ne + cp.restart.initialBL.length) return null;
  let row = ne;
  for (let i = ne + 1; i < residual.length; i++) if (Math.abs(residual[i]) > Math.abs(residual[row])) row = i;
  const id = Math.floor((row - ne) / 4), bl = result.boundaryLayer;
  const troubled = bl.surfaces.find(s => s.ids.includes(id));
  if (!troubled) return null;
  const natural = s => bl.transitions.some(t => t.body === s.body && t.side === s.side
    && t.kind === 'natural' && t.forced === false);
  const nearTransition = Math.abs(troubled.ids.indexOf(id) - troubled.transition) <= 2 && natural(troubled);
  // Preserve the transition-window recovery when onset is unresolved. Away
  // from onset, refine where the equations actually stall: a laminar or
  // turbulent interval can be under-resolved without a transition event.
  const selected = nearTransition ? bl.surfaces.filter(natural) : [troubled];
  const input = cp.restart.input, nx = input.outerLower.length - 1, counts = Array(nx).fill(1), surfaces = [];
  for (const s of selected) {
    const station = bl.stations[nearTransition ? s.ids[s.transition] : id], body = input.bodies[s.body], intervals = [];
    for (let i = Math.max(body.leadingIndex, station.i - 5); i < Math.min(body.trailingIndex, station.i + (nearTransition ? 3 : 4)); i++) {
      counts[i] = 2; intervals.push(i);
    }
    surfaces.push({ body: s.body, side: s.side, dominantStation: id,
      ...(nearTransition ? { transitionIndex: s.transition, transitionStation: station.i } : { station: station.i }), intervals });
  }
  const refinedNx = counts.reduce((sum, n) => sum + n, 0);
  const nodeCount = (refinedNx + 1) * input.weights.reduce((sum, row) => sum + row.length + 1, 0);
  if (refinedNx === nx || nodeCount > maxNodes) return null;
  return { reason: nearTransition ? 'stalled-natural-transition-resolution' : 'stalled-boundary-layer-resolution', source: 'accepted-complete-checkpoint',
    dominantStation: id, dominantResidual: residual[row], surfaces, parentNx: nx, refinedNx,
    streamwiseSubdivisions: counts, normalFactor: 1, maxNodes, nodeCount,
    selection: { kind: 'bounded-work-heuristic', maximumRefinements: 1, equationsChanged: false } };
}

// Keep one accepted complete state, never an independently interpolated BL or
// a rejected proposal. This is selection of bounded recovery work, not a new
// convergence criterion. The driver already enforces gas/BL/convexity gates.
export function retainBestCoupledCheckpoint(best, checkpoint, { iteration, initialRedistribution } = {}) {
  if (initialRedistribution?.accepted !== true || !Number.isInteger(iteration) || iteration < 0
    || familyMaximum(checkpoint?.families) >= familyMaximum(best?.checkpoint?.families)
    || !completeCheckpoint(checkpoint)) return best;
  return { iteration, checkpoint: structuredClone(checkpoint) };
}

function boundaryLayerStallPlan(result, bestState, tolerance, maxNodes) {
  if (!completeCheckpoint(bestState?.checkpoint) || !Number.isInteger(bestState.iteration)
    || bestState.iteration < 0 || ((result.solverStopReason ?? result.reason) !== 'iteration limit'
      && !(result.progressControl?.stopReason && result.progressControl.stopReason === result.reason))) return null;
  const best = bestState.checkpoint, final = result.checkpoint, history = result.history;
  const score = familyMaximum(best.families), endScore = familyMaximum(final.families), last = history.at(-1);
  if (!Number.isFinite(score) || !Number.isFinite(endScore) || score <= 10 * tolerance
    || best.families.boundaryLayer < 2 * Math.max(best.families.euler, best.families.edgeMatching)
    || final.families.boundaryLayer < 2 * Math.max(final.families.euler, final.families.edgeMatching)
    || !Number.isInteger(last?.iteration) || last.iteration - bestState.iteration < 8
    || endScore < 2 * score) return null;
  // A saved best state can only seed the same operating point and grid.
  const withoutPhase = options => ({ ...options, transitionState: undefined });
  if (JSON.stringify(best.restart.input) !== JSON.stringify(final.restart.input)
    || JSON.stringify(withoutPhase(best.restart.options)) !== JSON.stringify(withoutPhase(final.restart.options))) return null;
  const recent = history.filter(h => Number.isInteger(h.iteration) && h.iteration > 0).slice(-8);
  if (recent.length !== 8 || recent.some((h, i) => !Number.isFinite(h.step) || h.step <= 0 || h.step > 1
    || i && h.iteration !== recent[i - 1].iteration + 1)) return null;
  const input = best.restart.input, nx = input.outerLower.length - 1, subdivisions = Array(nx).fill(1), surfaces = [];
  for (const [index, surface] of result.boundaryLayer.surfaces.entries()) {
    const phase = best.restart.options.transitionState[index], ids = surface.ids;
    if (!Array.isArray(ids) || !Number.isInteger(phase) || phase < 1 || phase >= ids.length - 1) continue;
    const transition = result.boundaryLayer.transitions?.find(t => t.body === surface.body && t.side === surface.side);
    if (!transition || transition.forced !== false || transition.kind !== 'natural') continue;
    const evidence = recent.flatMap(h => {
      const limiter = h.viscousLimiter;
      if (limiter?.kind !== 'xfoil-bl-update' || !Number.isInteger(limiter.station)
        || !Number.isFinite(limiter.normalizedIncrement) || !Number.isFinite(limiter.bound)
        || limiter.normalizedIncrement === 0 || limiter.bound === 0
        || limiter.bound / limiter.normalizedIncrement <= 0
        || limiter.bound / limiter.normalizedIncrement >= 1) return [];
      const stationIndex = ids.indexOf(limiter.station);
      return stationIndex >= 0 && Math.abs(stationIndex - phase) <= 2
        ? [{ iteration: h.iteration, station: limiter.station, stationIndex, variable: limiter.variable,
          normalizedIncrement: limiter.normalizedIncrement, bound: limiter.bound }] : [];
    });
    if (evidence.length < 4) continue;
    const body = input.bodies[surface.body], intervals = new Set(), stations = [];
    for (let k = Math.max(0, phase - 6); k <= Math.min(ids.length - 1, phase + 1); k++) {
      const station = result.boundaryLayer.stations[ids[k]];
      if (!station || !Number.isInteger(station.i) || station.body !== surface.body || station.side !== surface.side) return null;
      stations.push(station.id);
      for (const i of [station.i - 1, station.i]) if (i >= body.leadingIndex && i < body.trailingIndex) {
        subdivisions[i] = 2; intervals.add(i);
      }
    }
    surfaces.push({ body: surface.body, side: surface.side, transitionIndex: phase,
      transitionStation: result.boundaryLayer.stations[ids[phase]].i,
      stations, intervals: [...intervals].sort((a, b) => a - b), evidence });
  }
  if (!surfaces.length) return null;
  const refinedNx = subdivisions.reduce((sum, n) => sum + n, 0);
  const nodeCount = (refinedNx + 1) * input.weights.reduce((sum, row) => sum + row.length + 1, 0);
  if (nodeCount > maxNodes) return null;
  return { reason: 'boundary-layer-transition-stall', source: 'best-accepted-checkpoint',
    bestIteration: bestState.iteration, terminalIteration: last.iteration,
    bestFamilies: { ...best.families }, terminalFamilies: { ...final.families },
    selection: { kind: 'bounded-work-heuristic', minimumStalledUpdates: 8, minimumLimitedUpdates: 4,
      maximumTransitionStationDistance: 2, upstreamStations: 6, downstreamStations: 1,
      acceptanceChanged: false },
    surfaces, streamwiseSubdivisions: subdivisions, normalFactor: 1, maxNodes, nodeCount, parentNx: nx, refinedNx };
}

export function transitionRecoveryPlan(result, { maxNodes = 50000, bestState, tolerance = 1e-10 } = {}) {
  if (!Number.isInteger(maxNodes) || maxNodes <= 0 || maxNodes > 50000) throw new Error('Invalid transition-recovery node budget.');
  if (!Number.isFinite(tolerance) || tolerance <= 0) throw new Error('Invalid transition-recovery tolerance.');
  if (result?.converged || result?.conditions?.transitionMode !== 'automatic'
    || result?.automaticRefinement || !result.mesh?.quality?.valid || result.checkpoint?.version !== 1) return null;
  const input = result.checkpoint.restart?.input, history = result.history;
  if (!history?.length || !input?.weights || !input.outerLower || !result.boundaryLayer?.surfaces) return null;
  const nx = input.outerLower.length - 1, subdivisions = Array(nx).fill(1), surfaces = [];
  for (const surface of result.boundaryLayer.surfaces) {
    const events = history.flatMap(h => (h.changes ?? []).filter(c => c.body === surface.body && c.side === surface.side
      && c.from !== c.to).map(c => ({ ...c, iteration: h.iteration })));
    // Two complete reversals, still active in the last four updates. These
    // are work-selection controls, not a modified transition criterion.
    const recent = events.slice(-4), last = recent.at(-1);
    if (recent.length < 4 || last.to !== surface.transition || last.iteration < history.at(-1).iteration - 4
      || recent.some((c, j) => c.kind !== 'natural' || Math.abs(c.to - c.from) !== 1
        || j && (c.from !== recent[j - 1].to || c.to !== recent[j - 1].from))) continue;
    const station = result.boundaryLayer.stations[surface.ids[surface.transition]];
    if (!station || !Number.isInteger(station.i)) continue;
    const body = input.bodies[surface.body], intervals = [];
    for (let i = Math.max(body.leadingIndex + 1, station.i - 3); i <= Math.min(body.trailingIndex, station.i + 1); i++) {
      subdivisions[i - 1] = 4; intervals.push(i);
    }
    surfaces.push({ body: surface.body, side: surface.side, intervals,
      transitionIndex: surface.transition, transitionStation: station.i,
      events: recent.map(({ from, to, iteration }) => ({ from, to, iteration })) });
  }
  if (!surfaces.length) return boundaryLayerStallPlan(result, bestState, tolerance, maxNodes);
  const refinedNx = subdivisions.reduce((sum, n) => sum + n, 0);
  const nodeCount = (refinedNx + 1) * input.weights.reduce((sum, row) => sum + row.length + 1, 0);
  if (nodeCount > maxNodes) return null;
  return { reason: 'repeated-adjacent-natural-transition', surfaces, streamwiseSubdivisions: subdivisions,
    normalFactor: 1, maxNodes, nodeCount, parentNx: nx, refinedNx };
}

export function recoverCoupledTransition(result, { plan = transitionRecoveryPlan(result), maxIterations = 40,
  tolerance = 1e-10, blPredictor = 'interpolate', onIteration, onMesh, onStage, onIterationCheckpoint, onFlow, normalization, startupAttempt } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || tolerance <= 0)
    throw new Error('Invalid transition-recovery iteration controls.');
  if (!['interpolate', 'xfoil-mrchdu'].includes(blPredictor)) throw new Error('Unknown transition-refinement BL predictor.');
  if (!plan) return null;
  const f = result.checkpoint.restart;
  onStage?.({ stage: 'transition-refinement', transitionRecovery: plan });
  const parent = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const parentValue = parent.evaluate(parent.initial);
  if (Object.keys(result.families).some(k => parentValue.families[k] !== result.families[k]))
    throw new Error('Transition-recovery parent residual does not replay exactly.');
  if (plan.reason === 'boundary-layer-transition-stall' && plan.surfaces.some(surface =>
    !result.boundaryLayer.transitions.some(t => t.body === surface.body && t.side === surface.side
      && t.kind === 'natural' && t.forced === false)))
    throw new Error('BL-stall recovery requires natural transition in the retained best state.');
  const mapped = refineCoupledStreamtubeBody(f.input, parent, plan);
  const smoothing = result.mesh.initialization?.gridSmoothing;
  const automaticRefinement = { kind: plan.reason === 'stalled-boundary-layer-resolution' ? 'boundary-layer-local' : 'transition-local', ...plan, parentUnknowns: parent.n,
    unknowns: mapped.system.n, parentFamilies: { ...result.families }, initialization: mapped.diagnostics,
    ...(blPredictor === 'xfoil-mrchdu' ? { parentNx: parent.euler.layout.nx, refinedNx: mapped.system.euler.layout.nx } : {}) };
  let preparedCheckpoint;
  if (blPredictor === 'xfoil-mrchdu') {
    if (!result.converged || familyMaximum(result.families) > tolerance)
      throw new Error('Native startup refinement requires a complete converged parent.');
    // Establish the new grid's initial SMOVE/maintenance chart once, then
    // prepare native BL profiles before any Newton update on that grid.
    const c = result.checkpoint.continuation;
    const initialized = solveCoupledStreamtubeIses(mapped.input, { ...mapped.options,
      initialEuler: mapped.initialEuler, initialBL: mapped.initialBL, maxIterations: 0, tolerance,
      iterationGeometry: c.iterationGeometry, stepAcceptance: c.stepAcceptance, stagnationLimiter: c.stagnationLimiter,
      blUpdate: c.blUpdate, projectionGeometry: c.projectionGeometry, linearOrdering: c.linearOrdering ?? 'auto',
      eliminateAmplification: c.eliminateAmplification ?? false,
      ...(c.shearCoordinate === undefined ? {} : { shearCoordinate: c.shearCoordinate }) });
    if (!initialized.checkpoint || initialized.initialRedistribution?.accepted !== true || !initialized.mesh?.quality?.valid)
      throw Object.assign(new Error(`Native startup refinement could not initialize its new grid: ${initialized.reason}`), {
        code: 'coupled-refinement-initialization', diagnostics: structuredClone({
          reason: initialized.reason, families: initialized.families,
          initialRedistribution: initialized.initialRedistribution, lastRejectedStep: initialized.lastRejectedStep,
          meshQuality: initialized.mesh?.quality, checkpointAvailable: !!initialized.checkpoint,
        }),
      });
    const prepared = prepareRefinedCoupledGridProfile(initialized.checkpoint.restart,
      { wakeWidthIncrement: true, reinitializeAmplification: true });
    preparedCheckpoint = { ...initialized.checkpoint, families: { ...prepared.value.families }, restart: {
      input: prepared.input, options: prepared.options, initialEuler: prepared.initialEuler, initialBL: prepared.initialBL } };
    automaticRefinement.profilePreparation = prepared.transfer.profilePreparation;
  }
  let checkpoint;
  const meshCallback = state => {
    const mesh = streamtubeMeshSnapshot(state);
    mesh.initialization.gridRefinement = automaticRefinement;
    if (smoothing) mesh.initialization.gridSmoothing = smoothing;
    onMesh?.(mesh, state.iteration.iteration === 0 ? 'initial' : 'solving', 'coupled');
    if (onFlow && checkpoint) onFlow(structuredClone({ checkpoint,
      flow: observableFlow(state.flow),
      bl: observableBL(mapped.system.bl),
      bodies: state.system.layout.bodies, normalization, iteration: state.iteration,
      mach: checkpoint.restart.input.mach, actualNcrit: checkpoint.restart.options.ncrit,
      stage: 'coupled', startupAttempt, automaticRefinement }));
  };
  onStage?.({ stage: 'coupled', transitionRecovery: automaticRefinement });
  const solved = solveCoupledStreamtubeIses(preparedCheckpoint ? undefined : mapped.input, {
    convergence: result.checkpoint.continuation.convergence,
    ...(preparedCheckpoint ? { resume: preparedCheckpoint } : { ...mapped.options,
      initialEuler: mapped.initialEuler, initialBL: mapped.initialBL }), maxIterations, tolerance,
    ...(preparedCheckpoint ? { iterationGeometry: preparedCheckpoint.continuation.iterationGeometry,
      stepAcceptance: preparedCheckpoint.continuation.stepAcceptance,
      stagnationLimiter: preparedCheckpoint.continuation.stagnationLimiter,
      blUpdate: preparedCheckpoint.continuation.blUpdate ?? 'giles' }
      : { iterationGeometry: 'ises-sampled',
        stepAcceptance: ['armijo', 'event-armijo'].includes(result.checkpoint.continuation.stepAcceptance)
          ? result.checkpoint.continuation.stepAcceptance : 'admissible', blUpdate: 'xfoil' }),
    // Keep the selected source ordering, including legacy omission as auto.
    // This new grid starts fresh factorization history, so a prior stage's
    // stationFallback marker is deliberately not carried over.
    linearOrdering: result.checkpoint.continuation?.linearOrdering ?? 'auto',
    directionRecovery: result.checkpoint.continuation?.directionRecovery ?? false,
    eliminateAmplification: result.checkpoint.continuation?.eliminateAmplification ?? false,
    ...(result.checkpoint.continuation?.shearCoordinate === undefined ? {}
      : { shearCoordinate: result.checkpoint.continuation.shearCoordinate }),
    // Refinement starts a new chart, but retains the source update policy.
    // Old checkpoints without this field keep the lower-level fixed default.
    ...(result.checkpoint.continuation?.projectionGeometry === undefined ? {}
          : { projectionGeometry: result.checkpoint.continuation.projectionGeometry }), onIteration,
    onCheckpoint: (value, details) => {
      checkpoint = value;
      onIterationCheckpoint?.(structuredClone(value), structuredClone({ ...details, stage: 'coupled',
        startupAttempt, automaticRefinement }));
    },
    onMesh: meshCallback });
  solved.mesh.initialization.gridRefinement = automaticRefinement;
  if (smoothing) solved.mesh.initialization.gridSmoothing = smoothing;
  return { ...solved, automaticRefinement };
}
