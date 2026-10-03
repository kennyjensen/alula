// SPDX-License-Identifier: GPL-2.0-or-later
// One complete-state nested initial guess from an accepted coupled state.
// The final grid replaces the independently initialized requested grid;
// every Euler/grid/BL/transition/wake equation still requires a new solve.
import { coupledCheckpointShearCoordinate } from './streamtube-coupled-shear-policy.js';
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { coupledConvergenceSatisfied } from './streamtube-coupled-convergence.js';
import { refineCoupledStreamtubeBody } from './streamtube-coupled-refinement.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { initialStreamtubeDisplacement } from './streamtube-geometry.js';
import { prepareRefinedCoupledGridProfile } from './streamtube-coupled-grid-profile.js';

export const COUPLED_GRID_SEQUENCE_MAX_NODES = 50000;
const require = (ok, message) => { if (!ok) throw new Error(message); };
const same = (a, b) => {
  if (a === b) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  if (Array.isArray(a) || ArrayBuffer.isView(a)) return (Array.isArray(b) || ArrayBuffer.isView(b))
    && a.length === b.length && Array.from(a).every((v, i) => same(v, b[i]));
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};
const finiteVector = v => (Array.isArray(v) || ArrayBuffer.isView(v)) && v.length > 0
  && Array.from(v).every(Number.isFinite);
const finiteGrid = nodes => Array.isArray(nodes) && nodes.length > 0 && nodes.every(grid =>
  Array.isArray(grid) && grid.length >= 2 && grid.every(row => Array.isArray(row) && row.length >= 2
    && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))));
const finiteFamilies = f => ['euler', 'boundaryLayer', 'edgeMatching'].every(k => Number.isFinite(f?.[k]) && f[k] >= 0);
// An opt-in partial seed must come from ordinary accepted iterations, not
// a failed update, cancellation, stale chart, or fabricated convergence flag.
// Maintenance is checked for provenance but is not inherited by the child.
const acceptedPartialSource = result => {
  const cp = result?.checkpoint, saved = cp?.restart, c = cp?.continuation;
  const last = result?.history?.at(-1), bodies = saved?.input?.bodies;
  return result?.converged === false && result.reason === 'iteration limit' && result.status === 'unconverged'
    && result.lastRejectedStep === null && result.initialRedistribution?.accepted === true
    && cp?.version === 1 && saved?.options?.transitionMode === 'automatic'
    && finiteVector(result.x) && finiteVector(result.residual) && finiteFamilies(result.families)
    && same(cp.families, result.families)
    && Number.isInteger(last?.iteration) && last.iteration > 0
    && ['euler', 'boundaryLayer', 'edgeMatching'].every(k => last[k] === result.families[k])
    && Number.isFinite(last.residual) && last.residual === result.residual.reduce((m, v) => Math.max(m, Math.abs(v)), 0)
    && finiteVector(saved?.initialEuler?.x) && finiteVector(saved?.initialBL)
    && finiteGrid(saved?.initialEuler?.nodes) && finiteGrid(saved?.initialEuler?.undisplacedNodes)
    && Array.isArray(bodies) && bodies.length > 0
    && Array.isArray(saved.options.transitionState) && saved.options.transitionState.length === 2 * bodies.length
    && saved.options.transitionState.every(n => Number.isInteger(n) && n >= 0)
    && same(saved.options.transitionState, result.boundaryLayer?.transitionState)
    && ['convex', 'ises-sampled'].includes(c?.iterationGeometry) && ['listing', 'admissible', 'armijo', 'event-armijo'].includes(c?.stepAcceptance)
    && ['listing', 'prose'].includes(c?.stagnationLimiter) && ['amd', 'colamd'].includes(c?.preferredOrdering)
    && [.001, 1].includes(c?.pivotTolerance) && ['giles', 'xfoil'].includes(c?.blUpdate ?? 'giles')
    && ['auto', 'station', 'station-auto', 'aligned-auto'].includes(c?.linearOrdering ?? 'auto')
    && (c?.stationFallback === undefined || typeof c.stationFallback === 'boolean')
    && (c?.stationFallback !== true || c.linearOrdering === 'station-auto')
    && ['fixed', 'boundary-increment'].includes(c?.projectionGeometry ?? 'fixed')
    && (c?.projectionGeometry !== 'boundary-increment' || c?.blUpdate === 'xfoil')
    && finiteVector(c?.lastRedistributedStagnation) && c.lastRedistributedStagnation.length === bodies.length
    && Array.isArray(c?.fractions) && c.fractions.length === bodies.length && c.fractions.every((row, b) =>
      Array.isArray(row) && row.length === bodies[b].leadingIndex + 1 && row[0] === 0 && row.at(-1) === 1
      && row.every((v, i) => Number.isFinite(v) && (i === 0 || v > row[i - 1])));
};
const tubeRows = input => Array.isArray(input?.weights) && input.weights.length === input.bodies?.length + 1
  && input.weights.every(row => Array.isArray(row) && row.length >= 2 && row.every(w => Number.isFinite(w) && w > 0));

export function planCoupledGridSequence({ sourceInput, input, requestedGridIntervals, sourceGridIntervals } = {},
  { maxNodes = COUPLED_GRID_SEQUENCE_MAX_NODES } = {}) {
  require(Number.isInteger(maxNodes) && maxNodes > 0 && maxNodes <= COUPLED_GRID_SEQUENCE_MAX_NODES,
    'Invalid coupled grid-sequence node budget.');
  require([requestedGridIntervals, sourceGridIntervals].every(n => Number.isInteger(n) && n >= 8 && n <= 128),
    'Grid sequencing requires explicit valid source and requested surface interval counts.');
  require(tubeRows(sourceInput) && tubeRows(input) && sourceInput.bodies.length === input.bodies.length
    && Array.isArray(sourceInput.outerLower) && sourceInput.outerLower.length >= 3
    && Array.isArray(input.outerLower) && input.outerLower.length >= 3,
  'Grid sequencing requires complete source and requested passage dimensions.');
  const streamwiseFactor = Math.max(1, Math.ceil(requestedGridIntervals / sourceGridIntervals));
  require(streamwiseFactor <= 4, 'Requested grid sequencing needs more than four subdivisions per interval; a converged intermediate level is required.');
  const normalSubdivisions = sourceInput.weights.map((row, g) => {
    const target = Math.max(row.length, input.weights[g].length), counts = row.map(() => 1);
    require(target <= 4 * row.length, 'Requested grid sequencing needs more than four subdivisions per parent tube; a converged intermediate level is required.');
    // Each count partitions the SAME captured parent mass. Split its
    // largest remaining child mass, with stable upstream-index tie breaking.
    for (let total = row.length; total < target; total++) {
      let best = -1;
      for (let j = 0; j < row.length; j++) if (counts[j] < 4
        && (best < 0 || row[j] / counts[j] > row[best] / counts[best])) best = j;
      require(best >= 0, 'Grid sequencing exhausted its bounded normal subdivisions.');
      counts[best]++;
    }
    return counts;
  });
  require(streamwiseFactor > 1 || normalSubdivisions.some(row => row.some(n => n > 1)), 'Requested grid sequencing does not refine its parent.');
  const nx = sourceInput.outerLower.length - 1, tubes = normalSubdivisions.map(row => row.reduce((a, b) => a + b, 0));
  const nodeCount = (streamwiseFactor * nx + 1) * tubes.reduce((sum, n) => sum + n + 1, 0);
  require(Number.isSafeInteger(nodeCount) && nodeCount <= maxNodes,
    `Coupled grid sequencing needs ${nodeCount} nodes, exceeding its ${maxNodes} node budget.`);
  return { streamwiseFactor, normalSubdivisions, normalInterpolation: 'streamfunction-quadratic',
    streamwiseInterpolation: 'surface-pchip', maxNodes, nodeCount,
    sourceGridIntervals, requestedGridIntervals, nominalGridIntervals: streamwiseFactor * sourceGridIntervals,
    parent: { nx, tubes: sourceInput.weights.map(row => row.length) },
    requested: { nx: input.outerLower?.length - 1, tubes: input.weights.map(row => row.length) },
    child: { nx: streamwiseFactor * nx, tubes }, requestedGridRetained: false };
}

export function prepareCoupledGridSequence({ sourceSystem, sourceResult, input, options,
  requestedGridIntervals, sourceGridIntervals } = {}, { tolerance = 1e-10, maxNodes = COUPLED_GRID_SEQUENCE_MAX_NODES,
  allowPartialSource = false, blPredictor = 'interpolate' } = {}) {
  require(Number.isFinite(tolerance) && tolerance > 0, 'Invalid coupled grid-sequence tolerance.');
  require(['interpolate', 'xfoil-mrchdu'].includes(blPredictor), 'Invalid coupled grid-sequence BL predictor.');
  require(typeof allowPartialSource === 'boolean', 'Invalid partial coupled grid-sequence policy.');
  const partialSource = allowPartialSource && acceptedPartialSource(sourceResult);
  const sourceConverged = sourceResult?.converged === true && coupledConvergenceSatisfied(sourceResult, tolerance);
  require((sourceConverged || partialSource) && sourceResult.mesh?.quality?.valid === true && sourceSystem?.bl && sourceSystem.euler,
    allowPartialSource ? 'Grid sequencing requires a converged or complete accepted iteration-limit coupled source.'
      : 'Grid sequencing requires a converged accepted coupled source.');
  const shearCoordinate = coupledCheckpointShearCoordinate(sourceResult.checkpoint);
  const sourceInput = sourceResult.checkpoint?.restart?.input ?? sourceResult.solverInput;
  const plan = planCoupledGridSequence({ sourceInput, input, requestedGridIntervals, sourceGridIntervals }, { maxNodes });
  require(input.displacement === undefined && sourceInput.displacement === undefined,
    'Coupled grid sequencing owns its BL displacement.');
  require(sourceSystem.bl.transitionMode === 'automatic' && options?.transitionMode === 'automatic'
    && sourceSystem.bl.trips.every(pair => pair.every(t => t === 1)),
  'This grid-sequence startup requires automatic transition with terminal material trips.');
  const allowedGridKeys = new Set(['weights', 'gridSpacing', 'bodies', 'outerLower', 'outerUpper', 'cutPaths']);
  for (const key of new Set([...Object.keys(sourceInput), ...Object.keys(input)])) if (!allowedGridKeys.has(key))
    require(same(sourceInput[key], input[key]), `Grid sequencing changed physical/model input ${key}.`);
  const elements = new Set();
  sourceInput.bodies.forEach((body, b) => {
    const target = input.bodies[b];
    require(Number.isInteger(body.element) && body.element >= 0 && !elements.has(body.element)
      && target.element === body.element, 'Grid sequencing changed explicit body identities or passage order.');
    elements.add(body.element);
    require(same(body.points, target.points) && same(body.trailingEdge, target.trailingEdge),
      'Grid sequencing changed a physical solid contour or trailing-edge topology.');
    require(same(body.points, sourceSystem.euler.layout.bodies[b]?.points)
      && same(body.trailingEdge, sourceSystem.euler.layout.bodies[b]?.trailingEdge)
      && body.element === sourceSystem.euler.layout.bodies[b]?.element,
    'Coupled source system and checkpoint have different material bodies.');
  });
  // Compare effective target conditions, including normalization and gas
  // equations, without evaluating its flow or reinitializing any BL profile.
  const target = createStreamtubeBodySystem({ ...input, displacement: initialStreamtubeDisplacement({
    bodies: input.bodies, nx: input.outerLower.length - 1 }, sourceSystem.euler.baseGeometry) });
  require(same(target.conditions, sourceSystem.euler.conditions), 'Grid sequencing changed effective Euler conditions or normalization.');
  const sourceOptions = sourceResult.checkpoint?.restart?.options ?? sourceResult.coupledOptions;
  require(sourceOptions, 'Grid sequencing requires the accepted source physical options.');
  for (const key of ['reynolds', 'ncrit', 'edgeMatching', 'transitionMode', 'blThermodynamics'])
    require(same(options[key], sourceOptions[key]), `Grid sequencing changed coupled physical/model option ${key}.`);
  require((options.hkFloorLinearization ?? 'exact') === (sourceOptions.hkFloorLinearization ?? 'exact')
    && (sourceSystem.conditions.hkFloorLinearization ?? 'exact') === (sourceOptions.hkFloorLinearization ?? 'exact'),
  'Grid sequencing changed the explicit Hk-floor linearization policy.');
  require((options.geometryReplay ?? 'legacy') === (sourceOptions.geometryReplay ?? 'legacy')
    && (sourceSystem.conditions.geometryReplay ?? 'legacy') === (sourceOptions.geometryReplay ?? 'legacy'),
  'Grid sequencing changed the explicit geometry replay policy.');
  require(same(options.tripFractions ?? input.bodies.map(() => [1, 1]), sourceSystem.bl.trips),
    'Grid sequencing moved its material trips.');
  require(sourceSystem.euler.layout.independentWakeBanks && target.layout.independentWakeBanks,
    'This grid-sequence startup requires independent wake banks.');
  // Reconstruct the saved source independently so an unrelated source
  // object cannot supply a valid-looking convergence envelope.
  const saved = sourceResult.checkpoint?.restart;
  require(saved?.initialEuler && saved.initialBL && saved.options, 'Grid sequencing requires a complete accepted source checkpoint.');
  const replay = createCoupledStreamtubeBody(sourceInput, { ...saved.options,
    initialEuler: saved.initialEuler, initialBL: saved.initialBL });
  const before = replay.evaluate(replay.initial);
  require(same(before.residual, sourceResult.residual) && same(before.families, sourceResult.families)
    && before.residual.every(Number.isFinite)
    && (partialSource || coupledConvergenceSatisfied({ ...sourceResult, families: before.families }, tolerance)) && replay.admissible(replay.initial),
  'Accepted coupled source does not replay exactly.');
  if (partialSource) require(same(replay.initial, sourceResult.x)
    && same(replay.bl.snapshotActive(), saved.options.transitionState)
    && same(sourceSystem.bl.snapshotActive(), saved.options.transitionState)
    && same(before.outer.nodes, saved.initialEuler.nodes),
  'Accepted partial coupled source state, transition map or physical grid does not replay exactly.');
  require(same(replay.initial, sourceSystem.initial) && same(replay.euler.decode(replay.initial.slice(0, replay.ne)).nodes,
    sourceSystem.euler.decode(sourceSystem.initial.slice(0, sourceSystem.ne)).nodes),
  'Grid sequence source system does not match the accepted source checkpoint.');
  const sourceMesh = streamtubeMeshSnapshot({ system: replay.euler, nodes: before.outer.nodes });
  require(sourceMesh.quality.valid, 'Replayed grid-sequence source mesh is invalid.');
  const refined = refineCoupledStreamtubeBody(sourceInput, replay, plan), system = refined.system;
  const value = system.evaluate(system.initial), mesh = streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes });
  require(mesh.quality.valid && system.admissible(system.initial) && value.residual.every(Number.isFinite),
    'Refined complete-state seed failed its physical/grid domain.');
  require(same(system.conditions, replay.conditions) && same(system.euler.conditions, replay.euler.conditions)
    && same(system.bl.trips, replay.bl.trips), 'Complete-state refinement changed effective physics or material trips.');
  const surfaceIntervals = system.euler.layout.bodies.map(b => ({ element: b.element, intervals: b.trailingIndex - b.leadingIndex }));
  // The refiner's topology is authoritative. Keep element identities rather
  // than relabeling sorted passages with their array indices.
  refined.input.gridSpacing.surfaceIntervalsByElement = structuredClone(surfaceIntervals);
  const transfer = { method: partialSource ? 'partial-coupled-grid-sequence' : 'converged-coupled-grid-sequence',
    ...(partialSource ? { sourceConverged: false, sourceReason: sourceResult.reason,
      sourceIteration: sourceResult.history.at(-1).iteration } : {}), initialGuessOnly: true, equationsChanged: false,
    requestedGridRetained: false, nominalGridIntervals: plan.nominalGridIntervals,
    requestedGridIntervals, sourceGridIntervals, parent: plan.parent, requested: plan.requested,
    child: { ...plan.child, cells: mesh.cells.length, nodes: plan.nodeCount, surfaceIntervals },
    streamwiseFactor: plan.streamwiseFactor, normalSubdivisions: plan.normalSubdivisions,
    normalInterpolation: plan.normalInterpolation, streamwiseInterpolation: plan.streamwiseInterpolation,
    sourceFamilies: { ...before.families }, targetFamilies: { ...value.families },
    sourceTransitionState: replay.bl.snapshotActive(), targetTransitionState: system.bl.snapshotActive(),
    refinement: refined.diagnostics, thicknessFactor: 1, maintenanceHistoryInherited: false,
    operations: { sourceResidualReplays: 1, nativeBLInitializations: 0, globalNewtonUpdates: 0, globalLinearSolves: 0 } };
  mesh.initialization.gridRefinement = structuredClone(transfer);
  const prepared = { system, input: refined.input, options: { ...refined.options, transitionState: system.bl.snapshotActive() },
    // A new grid resets its maintenance chart, not its source update policy.
    iterationControls: { shearCoordinate },
    initialEuler: refined.initialEuler, initialBL: system.initial.slice(system.ne), mesh, value, transfer,
    initialization: { method: partialSource ? 'Complete-state nested grid sequencing from an unconverged accepted coupled state'
      : 'Complete-state nested grid sequencing from a converged coupled root', thicknessFactor: 1,
      history: [{ attempt: 0, accepted: true, thicknessFactor: 1, quality: mesh.quality, families: value.families }],
      equationsChanged: false, flowSolved: false, transfer: structuredClone(transfer) } };
  return blPredictor === 'xfoil-mrchdu' ? prepareRefinedCoupledGridProfile(prepared) : prepared;
}
