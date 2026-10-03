// SPDX-License-Identifier: GPL-2.0-or-later
// Bounded complete-state grid preparation. Only intermediate levels
// are solved here; the requested final grid returns an initial guess.
import { coupledCheckpointShearCoordinate, coupledFreshShearCoordinate } from './streamtube-coupled-shear-policy.js';
import { planCoupledGridSequence, prepareCoupledGridSequence, COUPLED_GRID_SEQUENCE_MAX_NODES } from './streamtube-coupled-grid-sequence.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { observableFlow, observableBL } from './streamtube-flow-preview.js';
import { solveCoupledStreamtubeIses } from './streamtube-coupled-ises.js';
import { coupledConvergenceSatisfied } from './streamtube-coupled-convergence.js';

const require = (ok, message) => { if (!ok) throw new Error(message); };
const copy = structuredClone;
const same = (a, b) => {
  if (a === b) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  if (Array.isArray(a) || ArrayBuffer.isView(a)) return (Array.isArray(b) || ArrayBuffer.isView(b))
    && a.length === b.length && Array.from(a).every((value, i) => same(value, b[i]));
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(key => Object.hasOwn(b, key) && same(a[key], b[key]));
};
const dimensions = input => ({ nx: input.outerLower.length - 1, tubes: input.weights.map(row => row.length) });
const finiteVector = value => (Array.isArray(value) || ArrayBuffer.isView(value))
  && value.length > 0 && value.every(Number.isFinite);

// The target is only a requested passage-count template for the existing
// nested preparer. Aggregation preserves its positive total weight without
// selecting or changing any captured source mass. Actual child weights come
// exclusively from the original planner's parent-mass subdivisions.
const aggregateWeights = (weights, count) => {
  if (weights.length <= count) return weights.slice();
  const reduced = Array(count).fill(0);
  weights.forEach((value, j) => { reduced[Math.floor(j * count / weights.length)] += value; });
  return reduced;
};

export function planCoupledGridLevels({ sourceInput, input, sourceGridIntervals, requestedGridIntervals } = {},
  { maxNodes = COUPLED_GRID_SEQUENCE_MAX_NODES, maxStreamwiseFactor = 2 } = {}) {
  require(Number.isInteger(maxStreamwiseFactor) && maxStreamwiseFactor >= 2 && maxStreamwiseFactor <= 4,
    'Invalid coupled grid-level subdivision bound.');
  require([sourceGridIntervals, requestedGridIntervals].every(n => Number.isInteger(n) && n >= 8 && n <= 128)
    && requestedGridIntervals >= sourceGridIntervals,
  'Grid levels require ordered source and requested surface interval counts.');
  require(Number.isInteger(maxNodes) && maxNodes > 0 && maxNodes <= COUPLED_GRID_SEQUENCE_MAX_NODES,
    'Invalid coupled grid-level node budget.');
  const valid = value => Array.isArray(value?.bodies) && value.bodies.length > 0
    && Array.isArray(value.outerLower) && value.outerLower.length >= 3
    && Array.isArray(value.weights) && value.weights.length === value.bodies.length + 1
    && value.weights.every(row => Array.isArray(row) && row.length >= 2 && row.every(w => Number.isFinite(w) && w > 0));
  require(valid(sourceInput) && valid(input) && sourceInput.bodies.length === input.bodies.length,
    'Grid levels require complete matching source and target dimensions.');
  const requested = dimensions(input), source = dimensions(sourceInput), levels = [];
  let nominal = sourceGridIntervals, virtualSource = sourceInput;
  for (let index = 0; index < 8; index++) {
    const stageIntervals = Math.min(requestedGridIntervals, nominal * maxStreamwiseFactor);
    const targetWeights = input.weights.map((row, g) => aggregateWeights(row, 4 * virtualSource.weights[g].length));
    const plan = planCoupledGridSequence({ sourceInput: virtualSource, input: { ...input, weights: targetWeights },
      sourceGridIntervals: nominal, requestedGridIntervals: stageIntervals }, { maxNodes });
    const final = plan.nominalGridIntervals >= requestedGridIntervals
      && plan.child.tubes.every((n, g) => n >= requested.tubes[g]);
    levels.push({ index, final, ...plan, targetWeights });
    if (final) return { kind: 'coupled-grid-levels', sourceGridIntervals, requestedGridIntervals,
      finalNominalGridIntervals: plan.nominalGridIntervals, maxNodes, maxStreamwiseFactor,
      source, requested, child: copy(plan.child), levels, maximumNodes: Math.max(...levels.map(level => level.nodeCount)),
      globalSolvesBeforeFinalSeed: levels.length - 1, finalSeedOnly: true };
    nominal = plan.nominalGridIntervals;
    // Only length and weight data are read by the pure planner. These null
    // coordinates are deliberately not a geometry or a numerical-system input.
    virtualSource = { bodies: sourceInput.bodies, outerLower: Array(plan.child.nx + 1).fill(null),
      weights: virtualSource.weights.map((row, g) => row.flatMap((w, j) =>
        Array(plan.normalSubdivisions[g][j]).fill(w / plan.normalSubdivisions[g][j]))) };
  }
  throw new Error('Coupled grid-level planning exceeded its bounded level count.');
}

const physicalGridKeys = new Set(['weights', 'gridSpacing', 'bodies', 'outerLower', 'outerUpper', 'cutPaths']);
function assertSamePhysics(sourceInput, targetInput, sourceOptions, targetOptions) {
  for (const key of new Set([...Object.keys(sourceInput), ...Object.keys(targetInput)]))
    if (!physicalGridKeys.has(key)) require(same(sourceInput[key], targetInput[key]),
      `Coupled grid level changed physical/model input ${key}.`);
  require(sourceInput.bodies.length === targetInput.bodies.length, 'Coupled grid level changed body count.');
  sourceInput.bodies.forEach((body, b) => {
    const target = targetInput.bodies[b];
    require(body.element === target.element && same(body.points, target.points) && same(body.trailingEdge, target.trailingEdge),
      'Coupled grid level changed a solid contour, trailing edge or explicit body identity.');
  });
  for (const key of ['reynolds', 'ncrit', 'edgeMatching', 'transitionMode', 'blThermodynamics'])
    require(same(sourceOptions[key], targetOptions[key]), `Coupled grid level changed physical/model option ${key}.`);
  require((sourceOptions.hkFloorLinearization ?? 'exact') === (targetOptions.hkFloorLinearization ?? 'exact'),
    'Coupled grid level changed Hk-floor linearization.');
  require((sourceOptions.geometryReplay ?? 'legacy') === (targetOptions.geometryReplay ?? 'legacy'),
    'Coupled grid level changed geometry replay policy.');
  require(same(sourceOptions.tripFractions ?? sourceInput.bodies.map(() => [1, 1]),
    targetOptions.tripFractions ?? targetInput.bodies.map(() => [1, 1])), 'Coupled grid level moved material trips.');
}

function verifyRoot(result, tolerance, suppliedSystem) {
  const cp = result?.checkpoint, saved = cp?.restart;
  require(result?.converged === true && result.mesh?.quality?.valid === true && coupledConvergenceSatisfied(result, tolerance)
    && cp?.version === 1 && saved?.input && saved.options && saved.initialEuler && finiteVector(saved.initialEuler.x)
    && finiteVector(saved.initialBL) && finiteVector(result.x) && finiteVector(result.residual)
    && same(cp.families, result.families), 'Grid levels require a complete converged coupled checkpoint.');
  const system = suppliedSystem ?? createCoupledStreamtubeBody(saved.input, { ...saved.options,
    initialEuler: saved.initialEuler, initialBL: saved.initialBL });
  const value = system.evaluate(system.initial);
  const checks = { state: same(system.initial, result.x), residual: same(value.residual, result.residual),
    families: same(value.families, result.families), finite: value.residual.every(Number.isFinite),
    admissible: system.admissible(system.initial),
    nodes: gridReplayDeparture(value.outer.nodes, saved.initialEuler.nodes),
    undisplacedNodes: gridReplayDeparture(value.outer.undisplacedNodes, saved.initialEuler.undisplacedNodes) };
  const failed = Object.entries(checks).filter(([, check]) => typeof check === 'boolean' ? !check : !check.equivalent)
    .map(([key, check]) => typeof check === 'boolean' ? key : `${key} (maximum coordinate change ${check.maximum})`);
  if (failed.length)
    throw Object.assign(new Error(`Grid-level source state, residual or physical grid does not replay exactly: ${failed.join(', ')}.`),
      { code: 'GRID_LEVEL_REPLAY', diagnostics: checks });
  require(system.bl.transitionMode === 'automatic' && system.bl.trips.every(pair => pair.every(v => v === 1))
    && same(system.bl.snapshotActive(), saved.options.transitionState)
    && same(system.bl.snapshotActive(), result.boundaryLayer?.transitionState),
  'Grid-level source automatic transition phase or terminal trips are inconsistent.');
  require(streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes }).quality.valid,
    'Grid-level source physical mesh is invalid.');
  return system;
}

// solveIntermediateLevel is caller-owned ordinary coupled solving, with the
// caller's existing budgets, observers and cancellation. It returns a result
// (or {result, failure}); thrown values ALWAYS propagate unchanged. Expected
// numerical failure should be returned explicitly, never thrown by a wrapper
// that also handles observer cancellation. Final preparation never calls it.
export function prepareCoupledGridLevels({ sourceSystem, sourceResult, input, options,
  requestedGridIntervals, sourceGridIntervals } = {}, { tolerance = 1e-10,
  maxNodes = COUPLED_GRID_SEQUENCE_MAX_NODES, maxStreamwiseFactor = 2,
  blPredictor = 'interpolate', solveIntermediateLevel, onLevel } = {}) {
  require(Number.isFinite(tolerance) && tolerance > 0, 'Invalid coupled grid-level tolerance.');
  require(['interpolate', 'xfoil-mrchdu'].includes(blPredictor), 'Invalid coupled grid-level BL predictor.');
  require(onLevel === undefined || typeof onLevel === 'function', 'Invalid grid-level observer.');
  // Preflight the WHOLE requested hierarchy before even constructing a source
  // numerical system. A final oversized grid cannot waste intermediate solves.
  const plan = planCoupledGridLevels({ sourceInput: sourceResult?.checkpoint?.restart?.input,
    input, requestedGridIntervals, sourceGridIntervals }, { maxNodes, maxStreamwiseFactor });
  require(plan.levels.length === 1 || typeof solveIntermediateLevel === 'function',
    'Multiple grid levels require an intermediate coupled solve callback.');
  let retained = sourceResult, retainedNominal = sourceGridIntervals;
  let system = verifyRoot(retained, tolerance, sourceSystem);
  const saved = retained.checkpoint.restart;
  assertSamePhysics(saved.input, input, saved.options, options);
  const levels = [], emit = event => onLevel?.(copy(event));
  const finish = (prepared, failure) => ({ prepared, sourceResult: retained, sourceSystem: system,
    sourceGridIntervals: retainedNominal, diagnostics: { method: 'converged-coupled-grid-levels',
      requestedGridIntervals, actualGridIntervals: retainedNominal, actual: dimensions(retained.checkpoint.restart.input),
      finalSeedPrepared: prepared !== null, reachedTarget: false, stateConverged: true,
      equationsChanged: false, plan: copy(plan), levels: copy(levels), ...(failure ? { failure: copy(failure) } : {}) } });
  for (const planned of plan.levels) {
    const record = { index: planned.index, final: planned.final, sourceGridIntervals: retainedNominal,
      requestedGridIntervals: planned.requestedGridIntervals, nominalGridIntervals: planned.nominalGridIntervals,
      source: dimensions(retained.checkpoint.restart.input), prepared: false, converged: false };
    levels.push(record);
    emit({ stage: 'preparing', ...record });
    let prepared;
    try {
      prepared = prepareCoupledGridSequence({ sourceSystem: system, sourceResult: retained,
        input: { ...input, weights: copy(planned.targetWeights) }, options,
        requestedGridIntervals: planned.requestedGridIntervals, sourceGridIntervals: retainedNominal },
        { tolerance, maxNodes, ...(blPredictor === 'xfoil-mrchdu' ? { blPredictor } : {}) });
      require(prepared.transfer.nominalGridIntervals === planned.nominalGridIntervals,
        'Grid-level preparation changed planned nominal resolution.');
      record.prepared = true; record.actual = dimensions(prepared.input);
      record.transfer = copy(prepared.transfer);
    } catch (error) {
      const failure = { stage: 'preparation', level: planned.index, reason: error?.message ?? String(error),
        ...(error?.code ? { code: error.code } : {}),
        ...(error?.diagnostics === undefined ? {} : { diagnostics: copy(error.diagnostics) }) };
      emit({ stage: 'retained', ...record, failure });
      return finish(null, failure);
    }
    emit({ stage: 'prepared', ...record });
    if (planned.final) return finish(prepared);
    const output = solveIntermediateLevel(prepared, copy(record));
    const candidate = output?.result ?? output;
    let nextSystem;
    try {
      if (output?.failure) throw Object.assign(new Error(output.failure.reason ?? output.failure.message ?? 'Intermediate coupled solve failed.'),
        { code: output.failure.code, diagnostics: output.failure.diagnostics });
      if (!candidate?.converged) throw new Error(candidate?.reason ?? 'Intermediate coupled grid did not converge.');
      assertSamePhysics(saved.input, candidate.checkpoint?.restart?.input ?? {}, saved.options, candidate.checkpoint?.restart?.options ?? {});
      require(coupledCheckpointShearCoordinate(candidate.checkpoint) === coupledFreshShearCoordinate(prepared),
        'Intermediate coupled grid changed its inherited shear-coordinate policy.');
      nextSystem = verifyRoot(candidate, tolerance);
      require(same(nextSystem.conditions, prepared.system.conditions)
        && same(nextSystem.euler.conditions, prepared.system.euler.conditions),
      'Intermediate coupled grid changed effective normalized physics.');
      require(candidate.checkpoint.restart.input.outerLower.length === prepared.input.outerLower.length
        && same(candidate.checkpoint.restart.input.weights, prepared.input.weights),
      'Intermediate coupled grid changed its prepared topology or captured weights.');
    } catch (error) {
      const failure = { stage: 'intermediate-solve', level: planned.index, reason: error?.message ?? String(error),
        ...(error?.code ? { code: error.code } : {}),
        ...(error?.diagnostics === undefined ? {} : { diagnostics: copy(error.diagnostics) }) };
      record.reason = failure.reason; record.families = copy(candidate?.families);
      emit({ stage: 'retained', ...record, failure });
      return finish(null, failure);
    }
    retained = candidate; retainedNominal = prepared.transfer.nominalGridIntervals; system = nextSystem;
    record.converged = true; record.families = copy(candidate.families);
    record.iterations = Math.max(0, (candidate.history?.length ?? 1) - 1);
    emit({ stage: 'converged', ...record });
  }
  throw new Error('Grid-level preparation completed without a final seed.');
}

// Coordinate reconstruction can subtract and restore a BL displacement.
// Allow only arithmetic roundoff, never a state or residual tolerance.
export function gridReplayDeparture(actual, saved) {
  let maximum = 0, equivalent = true;
  const visit = (a, b) => {
    if (Array.isArray(a)) {
      if (!Array.isArray(b) || a.length !== b.length) { equivalent = false; return; }
      a.forEach((value, i) => visit(value, b[i]));
      return;
    }
    for (const key of ['x', 'y']) {
      const x = a?.[key], y = b?.[key];
      if (!Number.isFinite(x) || !Number.isFinite(y)) { equivalent = false; return; }
      const difference = Math.abs(x - y);
      maximum = Math.max(maximum, difference);
      if (difference > 8 * Number.EPSILON * Math.max(1, Math.abs(x), Math.abs(y))) equivalent = false;
    }
  };
  visit(actual, saved);
  return { equivalent, maximum };
}

// Solve one transferred level of the same coupled equations. Each published
// pressure field uses the checkpoint and station topology of that level.

export function solveCoupledGridLevel(prepared, level, { maxIterations, tolerance, normalization,
  startupAttempt, dissipationEnhancement = false, iterationRecovery, onStage, onIteration, onMesh, onFlow, onIterationCheckpoint } = {}) {
  const stage = 'coupled-grid-refinement';
  const labels = { stage, startupAttempt, gridLevel: level.nominalGridIntervals,
    requestedGridIntervals: level.targetGridIntervals ?? level.requestedGridIntervals,
    actualNcrit: prepared.options.ncrit, targetNcrit: level.targetNcrit ?? prepared.options.ncrit };
  let checkpoint, observerFailed = false, observerError;
  const observe = callback => (...args) => {
    try { return callback?.(...args); }
    catch (error) { observerFailed = true; observerError = error; throw error; }
  };
  const publishMesh = observe(onMesh), publishFlow = observe(onFlow);
  try {
    observe(onStage)(labels);
    const result = solveCoupledStreamtubeIses(prepared.input, { ...prepared.options,
      initialEuler: prepared.initialEuler, initialBL: prepared.initialBL,
      maxIterations, tolerance, dissipationEnhancement, iterationRecovery, iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible',
      blUpdate: 'xfoil', projectionGeometry: 'boundary-increment',
      shearCoordinate: coupledFreshShearCoordinate(prepared),
      onIteration: observe(h => onIteration?.({ ...h, ...labels })),
      onCheckpoint: observe((value, details) => {
        checkpoint = value;
        onIterationCheckpoint?.(structuredClone(value), { ...details, ...labels });
      }),
      onMesh: state => {
        const mesh = streamtubeMeshSnapshot(state);
        mesh.initialization.gridRefinement = prepared.transfer;
        publishMesh(mesh, state.iteration.iteration === 0 ? 'initial' : 'solving', stage, state);
        if (onFlow && checkpoint) publishFlow(structuredClone({ checkpoint,
          flow: observableFlow(state.flow),
          bl: observableBL(prepared.system.bl),
          bodies: state.system.layout.bodies, normalization, iteration: state.iteration,
          mach: checkpoint.restart.input.mach, ...labels }));
      } });
    result.mesh.initialization.gridRefinement = prepared.transfer;
    return result;
  } catch (error) {
    if (observerFailed) throw observerError;
    return { result: null, failure: { reason: error?.message ?? String(error),
      code: error?.code, diagnostics: error?.diagnostics } };
  }
}
