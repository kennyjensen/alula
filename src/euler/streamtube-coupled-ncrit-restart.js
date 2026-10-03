// SPDX-License-Identifier: GPL-2.0-or-later
// Browser-compatible complete-state
// Ncrit continuation preparation; never a coupled solve or convergence claim.
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { coupledConvergenceSatisfied } from './streamtube-coupled-convergence.js';
import { assertConvexStreamtubeGrid } from '../geometry/streamtube-convex-step.js';
import { prepareCoupledMrchduProfiles } from './streamtube-coupled-mrchdu-predictor.js';
import { extendWarmBoundaryIncrements } from './streamtube-displacement.js';

export class CoupledNcritRestartError extends Error {
  constructor(message, diagnostics) {
    super(message); this.name = 'CoupledNcritRestartError';
    this.code = 'COUPLED_NCRIT_RESTART'; this.diagnostics = diagnostics;
  }
}
const vector = a => Array.isArray(a) || ArrayBuffer.isView(a) && !(a instanceof DataView);
const finiteVector = a => vector(a) && a.length > 0 && Array.from(a).every(Number.isFinite);
const same = (a, b) => {
  if (a === b) return true; // JSON does not retain signed zero.
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  if (vector(a) || vector(b)) return vector(a) && vector(b) && a.length === b.length
    && Array.from(a).every((v, i) => same(v, b[i]));
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};
const maximum = a => { let v = 0; for (const x of a) v = Math.max(v, Math.abs(x)); return v; };

/** Prepare a strictly increasing Ncrit initial guess at unchanged physics.
 * Ncrit affects the amplification criterion, not gas reference or displacement.
 * Consequently theta, delta*, Ue, all Euler coordinates and maintenance history
 * remain exact. Only the existing automatic-transition N/Ctau handoff may alter
 * auxiliaries by default. Explicit xfoil-mrchdu instead makes one native profile
 * prediction and extends only its boundary increments into the existing grid.
 * That candidate may change physical profiles and geometry, and native local
 * warnings remain explicit. Either candidate must be solved/qualified before
 * being used as the next converged source.
 */
export function prepareCoupledNcritRestart(targetNcrit, checkpoint, {
  requestedNcrit = targetNcrit, tolerance = 1e-10, blPredictor = 'preserve',
} = {}) {
  let stage = 'source validation';
  const diagnostics = { sourceConverged: false, targetConverged: false,
    selectedNcrit: targetNcrit, requestedNcrit, intermediate: targetNcrit !== requestedNcrit,
    provisional: true, physicalAcceptance: false,
    method: blPredictor === 'xfoil-mrchdu' ? 'xfoil-mrchdu' : 'preserve-physical-profile', warnings: [],
    operations: { constructors: 0, fullCoupledEvaluations: 0, strictDomainChecks: 0,
      nativeProfileCalls: 0, globalJacobians: 0, globalLinearSolves: 0, globalNewtonUpdates: 0 } };
  const require = (condition, message) => {
    if (!condition) throw new CoupledNcritRestartError(message, { ...diagnostics, stage });
  };
  const equal = (a, b, label) => require(same(a, b), `Ncrit restart changed ${label}.`);
  const construct = r => {
    diagnostics.operations.constructors++;
    return createCoupledStreamtubeBody(r.input, { ...r.options,
      initialEuler: r.initialEuler, initialBL: r.initialBL });
  };
  const evaluate = (system, state) => {
    diagnostics.operations.fullCoupledEvaluations++; return system.evaluate(state);
  };
  const domain = (system, state) => {
    diagnostics.operations.strictDomainChecks++;
    let failure;
    const accepted = system.admissible(state, { requireConvex: true, onFailure: f => { failure = f; } });
    if (!accepted) diagnostics.domainFailure = failure;
    require(accepted, 'Ncrit restart fails complete physical or convex-grid admissibility.');
  };
  try {
    require(['preserve', 'xfoil-mrchdu'].includes(blPredictor), 'Unknown Ncrit restart BL predictor.');
    require(Number.isFinite(tolerance) && tolerance > 0, 'Invalid Ncrit restart residual tolerance.');
    require(checkpoint?.version === 1 && checkpoint.restart && checkpoint.continuation && checkpoint.families,
      'Ncrit restart requires a complete version-1 coupled ISES checkpoint.');
    const saved = structuredClone(checkpoint), f = saved.restart, c = saved.continuation;
    require(f.input && f.options && f.initialEuler && finiteVector(f.initialEuler.x) && finiteVector(f.initialBL)
      && Array.isArray(f.initialEuler.nodes) && Array.isArray(f.initialEuler.undisplacedNodes),
    'Ncrit restart requires complete finite Euler/BL and physical/undisplaced grid arrays.');
    require(Number.isFinite(f.options.ncrit) && Number.isFinite(targetNcrit) && Number.isFinite(requestedNcrit)
      && targetNcrit > f.options.ncrit && targetNcrit <= requestedNcrit,
    'Ncrit restart must increase the source criterion toward the unchanged requested target.');
    diagnostics.sourceNcrit = f.options.ncrit;
    require(f.options.transitionMode === 'automatic' && Array.isArray(f.options.transitionState)
      && Array.isArray(f.options.tripFractions) && f.options.tripFractions.length === f.input.bodies?.length
      && f.options.tripFractions.every(p => Array.isArray(p) && p.length === 2 && p.every(t => t === 1)),
    'Ncrit restart requires complete automatic transition with terminal material trips.');
    const allowed = ['reynolds', 'ncrit', 'tripFractions', 'edgeMatching', 'transitionMode',
      'transitionState', 'blThermodynamics', 'hkFloorLinearization', 'geometryReplay'];
    require(Object.keys(f.options).every(k => allowed.includes(k)), 'Unsupported coupled checkpoint options.');
    require(coupledConvergenceSatisfied(saved, tolerance), 'Ncrit source must already satisfy its coupled convergence criterion.');
    require(['convex', 'ises-sampled'].includes(c.iterationGeometry)
      && ['listing', 'admissible', 'armijo', 'event-armijo'].includes(c.stepAcceptance) && ['listing', 'prose'].includes(c.stagnationLimiter)
      && ['amd', 'colamd'].includes(c.preferredOrdering) && [.001, 1].includes(c.pivotTolerance)
      && finiteVector(c.lastRedistributedStagnation) && c.lastRedistributedStagnation.length === f.input.bodies.length
      && Array.isArray(c.fractions) && c.fractions.length === f.input.bodies.length
      && c.fractions.every((a, b) => Array.isArray(a) && a.length === f.input.bodies[b].leadingIndex + 1
        && a[0] === 0 && a.at(-1) === 1 && a.every((x, i) => Number.isFinite(x) && (!i || x > a[i - 1]))),
    'Ncrit restart requires complete ISES maintenance history.');
    require(f.input.stagnationMotion === 'walls-only' && f.input.normalStencil !== undefined
      && f.input.geometryDomain === (c.iterationGeometry === 'convex' ? 'convex' : 'positive-simple'),
    'Ncrit restart requires normalized ISES geometry controls.');
    assertConvexStreamtubeGrid(f.initialEuler.nodes);
    const source = construct(f), before = evaluate(source, source.initial); domain(source, source.initial);
    equal(source.initial, [...f.initialEuler.x, ...f.initialBL], 'source packed state during replay');
    equal(before.families, saved.families, 'source residual families during replay');
    equal(before.outer.nodes, f.initialEuler.nodes, 'source physical grid during replay');
    equal(source.bl.snapshotActive(), f.options.transitionState, 'source transition map during replay');
    require(coupledConvergenceSatisfied({ ...saved, families: before.families }, tolerance), 'Replayed Ncrit source has not converged.');
    diagnostics.sourceConverged = true; diagnostics.sourceFamilies = { ...before.families };
    diagnostics.sourcePhase = source.bl.snapshotActive();

    stage = 'target auxiliary phase preparation';
    const restart = structuredClone(f); restart.options.ncrit = targetNcrit;
    let target = construct(restart), initial = target.initial.slice();
    equal(initial, source.initial, 'physical seed during target kernel construction');
    if (blPredictor === 'xfoil-mrchdu') {
      stage = 'target native profile preparation';
      const prediction = prepareCoupledMrchduProfiles({ bl: target.bl, states: before.layers.states,
        initialBL: f.initialBL, targetMach: f.input.mach });
      diagnostics.operations.nativeProfileCalls = prediction.diagnostics.operations.translatedMRCHDUCalls;
      diagnostics.prediction = prediction.diagnostics;
      diagnostics.warnings = prediction.diagnostics.bodies.flatMap(b => b.localConvergenceWarnings);
      diagnostics.warningFree = diagnostics.warnings.length === 0;
      restart.options.transitionState = prediction.transitionState;
      restart.initialBL = prediction.initialBL;
      target = construct(restart); initial = target.initial.slice();
    }
    stage = 'target auxiliary phase preparation';
    const beforePhase = target.bl.snapshotActive();
    const targets = target.bl.activeTargets(initial.subarray(0, target.ne), initial.subarray(target.ne));
    const reconcileAmplificationSurfaces = targets.flatMap((t, k) => t.amplificationReconciliation ? [k] : []);
    diagnostics.targetPhaseInitialization = {
      method: 'existing-automatic-transition-transfer', targets,
      ...target.bl.updateActive(initial.subarray(target.ne), initial.subarray(0, target.ne),
        blPredictor === 'preserve' ? { reconcileAmplificationSurfaces } : {}),
      before: beforePhase, after: target.bl.snapshotActive(),
    };
    restart.options.transitionState = target.bl.snapshotActive();
    restart.initialBL = Array.from(initial.subarray(target.ne));
    diagnostics.auxiliaryChanges = target.bl.stations.flatMap(({ id }) => {
      const before = source.initial[source.ne + 4 * id], after = initial[target.ne + 4 * id];
      return before === after ? [] : [{ id, before, after }];
    });
    if (blPredictor === 'xfoil-mrchdu') {
      stage = 'native profile boundary-increment transfer';
      target.euler.setDisplacement(target.bl.thicknesses(initial.subarray(target.ne)));
      const decoded = target.euler.decode(initial.subarray(0, target.ne));
      const nodes = extendWarmBoundaryIncrements({ sourceNodes: before.outer.nodes, targetNodes: decoded.nodes,
        masses: before.outer.allocation.groups.map(g => g.map(t => t.massFlow)) });
      initial.set(target.euler.adoptGeometry(initial.subarray(0, target.ne), nodes));
      diagnostics.boundaryIncrementTransfer = { method: 'existing-mass-coordinate-boundary-increments',
        fullDisplacementReapplied: false, independentWakeBanksAveraged: false };
    }
    const preserved = (system, state, value) => {
      equal(system.euler.conditions, source.euler.conditions, 'Euler equation/gas conditions');
      equal(system.conditions, { ...source.conditions, ncrit: targetNcrit }, 'non-Ncrit coupled conditions');
      equal(system.bl.kernel.parameters, { ...source.bl.kernel.parameters, ncrit: targetNcrit }, 'non-Ncrit integral closure parameters');
      for (const key of ['allocation', 'captured', 'stagnation', 'strengths'])
        equal(value.outer[key], before.outer[key], `physical Euler ${key}`);
      equal(value.outer.sections.map(row => row.map(g => g.map(s => s.rho))),
        before.outer.sections.map(row => row.map(g => g.map(s => s.rho))), 'physical Euler densities');
      require(value.layers.states.length === before.layers.states.length, 'Ncrit restart changed the BL topology.');
      if (blPredictor === 'preserve') {
        equal(state.subarray(0, system.ne), source.initial.subarray(0, source.ne), 'stored Euler coordinates');
        for (const key of ['sections', 'nodes', 'undisplacedNodes', 'residual'])
          equal(value.outer[key], before.outer[key], `physical Euler ${key}`);
        value.layers.states.forEach((s, id) => {
          for (const key of ['theta', 'deltaStar', 'ue', 's', 'wakeGap']) equal(s[key], before.layers.states[id][key], `BL ${key} at ${id}`);
        });
      }
    };
    stage = 'strict target evaluation'; domain(target, initial);
    const value = evaluate(target, initial); preserved(target, initial, value);
    stage = 'canonical checkpoint replay';
    if (blPredictor === 'xfoil-mrchdu') restart.initialEuler = { x: initial.slice(0, target.ne),
      nodes: value.outer.nodes, undisplacedNodes: value.outer.undisplacedNodes };
    const resultCheckpoint = { ...saved, families: { ...value.families }, restart };
    delete resultCheckpoint.convergence; // New physics needs fresh change evidence.
    const replay = construct(structuredClone(restart)); domain(replay, replay.initial);
    const replayValue = evaluate(replay, replay.initial);
    equal(replay.initial, initial, 'canonical packed state'); equal(replayValue.residual, value.residual, 'canonical full residual');
    preserved(replay, replay.initial, replayValue);
    if (blPredictor === 'xfoil-mrchdu') {
      // Match the measured native preparation's canonical restart convention.
      // Its decoded reference geometry can differ by roundoff from the input
      // chart; retain the canonical chart and prove the complete replay again.
      restart.initialEuler.nodes = replayValue.outer.nodes;
      restart.initialEuler.undisplacedNodes = replayValue.outer.undisplacedNodes;
      resultCheckpoint.families = { ...replayValue.families };
      const canonical = construct(structuredClone(restart));
      const canonicalValue = evaluate(canonical, canonical.initial);
      equal(canonical.initial, replay.initial, 'second canonical packed state');
      equal(canonicalValue.residual, replayValue.residual, 'second canonical residual');
      equal(canonicalValue.outer.nodes, replayValue.outer.nodes, 'second canonical physical grid');
      preserved(canonical, canonical.initial, canonicalValue);
    }
    const quality = assertConvexStreamtubeGrid(replayValue.outer.nodes);
    diagnostics.targetFamilies = { ...replayValue.families }; diagnostics.targetResidual = maximum(replayValue.residual);
    diagnostics.residualSatisfiesTolerance = diagnostics.targetResidual <= tolerance;
    diagnostics.targetPhase = replay.bl.snapshotActive(); diagnostics.quality = quality;
    diagnostics.physicalBLPreserved = same(initial.subarray(target.ne).filter((_, k) => k % 4 !== 0),
      source.initial.subarray(source.ne).filter((_, k) => k % 4 !== 0));
    diagnostics.physicalEulerPreserved = same(initial.subarray(0, target.ne), source.initial.subarray(0, source.ne))
      && same(replayValue.outer.nodes, before.outer.nodes);
    if (blPredictor === 'xfoil-mrchdu') {
      diagnostics.physicalDensityPreserved = true; diagnostics.physicalMassPreserved = true;
      diagnostics.physicalGlobalsPreserved = true;
    }
    diagnostics.exactReplay = true; diagnostics.maintenanceHistoryPreserved = true;
    diagnostics.prepared = true;
    return { checkpoint: resultCheckpoint,
      seed: { ...structuredClone(restart), checkpoint: structuredClone(resultCheckpoint), families: { ...replayValue.families },
        sourceNcrit: f.options.ncrit, selectedNcrit: targetNcrit, requestedNcrit,
        intermediate: diagnostics.intermediate, sourceConverged: true, converged: false, physicalAcceptance: false },
      diagnostics };
  } catch (error) {
    if (error instanceof CoupledNcritRestartError) throw error;
    throw new CoupledNcritRestartError(error.message, { ...diagnostics, stage,
      cause: { name: error.name, code: error.code, diagnostics: error.diagnostics } });
  }
}
