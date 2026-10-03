// SPDX-License-Identifier: GPL-2.0-or-later
// Native profile preparation after complete grid transfer.
// This supplies a new initial guess; it performs no coupled Newton solve.
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { prepareCoupledMrchduProfiles } from './streamtube-coupled-mrchdu-predictor.js';
import { extendWarmBoundaryIncrements } from './streamtube-displacement.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { assertConvexStreamtubeGrid } from '../geometry/streamtube-convex-step.js';
import { incrementIndependentWakeWidths } from './streamtube-wake-geometry.js';

const same = (a, b) => {
  if (a === b) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  if (Array.isArray(a) || ArrayBuffer.isView(a)) return (Array.isArray(b) || ArrayBuffer.isView(b))
    && a.length === b.length && Array.from(a).every((value, i) => same(value, b[i]));
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(key => Object.hasOwn(b, key) && same(a[key], b[key]));
};
export function prepareCoupledGridProfile(prepared, { wakeWidthIncrement = false, reinitializeAmplification = false } = {}) {
  let stage = 'source replay';
  const diagnostics = { method: 'xfoil-mrchdu-after-grid-transfer', initialGuessOnly: true, flowSolved: false, converged: false,
    physicalAcceptance: false, equationsChanged: false, warnings: [],
    operations: { constructors: 0, fullCoupledEvaluations: 0, strictDomainChecks: 0,
      nativeProfileCalls: 0, globalNewtonUpdates: 0, globalJacobians: 0, globalLinearSolves: 0 } };
  const require = (ok, message) => {
    if (!ok) throw Object.assign(new Error(message), { code: 'COUPLED_GRID_PROFILE_PREPARATION',
      diagnostics: { ...structuredClone(diagnostics), stage } });
  };
  const equal = (a, b, label) => {
    if (!same(a, b)) { diagnostics.failedInvariant = label;
      if (label === 'source physical nodes') {
        // This compares reconstruction with the supplied chart; it is not
        // evidence of mutation of the caller's (cloned) source arrays.
        let changedCoordinates = 0, maximumAbsoluteDifference = 0, firstDifference;
        a.forEach((group, g) => group.forEach((row, i) => row.forEach((point, j) => {
          for (const coordinate of ['x', 'y']) {
            const supplied = b[g]?.[i]?.[j]?.[coordinate], replayed = point[coordinate];
            if (supplied !== replayed) {
              changedCoordinates++;
              firstDifference ??= { group: g, station: i, node: j, coordinate, supplied, replayed };
              maximumAbsoluteDifference = Math.max(maximumAbsoluteDifference, Math.abs(replayed - supplied));
            }
          }
        })));
        diagnostics.geometryReplayDifference = { changedCoordinates, maximumAbsoluteDifference, firstDifference,
          comparison: 'reconstructed versus supplied physical grid', callerSourceCloned: true };
      }
      require(false, `Grid profile preparation changed ${label}.`); }
  };
  const construct = f => { diagnostics.operations.constructors++;
    return createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL }); };
  const evaluate = (system, state) => { diagnostics.operations.fullCoupledEvaluations++; return system.evaluate(state); };
  const domain = (system, state) => {
    diagnostics.operations.strictDomainChecks++;
    let failure;
    require(system.admissible(state, { requireConvex: true, onFailure: value => { failure = value; } }),
      `Grid profile state is outside the physical/convex domain${failure ? `: ${JSON.stringify(failure)}` : '.'}`);
  };
  try {
    require(typeof wakeWidthIncrement === 'boolean' && typeof reinitializeAmplification === 'boolean',
      'Grid-profile preparation controls must be boolean.');
    require(prepared?.input && prepared.options?.transitionMode === 'automatic'
      && prepared.initialEuler?.x && prepared.initialEuler.nodes && prepared.initialBL,
    'Grid profile preparation requires a complete automatic coupled initial guess.');
    require(!wakeWidthIncrement || prepared.input.wakeGeometry === 'independent-banks'
      && prepared.input.wakeDisplacementMotion === 'te-center',
    'Wake-width increments require independent TE-center wake banks.');
    const original = { input: structuredClone(prepared.input), options: structuredClone(prepared.options),
      initialEuler: structuredClone(prepared.initialEuler), initialBL: structuredClone(prepared.initialBL) };
    const source = construct(original), before = evaluate(source, source.initial); domain(source, source.initial);
    equal(source.initial, [...original.initialEuler.x, ...original.initialBL], 'source packed state');
    equal(before.outer.nodes, original.initialEuler.nodes, 'source physical nodes');
    if (prepared.value) equal(before.residual, prepared.value.residual, 'source residual');
    diagnostics.beforeFamilies = { ...before.families }; diagnostics.beforeTransitionState = source.bl.snapshotActive();

    stage = 'native profile preparation';
    const prediction = prepareCoupledMrchduProfiles({ bl: source.bl, states: before.layers.states,
      initialBL: original.initialBL, targetMach: original.input.mach });
    diagnostics.prediction = prediction.diagnostics;
    diagnostics.operations.nativeProfileCalls = prediction.diagnostics.operations.translatedMRCHDUCalls;
    diagnostics.warnings = prediction.diagnostics.bodies.flatMap(body => body.localConvergenceWarnings);
    require(diagnostics.warnings.length === 0, 'Native grid-profile preparation did not converge locally.');
    const restart = { ...structuredClone(original), options: { ...structuredClone(original.options),
      transitionState: prediction.transitionState.slice() }, initialBL: Array.from(prediction.initialBL) };
    let target = construct(restart), state = target.initial.slice();
    stage = 'native phase reconciliation';
    diagnostics.phasePreparation = target.bl.updateActive(state.subarray(target.ne), state.subarray(0, target.ne),
      ...(reinitializeAmplification ? [{ reinitializeAmplification: true }] : []));
    restart.options.transitionState = target.bl.snapshotActive(); restart.initialBL = Array.from(state.subarray(target.ne));

    stage = 'boundary increment transfer';
    target.euler.setDisplacement(target.bl.thicknesses(state.subarray(target.ne)));
    const decoded = target.euler.decode(state.subarray(0, target.ne));
    const wake = wakeWidthIncrement ? incrementIndependentWakeWidths({ layout: target.euler.layout, nodes: decoded.nodes,
      beforeWidths: source.bl.thicknesses(source.initial.subarray(source.ne)).wakes,
      afterWidths: target.bl.thicknesses(state.subarray(target.ne)).wakes }) : null;
    const nodes = extendWarmBoundaryIncrements({ sourceNodes: before.outer.nodes, targetNodes: wake?.nodes ?? decoded.nodes,
      masses: before.outer.allocation.groups.map(group => group.map(tube => tube.massFlow)) });
    state.set(target.euler.adoptGeometry(state.subarray(0, target.ne), nodes));
    diagnostics.geometry = { method: 'existing-mass-coordinate-boundary-increments',
      fullDisplacementReapplied: false, independentWakeBanksAveraged: false,
      ...(wake ? { wakeWidthIncrement: wake.diagnostics } : {}) };

    const preserved = (system, value) => {
      equal(system.conditions, source.conditions, 'effective coupled conditions');
      equal(system.euler.conditions, source.euler.conditions, 'effective Euler conditions');
      equal(system.bl.kernel.parameters, source.bl.kernel.parameters, 'BL closure parameters');
      equal(system.bl.trips, source.bl.trips, 'material trips');
      for (const key of ['allocation', 'captured', 'stagnation', 'strengths'])
        equal(value.outer[key], before.outer[key], `physical Euler ${key}`);
      equal(value.outer.sections.map(row => row.map(group => group.map(section => section.rho))),
        before.outer.sections.map(row => row.map(group => group.map(section => section.rho))), 'physical Euler density');
      require(system.n === source.n && system.ne === source.ne && value.layers.states.length === before.layers.states.length,
        'Grid profile preparation changed coupled topology.');
      assertConvexStreamtubeGrid(value.outer.nodes);
    };
    stage = 'complete target domain'; domain(target, state);
    let value = evaluate(target, state); preserved(target, value);
    restart.initialEuler = { x: Array.from(state.slice(0, target.ne)), nodes: value.outer.nodes,
      undisplacedNodes: value.outer.undisplacedNodes };
    stage = 'canonical replay';
    let canonical = construct(structuredClone(restart)); domain(canonical, canonical.initial);
    let canonicalValue = evaluate(canonical, canonical.initial);
    equal(canonical.initial, state, 'canonical packed state');
    equal(canonicalValue.residual, value.residual, 'canonical full residual'); preserved(canonical, canonicalValue);
    restart.initialEuler.nodes = canonicalValue.outer.nodes;
    restart.initialEuler.undisplacedNodes = canonicalValue.outer.undisplacedNodes;
    target = construct(structuredClone(restart)); domain(target, target.initial);
    value = evaluate(target, target.initial);
    equal(target.initial, canonical.initial, 'second canonical packed state');
    equal(value.residual, canonicalValue.residual, 'second canonical full residual'); preserved(target, value);
    diagnostics.afterFamilies = { ...value.families }; diagnostics.afterTransitionState = target.bl.snapshotActive();
    diagnostics.canonicalReplayExact = true; diagnostics.physicalEulerInvariantsPreserved = true;
    const transfer = { ...structuredClone(prepared.transfer ?? {}), targetFamilies: { ...value.families },
      targetTransitionState: target.bl.snapshotActive(), profilePreparation: diagnostics };
    const mesh = streamtubeMeshSnapshot({ system: target.euler, nodes: value.outer.nodes });
    require(mesh.quality.valid, 'Prepared native grid-profile mesh is invalid.');
    mesh.initialization.gridRefinement = structuredClone(transfer);
    return { ...prepared, system: target, input: restart.input, options: restart.options,
      initialEuler: restart.initialEuler, initialBL: restart.initialBL, value, mesh, transfer,
      initialization: { ...structuredClone(prepared.initialization ?? {}), transfer: structuredClone(transfer),
        profilePreparation: structuredClone(diagnostics) } };
  } catch (error) {
    if (error?.code === 'COUPLED_GRID_PROFILE_PREPARATION') throw error;
    throw Object.assign(new Error(error?.message ?? String(error), { cause: error }), {
      code: 'COUPLED_GRID_PROFILE_PREPARATION', diagnostics: { ...structuredClone(diagnostics), stage,
        cause: { message: error?.message ?? String(error), code: error?.code, diagnostics: error?.diagnostics } } });
  }
}

// Only for a NEW locally refined initial guess, never a retained root.
// Preserve the successful preparation path exactly. If its saved physical
// grid cannot replay, try retaining the supplied raw chart while refreshing
// physical NCALC directions. All original exact replay/domain gates remain.

export function prepareRefinedCoupledGridProfile(prepared, options) {
  try { return prepareCoupledGridProfile(prepared, options); }
  catch (error) {
    if (error?.code !== 'COUPLED_GRID_PROFILE_PREPARATION'
      || error.diagnostics?.stage !== 'source replay'
      || error.diagnostics.failedInvariant !== 'source physical nodes'
      || prepared.options?.geometryReplay === 'preserve-undisplaced') throw error;
    const recovery = { policy: 'preserve-undisplaced', attempted: true, accepted: false,
      sourceFailure: { message: error.message, code: error.code, diagnostics: structuredClone(error.diagnostics) },
      exactSourceGeometryRequired: true, originalGatesChanged: false };
    // Prepared grid sequences also carry live system/flow methods. Clone
    // the restart data we change, not those non-serializable runtime objects.
    const candidate = { ...prepared, input: structuredClone(prepared.input),
      initialEuler: structuredClone(prepared.initialEuler), initialBL: structuredClone(prepared.initialBL),
      options: { ...structuredClone(prepared.options), geometryReplay: 'preserve-undisplaced' } };
    try {
      const result = prepareCoupledGridProfile(candidate, options);
      recovery.accepted = true;
      for (const diagnostics of [result.transfer.profilePreparation, result.initialization?.profilePreparation,
        result.initialization?.transfer?.profilePreparation, result.mesh?.initialization?.gridRefinement?.profilePreparation])
        if (diagnostics) diagnostics.geometryReplayRecovery = structuredClone(recovery);
      return result;
    } catch (retryError) {
      // Retain original diagnostics alongside the rejected exact retry.
      // This numerical helper accepts no observers and handles no callbacks.
      if (retryError?.code === 'COUPLED_GRID_PROFILE_PREPARATION')
        retryError.diagnostics = { ...retryError.diagnostics, geometryReplayRecovery: recovery };
      throw retryError;
    }
  }
}
