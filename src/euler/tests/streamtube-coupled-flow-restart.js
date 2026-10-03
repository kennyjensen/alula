// SPDX-License-Identifier: GPL-2.0-or-later
// Complete accepted ISES Euler/BL checkpoint transfer between Mach conditions.
// This changes the gas reference, not the physical density or stored BL state.
// No cold initializer, density inversion, SMOVE or Newton operation is used.
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';
import { prepareCoupledMrchduProfiles } from '../streamtube-coupled-mrchdu-predictor.js';
import { extendWarmBoundaryIncrements } from '../streamtube-displacement.js';

const require = (condition, message) => { if (!condition) throw new Error(message); };
const vector = value => Array.isArray(value) || ArrayBuffer.isView(value) && !(value instanceof DataView);
const same = (a, b) => {
  if (a === b) return true; // JSON checkpoints do not retain signed zero.
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  if (vector(a) || vector(b)) return vector(a) && vector(b) && a.length === b.length
    && Array.from(a).every((v, i) => same(v, b[i]));
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};
const maximum = values => values.reduce((peak, value) => Math.max(peak, Math.abs(value)), -Infinity);
const pointDeparture = (a, b) => {
  require(Array.isArray(a) && Array.isArray(b) && a.length === b.length, 'Incomplete checkpoint node groups.');
  let maximum = 0;
  for (let g = 0; g < a.length; g++) {
    require(Array.isArray(a[g]) && Array.isArray(b[g]) && a[g].length === b[g].length, 'Incomplete checkpoint node stations.');
    for (let i = 0; i < a[g].length; i++) {
      require(Array.isArray(a[g][i]) && Array.isArray(b[g][i]) && a[g][i].length === b[g][i].length, 'Incomplete checkpoint node banks.');
      for (let j = 0; j < a[g][i].length; j++) {
        const p = a[g][i][j], q = b[g][i][j];
        require(p && q && [p.x, p.y, q.x, q.y].every(Number.isFinite), 'Nonfinite checkpoint node coordinate.');
        maximum = Math.max(maximum, Math.hypot(p.x - q.x, p.y - q.y));
      }
    }
  }
  return maximum;
};

export function initializeCoupledStreamtubeFromFlow(targetMach, checkpoint, { tolerance = 1e-10, blPredictor, targetAlpha } = {}) {
  require(Number.isFinite(targetMach) && targetMach > 0 && targetMach < 1
    && Number.isFinite(tolerance) && tolerance > 0, 'Invalid coupled Mach-restart target or tolerance.');
  require(checkpoint?.version === 1 && checkpoint.restart && checkpoint.continuation && checkpoint.families,
    'Coupled Mach restart requires one complete version-1 ISES checkpoint.');
  require(targetAlpha === undefined || Number.isFinite(targetAlpha), 'Invalid coupled alpha target.');
  const saved = structuredClone(checkpoint), { input, options, initialEuler, initialBL } = saved.restart;
  require(blPredictor === undefined || blPredictor === 'preserve' || blPredictor === 'xfoil-mrchdu', 'Unknown coupled Mach BL predictor.');
  if (blPredictor === 'xfoil-mrchdu') require(options?.transitionMode === 'automatic'
    && Array.isArray(options.tripFractions) && options.tripFractions.every(p => Array.isArray(p) && p.length === 2 && p.every(x => x === 1)),
  'MRCHDU Mach prediction requires automatic transition with terminal material trips.');
  require(input && options && initialEuler && vector(initialEuler.x) && initialEuler.x.length > 0
    && Array.from(initialEuler.x).every(Number.isFinite) && vector(initialBL) && initialBL.length > 0
    && Array.from(initialBL).every(Number.isFinite) && Array.isArray(initialEuler.nodes)
    && Array.isArray(initialEuler.undisplacedNodes), 'Coupled Mach restart requires complete finite Euler/BL and actual/undisplaced node arrays.');
  require(input.streamwiseMode === 'hybrid' && input.upwind
    && (input.flowModel === undefined || input.flowModel === 'compressible')
    && options.blThermodynamics === 'historical-common-isentrope' && options.edgeMatching === 'section-velocity',
  'Coupled Mach restart requires explicit compressible hybrid Euler and historical common-isentrope BL with arithmetic section-velocity matching.');
  const allowed = ['reynolds', 'ncrit', 'tripFractions', 'edgeMatching', 'transitionMode', 'transitionState', 'blThermodynamics', 'hkFloorLinearization', 'geometryReplay'];
  require(Object.keys(options).every(k => allowed.includes(k)), 'Unsupported coupled checkpoint options.');
  const familyNames = ['euler', 'boundaryLayer', 'edgeMatching'];
  require(familyNames.every(k => Number.isFinite(saved.families[k]) && saved.families[k] >= 0
    && saved.families[k] <= tolerance), 'Coupled source checkpoint must already satisfy the residual tolerance.');
  const c = saved.continuation, bodies = input.bodies;
  require(Array.isArray(bodies) && ['convex', 'ises-sampled'].includes(c.iterationGeometry)
    && ['listing', 'admissible', 'armijo', 'event-armijo'].includes(c.stepAcceptance) && ['listing', 'prose'].includes(c.stagnationLimiter)
    && ['amd', 'colamd'].includes(c.preferredOrdering) && [.001, 1].includes(c.pivotTolerance)
    && Array.isArray(c.lastRedistributedStagnation) && c.lastRedistributedStagnation.length === bodies.length
    && c.lastRedistributedStagnation.every(Number.isFinite), 'Invalid coupled checkpoint maintenance history.');
  require(input.stagnationMotion === 'walls-only' && input.normalStencil !== undefined
    && input.geometryDomain === (c.iterationGeometry === 'convex' ? 'convex' : 'positive-simple'),
  'Coupled checkpoint is missing its normalized ISES geometry controls.');
  require(Array.isArray(c.fractions) && c.fractions.length === bodies.length && c.fractions.every((f, b) =>
    Array.isArray(f) && f.length === bodies[b].leadingIndex + 1 && f[0] === 0 && f.at(-1) === 1
    && f.every((v, i) => Number.isFinite(v) && (!i || v > f[i - 1]))), 'Invalid coupled checkpoint inlet fractions.');

  const reconstruct = caseInput => createCoupledStreamtubeBody(caseInput, { ...options,
    initialEuler: { ...initialEuler, x: Float64Array.from(initialEuler.x) }, initialBL: Float64Array.from(initialBL) });
  const source = reconstruct(input), sourceState = source.initial;
  require(initialEuler.x.length === source.ne && initialBL.length === 4 * source.bl.stations.length,
    'Coupled checkpoint encoded state lengths do not match the complete system.');
  require(streamtubeMeshSnapshot({ system: source.euler, nodes: initialEuler.nodes }).quality.valid,
    'Coupled source checkpoint requires a valid convex physical grid.');
  const before = source.evaluate(sourceState);
  require(source.admissible(sourceState), 'Coupled source checkpoint fails physical or convex-grid admissibility.');
  require(familyNames.every(k => before.families[k] === saved.families[k]),
    'Coupled source checkpoint residual does not replay exactly.');
  require(maximum(before.residual) <= tolerance, 'Reconstructed coupled source must satisfy the residual tolerance.');
  const sourceNodeDeparture = pointDeparture(before.outer.nodes, initialEuler.nodes);
  // Subtracting and restoring a displacement can round a reconstructed
  // coordinate by an ulp. Residual families above still replay exactly;
  // reject actual geometry edits, but do not reject arithmetic roundoff.
  const sourceNodesRoundoffEquivalent = before.outer.nodes.every((group, g) => group.every((row, i) => row.every((p, j) =>
    ['x', 'y'].every(key => Math.abs(p[key] - initialEuler.nodes[g][i][j][key])
      <= 8 * Number.EPSILON * Math.max(1, Math.abs(p[key]), Math.abs(initialEuler.nodes[g][i][j][key]))))));
  require(sourceNodesRoundoffEquivalent,
    'Coupled source checkpoint does not reconstruct its actual physical nodes within roundoff.');
  const undisplacedReconstructionDeparture = pointDeparture(before.outer.undisplacedNodes, initialEuler.undisplacedNodes);
  require(same(sourceState.subarray(source.ne), initialBL), 'Coupled source BL unknowns changed during reconstruction.');
  if (options.transitionMode === 'automatic') require(same(source.bl.snapshotActive(), options.transitionState),
    'Coupled source transition map does not replay exactly.');

  // Keep exact same-Mach replay and the original omitted-option arithmetic.
  // A changed-Mach mixed predictor is a new initial guess, not an assertion
  // that the old physical thickness or displaced geometry was preserved.
  if (blPredictor === 'xfoil-mrchdu' && targetAlpha !== undefined && targetAlpha !== (input.alpha ?? 0))
    throw new Error('Simultaneous Mach predictor and alpha transfer are unsupported.');
  if (blPredictor === 'xfoil-mrchdu' && targetMach !== source.euler.conditions.mach)
    return predictedTarget(targetMach,saved,source,sourceState,before,undisplacedReconstructionDeparture,tolerance);

  // Retain the original checkpoint seed, including undisplaced coordinates.
  // Rewriting the latter with a decoded copy can add one-ULP stagnation
  // roundoff and change the exact-resume convention. A new Mach can move
  // natural transition outside the old active interval. Prepare that new
  // phase with the existing criterion before evaluating its interval rows;
  // retain every physical BL variable, rather than reinitializing a profile.
  const targetInput = { ...structuredClone(input), mach: targetMach, ...(targetAlpha === undefined ? {} : { alpha: targetAlpha }) };
  const alphaChanged = targetAlpha !== undefined && targetAlpha !== (input.alpha ?? 0);
  const system = reconstruct(targetInput), initial = system.initial.slice();
  let targetPhaseInitialization;
  if ((alphaChanged || targetMach !== source.euler.conditions.mach) && options.transitionMode === 'automatic') {
    const targets = system.bl.activeTargets(initial.subarray(0, system.ne), initial.subarray(system.ne));
    // An unchanged index can still have an unusable packed-N mixed
    // interval. Honor the selector's existing auxiliary-only repair request
    // without changing physical profiles or reconciling valid N prefixes.
    const reconcileAmplificationSurfaces = targets.flatMap((target, k) => target.amplificationReconciliation ? [k] : []);
    if (targets.some(t => t.from !== t.to) || reconcileAmplificationSurfaces.length) {
      const previous = system.bl.snapshotActive();
      const transferred = system.bl.updateActive(initial.subarray(system.ne), initial.subarray(0, system.ne),
        { reconcileAmplificationSurfaces });
      const auxiliaryChanges = system.bl.stations.flatMap(({ id }) => {
        const before = sourceState[source.ne + 4 * id], after = initial[system.ne + 4 * id];
        return before === after ? [] : [{ id, before, after }];
      });
      targetPhaseInitialization = { method: 'existing-automatic-transition-transfer',
        before: previous, after: system.bl.snapshotActive(), changes: transferred.changes, auxiliaryChanges,
        physicalBLPreserved: true, materialTripsPreserved: true,
        interpretation: 'Target-Mach phase initialization only; no MRCHUE/MRCHDU, physical profile update, or Newton event-step limit.' };
      system.initial.set(initial);
    }
  }
  let value;
  try { value = system.evaluate(initial); }
  catch (error) {
    if (targetPhaseInitialization) error.coupledMachPhaseInitialization = structuredClone(targetPhaseInitialization);
    throw error;
  }
  require(system.admissible(initial), 'Transferred coupled target fails physical or convex-grid admissibility.');
  const physicalStatePreserved = same(initial.subarray(0, system.ne), sourceState.subarray(0, source.ne))
    && Array.from(initial.subarray(system.ne)).every((v, k) => k % 4 === 0 && targetPhaseInitialization || v === sourceState[source.ne + k]);
  require(system.n === source.n && system.ne === source.ne && physicalStatePreserved,
    'Coupled Mach transfer changed a stored density, global, geometry or physical BL coordinate.');
  require(pointDeparture(before.outer.nodes, value.outer.nodes) === 0
    && pointDeparture(before.outer.undisplacedNodes, value.outer.undisplacedNodes) === 0,
  'Coupled Mach transfer changed actual or undisplaced physical geometry.');
  for (const key of (alphaChanged ? ['stagnation', 'strengths'] : ['captured', 'stagnation', 'strengths'])) require(same(before.outer[key], value.outer[key]),
    `Coupled Mach transfer changed physical ${key}.`);
  require((targetPhaseInitialization || same(source.bl.snapshotActive(), system.bl.snapshotActive()))
    && same(source.bl.thicknesses(sourceState.subarray(source.ne)), system.bl.thicknesses(initial.subarray(system.ne))),
  'Coupled Mach transfer changed BL thickness or an unprepared transition mapping.');
  for (let i = 0; i < source.euler.layout.nx; i++) for (let g = 0; g < source.euler.layout.tubes.length; g++)
    for (let j = 0; j < source.euler.layout.tubes[g]; j++) require(before.outer.sections[i][g][j].rho === value.outer.sections[i][g][j].rho,
      'Coupled Mach transfer changed a physical section density.');
  if (!alphaChanged) for (let g = 0; g < source.euler.layout.tubes.length; g++) for (let j = 0; j < source.euler.layout.tubes[g]; j++)
    require(before.outer.allocation.groups[g][j].massFlow === value.outer.allocation.groups[g][j].massFlow,
      'Coupled Mach transfer changed a physical streamtube mass.');
  if (!alphaChanged && targetMach === source.euler.conditions.mach) require(same(value.residual, before.residual),
    'Same-Mach coupled checkpoint residual changed.');
  const targetCheckpoint = structuredClone({ ...saved, families: { ...value.families }, restart: { ...saved.restart, input: targetInput,
    ...(targetPhaseInitialization ? { options: { ...saved.restart.options, transitionState: system.bl.snapshotActive() },
      initialBL: Array.isArray(initialBL) ? Array.from(initial.subarray(system.ne)) : initial.subarray(system.ne).slice() } : {}) } });
  return { system, initial, value, checkpoint: targetCheckpoint, diagnostics: {
    sourceMach: source.euler.conditions.mach, targetMach, sameMach: targetMach === source.euler.conditions.mach,
    sourceFamilies: { ...before.families }, targetFamilies: { ...value.families }, targetResidual: maximum(value.residual),
    sourceConverged: true, targetConverged: maximum(value.residual) <= tolerance,
    physicalDensityPreserved: true, physicalMassPreserved: !alphaChanged,
    ...(alphaChanged ? { sourceAlpha: input.alpha ?? 0, targetAlpha,
      inletMassTreatment: 'recomputed from freestream normal flux at target incidence' } : {}), packedBLPreserved: same(initial.subarray(system.ne), initialBL),
    transitionMapPreserved: same(source.bl.snapshotActive(), system.bl.snapshotActive()),
    ...(targetPhaseInitialization ? { targetPhaseInitialization } : {}),
    physicalGeometryPreserved: true, maintenanceHistoryPreserved: true,
    ...(sourceNodeDeparture > 0 ? { sourceNodeReconstructionDeparture: sourceNodeDeparture } : {}),
    sourceUndisplacedReconstructionDeparture: undisplacedReconstructionDeparture,
    sourceH0: source.euler.conditions.h0, targetH0: system.euler.conditions.h0,
    sourceRhoTotal: source.euler.conditions.rhoTotal, targetRhoTotal: system.euler.conditions.rhoTotal,
    entropyInterpretation: 'Density and mass are retained; relative entropy changes with the target total enthalpy and gas reference.',
    operations: { newtonUpdates: 0, linearSolves: 0, redistributions: 0, coldInitializations: 0, boundaryLayerInitializations: 0 },
    physicalAcceptance: false } };
}

function predictedTarget(targetMach,saved,source,sourceState,before,undisplacedReconstructionDeparture,tolerance) {
  const f=saved.restart,targetInput={...structuredClone(f.input),mach:targetMach};
  const prediction=prepareCoupledMrchduProfiles({bl:source.bl,states:before.layers.states,initialBL:f.initialBL,targetMach});
  const options={...f.options,transitionState:prediction.transitionState};
  let system=createCoupledStreamtubeBody(targetInput,{...options,initialEuler:f.initialEuler,initialBL:prediction.initialBL});
  let initial=system.initial.slice();
  // The native march uses its original local stopping tolerance. Reconcile
  // only auxiliaries, if necessary, against the tightly resolved global
  // transition criterion on the newly PREPARED physical profile.
  const phase=system.bl.updateActive(initial.subarray(system.ne),initial.subarray(0,system.ne));
  options.transitionState=system.bl.snapshotActive();
  // MRCHDU changes physical displacement. The old solved interior must
  // follow the boundary increment before evaluating the new initial guess:
  // holding it fixed can let a surface overtake its first streamline.
  // Wake banks remain independent Euler unknowns. The simultaneous residual
  // still applies wall displacement only; this extension runs only here.
  system.euler.setDisplacement(system.bl.thicknesses(initial.subarray(system.ne)));
  const unextended=system.euler.decode(initial.subarray(0,system.ne));
  const extended=extendWarmBoundaryIncrements({sourceNodes:before.outer.nodes,targetNodes:unextended.nodes,
    masses:before.outer.allocation.groups.map(group=>group.map(tube=>tube.massFlow))});
  const interiorGeometryDeparture=pointDeparture(unextended.nodes,extended);
  initial.set(system.euler.adoptGeometry(initial.subarray(0,system.ne),extended));
  system.initial.set(initial);
  let value=system.evaluate(initial);
  require(system.admissible(initial),'MRCHDU target guess fails physical or convex-grid admissibility.');
  const seed={x:initial.slice(0,system.ne),nodes:value.outer.nodes,undisplacedNodes:value.outer.undisplacedNodes};
  const packed=initial.slice(system.ne);
  // Capture the constructor's actual chart, then establish the canonical
  // checkpoint once. Its undisplaced seed is kept unchanged thereafter:
  // repeatedly replacing it with decoded copies can introduce ULP drift.
  system=createCoupledStreamtubeBody(targetInput,{...options,initialEuler:seed,initialBL:packed});
  initial=system.initial.slice();value=system.evaluate(initial);
  require(system.admissible(initial),'Canonical MRCHDU target checkpoint fails physical or convex-grid admissibility.');
  require(value.residual.every(Number.isFinite),'MRCHDU target residual contains nonfinite rows.');
  require(system.n===source.n && system.ne===source.ne && same(initial.subarray(system.ne),packed),
    'MRCHDU target changed state dimensions or its prepared packed BL fields.');
  for(const key of ['captured','stagnation','strengths']) require(same(before.outer[key],value.outer[key]),
    `MRCHDU target changed the source physical ${key}.`);
  for(let i=0;i<source.euler.layout.nx;i++)for(let g=0;g<source.euler.layout.tubes.length;g++)
    for(let j=0;j<source.euler.layout.tubes[g];j++)require(before.outer.sections[i][g][j].rho===value.outer.sections[i][g][j].rho,
      'MRCHDU target changed a physical Euler density.');
  for(let g=0;g<source.euler.layout.tubes.length;g++)for(let j=0;j<source.euler.layout.tubes[g];j++)
    require(before.outer.allocation.groups[g][j].massFlow===value.outer.allocation.groups[g][j].massFlow,
      'MRCHDU target changed a physical streamtube mass.');
  const physicalGeometryDeparture=pointDeparture(before.outer.nodes,value.outer.nodes);
  const physicalBLPreserved=initial.subarray(system.ne).every((x,k)=>k%4===0||x===f.initialBL[k]);
  const initialEuler={x:seed.x,nodes:value.outer.nodes,undisplacedNodes:seed.undisplacedNodes};
  const initialBL=Array.isArray(f.initialBL)?Array.from(packed):packed;
  const targetCheckpoint=structuredClone({...saved,families:{...value.families},restart:{...f,input:targetInput,options,initialEuler,initialBL}});
  // A serialized copy must reproduce ALL equations exactly, not merely
  // appear positive. This is a frozen replay, without Newton or new meshing.
  const check=createCoupledStreamtubeBody(targetCheckpoint.restart.input,{...targetCheckpoint.restart.options,
    initialEuler:targetCheckpoint.restart.initialEuler,initialBL:targetCheckpoint.restart.initialBL});
  const replay=check.evaluate(check.initial);
  require(same(replay.residual,value.residual)&&same(check.initial,initial)
    &&pointDeparture(replay.outer.nodes,value.outer.nodes)===0,'MRCHDU target checkpoint does not replay exactly.');
  return {system,initial,value,checkpoint:targetCheckpoint,diagnostics:{
    sourceMach:source.euler.conditions.mach,targetMach,sameMach:false,
    sourceFamilies:{...before.families},targetFamilies:{...value.families},targetResidual:maximum(value.residual),
    sourceConverged:true,targetConverged:maximum(value.residual)<=tolerance,
    physicalDensityPreserved:true,physicalMassPreserved:true,packedBLPreserved:same(initial.subarray(system.ne),f.initialBL),
    physicalBLPreserved,transitionMapPreserved:same(source.bl.snapshotActive(),system.bl.snapshotActive()),
    targetBLPrediction:{...prediction.diagnostics,targetTransitionMap:system.bl.snapshotActive(),
      postMarchPhaseChanges:phase.changes,geometryRebuilt:true,physicalGeometryDeparture,checkpointReplayExact:true,
      interiorGridExtension:{method:'boundary-increments-in-physical-mass-coordinates',interiorGeometryDeparture}},
    physicalGeometryPreserved:physicalGeometryDeparture===0,maintenanceHistoryPreserved:true,
    sourceUndisplacedReconstructionDeparture:undisplacedReconstructionDeparture,
    sourceH0:source.euler.conditions.h0,targetH0:system.euler.conditions.h0,
    sourceRhoTotal:source.euler.conditions.rhoTotal,targetRhoTotal:system.euler.conditions.rhoTotal,
    entropyInterpretation:'Euler density and mass are retained; relative entropy changes with target enthalpy, reference gas and the displaced geometry. No isentropic density inversion is used.',
    operations:{newtonUpdates:0,linearSolves:0,redistributions:0,coldInitializations:0,boundaryLayerInitializations:1,
      translatedMRCHDUCalls:prediction.diagnostics.operations.translatedMRCHDUCalls,interiorGridExtensions:1,checkpointReplays:1},physicalAcceptance:false}};
}
