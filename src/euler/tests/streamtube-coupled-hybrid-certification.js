// SPDX-License-Identifier: GPL-2.0-or-later
// Explicit same-Mach model certification, not a generic restart or initializer.
// The research controls must be supplied. A single nonzero floating-point
// transport-speed change or changed residual row rejects the conversion.
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { initializeCoupledStreamtubeFromFlow } from './streamtube-coupled-flow-restart.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';

const require = (condition, message) => { if (!condition) throw new Error(message); };
const vector = v => Array.isArray(v) || ArrayBuffer.isView(v) && !(v instanceof DataView);
const same = (a, b) => {
  if (a === b) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  if (vector(a) || vector(b)) return vector(a) && vector(b) && a.length === b.length && Array.from(a).every((v, i) => same(v, b[i]));
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};
const maximum = x => x.reduce((peak, value) => Math.max(peak, Math.abs(value)), -Infinity);

export function certifyCoupledStreamtubeHybrid(checkpoint, { epsilonP, upwind, tolerance = 1e-10 } = {}) {
  require(Number.isFinite(epsilonP) && epsilonP > 0 && Number.isFinite(tolerance) && tolerance > 0,
    'Hybrid certification requires explicit positive epsilonP and tolerance.');
  require(upwind && typeof upwind === 'object' && !Array.isArray(upwind)
    && Object.keys(upwind).every(k => ['mucon', 'mcrit', 'boundary'].includes(k))
    && Number.isFinite(upwind.mucon) && Number.isFinite(upwind.mcrit) && upwind.mcrit >= 0 && upwind.mcrit <= 1
    && upwind.boundary?.kind === 'unfiltered-first-two' && Object.keys(upwind.boundary).length === 1,
  'Hybrid certification requires explicit normalized mucon, mcrit and unfiltered-first-two controls.');
  require(checkpoint?.version === 1 && checkpoint.restart && checkpoint.continuation && checkpoint.families,
    'Hybrid certification requires a complete version-1 coupled ISES checkpoint.');
  const saved = structuredClone(checkpoint), controls = structuredClone(upwind);
  const { input, options, initialEuler, initialBL } = saved.restart;
  require(input && options && initialEuler && vector(initialEuler.x) && initialEuler.x.length > 0
    && Array.from(initialEuler.x).every(Number.isFinite) && vector(initialBL) && initialBL.length > 0
    && Array.from(initialBL).every(Number.isFinite) && Array.isArray(initialEuler.nodes)
    && Array.isArray(initialEuler.undisplacedNodes), 'Hybrid certification requires complete finite supplied Euler/BL and actual/undisplaced nodes.');
  require(input.streamwiseMode === 'isentropic' && input.upwind === undefined && input.hybrid === undefined
    && (input.flowModel === undefined || input.flowModel === 'compressible')
    && options.edgeMatching === 'section-velocity' && options.blThermodynamics === undefined,
  'Hybrid certification requires baseline compressible isentropic Euler and arithmetic section-velocity BL, without existing hybrid/upwind/thermodynamic overrides.');
  const allowed = ['reynolds', 'ncrit', 'tripFractions', 'edgeMatching', 'transitionMode', 'transitionState', 'hkFloorLinearization', 'geometryReplay'];
  require(Object.keys(options).every(k => allowed.includes(k)), 'Unsupported source checkpoint options for hybrid certification.');
  const names = ['euler', 'boundaryLayer', 'edgeMatching'];
  require(Object.keys(saved.families).length === names.length && names.every(k => Number.isFinite(saved.families[k])
    && saved.families[k] >= 0 && saved.families[k] <= tolerance), 'Hybrid certification source must already satisfy the residual tolerance.');
  const c = saved.continuation;
  require(input.stagnationMotion === 'walls-only' && input.normalStencil !== undefined
    && ['convex', 'ises-sampled'].includes(c.iterationGeometry)
    && input.geometryDomain === (c.iterationGeometry === 'convex' ? 'convex' : 'positive-simple')
    && Number.isFinite(input.mach) && input.mach > 0 && input.mach < 1,
  'Hybrid certification requires normalized ISES geometry and Mach controls.');
  const construct = (caseInput, caseOptions) => createCoupledStreamtubeBody(caseInput, { ...caseOptions,
    initialEuler: { ...initialEuler, x: Float64Array.from(initialEuler.x) }, initialBL: Float64Array.from(initialBL) });
  let source = construct(input, options), before = source.evaluate(source.initial);
  require(source.ne === initialEuler.x.length && initialBL.length === 4 * source.bl.stations.length,
    'Hybrid certification source state lengths are incomplete.');
  require(source.initialization.suppliedBL && source.admissible(source.initial)
    && streamtubeMeshSnapshot({ system: source.euler, nodes: initialEuler.nodes }).quality.valid,
  'Hybrid certification source must have admissible BL and valid convex physical geometry.');
  require(names.every(k => before.families[k] === saved.families[k]) && maximum(before.residual) <= tolerance,
    'Hybrid certification source residual does not replay exactly or is not converged.');
  require(same(before.outer.nodes, initialEuler.nodes) && same(source.initial.subarray(source.ne), initialBL),
    'Hybrid certification source physical nodes or packed BL do not replay exactly.');
  if (options.transitionMode === 'automatic') require(same(source.bl.snapshotActive(), options.transitionState),
    'Hybrid certification source transition map does not replay exactly.');
  const gas = source.euler.conditions;
  const subsonic = ue => Number.isFinite(ue) && ue > 0
    && ue * ue < (gas.gamma - 1) * (gas.h0 - .5 * ue * ue);
  require(before.layers.states.every(s => subsonic(s.ue))
    && before.edges.every(e => subsonic(e.ue) && (!e.sides || e.sides.every(s => subsonic(s.ue)))),
  'Hybrid certification cannot reinterpret sonic or supersonic source BL/bank states.');
  let sectionCount = 0, maximumMach = 0;
  for (const row of before.outer.sections) for (const group of row) for (const section of group) {
    require(Number.isFinite(section.machSquared) && section.machSquared < 1,
      'Hybrid certification cannot reinterpret a sonic or supersonic Euler source.');
    sectionCount++; maximumMach = Math.max(maximumMach, Math.sqrt(section.machSquared));
  }
  const converted = structuredClone(saved);
  converted.restart.input = { ...converted.restart.input, streamwiseMode: 'hybrid', hybrid: { epsilonP }, upwind: controls };
  converted.restart.options = { ...converted.restart.options, blThermodynamics: 'historical-common-isentrope' };
  let target = construct(converted.restart.input, converted.restart.options), after = target.evaluate(target.initial);
  require(target.initialization.suppliedBL && target.admissible(target.initial), 'Converted hybrid state is not admissible.');
  require(source.n === target.n && source.ne === target.ne && same(source.initial, target.initial),
    'Hybrid certification changed a stored Euler or packed BL coordinate.');
  for (let i = 0; i < source.euler.layout.nx; i++) for (let g = 0; g < source.euler.layout.tubes.length; g++)
    for (let j = 0; j < source.euler.layout.tubes[g]; j++) {
      const a = before.outer.sections[i][g][j], b = after.outer.sections[i][g][j];
      if (after.outer.transportSpeeds[i][g][j] !== b.q) throw Object.assign(
        new Error(`Hybrid certification rejects nonzero speed bias at section ${i}, group ${g}, tube ${j}.`),
        { code: 'coupled-hybrid-certification-speed-bias', diagnostics: { section: i, group: g, tube: j,
          mach: gas.mach, physicalSpeed: b.q, transportSpeed: after.outer.transportSpeeds[i][g][j],
          bias: after.outer.transportSpeeds[i][g][j] - b.q } });
      require(['rho', 'q', 'p', 'enthalpy', 'machSquared'].every(k => a[k] === b[k]),
        'Hybrid certification changed physical Euler gas.');
    }
  require(same(before.residual, after.residual), 'Hybrid certification requires every complete residual row to be exactly equal.');
  require(same(before.outer.nodes, after.outer.nodes) && same(before.outer.undisplacedNodes, after.outer.undisplacedNodes),
    'Hybrid certification changed actual or undisplaced geometry.');
  for (const key of ['captured', 'stagnation', 'strengths']) require(same(before.outer[key], after.outer[key]),
    `Hybrid certification changed physical ${key}.`);
  for (let g = 0; g < source.euler.layout.tubes.length; g++) for (let j = 0; j < source.euler.layout.tubes[g]; j++)
    require(before.outer.allocation.groups[g][j].massFlow === after.outer.allocation.groups[g][j].massFlow,
      'Hybrid certification changed physical streamtube mass.');
  require(same(source.bl.snapshotActive(), target.bl.snapshotActive())
    && same(source.bl.thicknesses(source.initial.subarray(source.ne)), target.bl.thicknesses(target.initial.subarray(target.ne)))
    && same(before.layers.states, after.layers.states), 'Hybrid certification changed BL, wake or transition state.');
  const residual = after.residual.slice(), diagnostics = { certified: true,
    kind: 'exact-same-mach-isentropic-to-hybrid', mach: gas.mach, epsilonP, upwind: structuredClone(controls),
    sourceFamilies: { ...before.families }, targetFamilies: { ...after.families },
    residualRowsExactlyEqual: source.n, physicalSectionsExactlyEqual: sectionCount,
    sourceMaximumMach: maximumMach, maximumSpeedBias: 0,
    boundaryLayers: source.bl.surfaces.length, wakes: source.bl.wakes.length,
    packedBLPreserved: true, physicalDensityPreserved: true, physicalMassPreserved: true,
    actualAndUndisplacedGeometryPreserved: true, transitionMapPreserved: true, maintenanceHistoryPreserved: true,
    physicalAcceptance: false,
    interpretation: 'Explicit numerical same-state certification only. The historical common-isentrope BL remains a stated transonic approximation; no postshock physical validation is implied.' };
  // Release the first two complete numerical reconstructions before the
  // existing strict helper performs its independent structural replay.
  source = null; target = null; before = null; after = null;
  const verified = initializeCoupledStreamtubeFromFlow(gas.mach, converted, { tolerance });
  require(same(verified.value.residual, residual) && same(verified.checkpoint, converted),
    'Certified hybrid checkpoint failed the strict complete-flow replay.');
  return { checkpoint: structuredClone(verified.checkpoint), diagnostics: {
    ...diagnostics, strictCompleteFlowReplay: true,
    sourceUndisplacedReconstructionDeparture: verified.diagnostics.sourceUndisplacedReconstructionDeparture,
    operations: { sourceTargetReconstructions: 2, strictCompleteFlowHelperCalls: 1,
      newtonUpdates: 0, linearSolves: 0, redistributions: 0, coldInitializations: 0, boundaryLayerInitializations: 0 } } };
}
