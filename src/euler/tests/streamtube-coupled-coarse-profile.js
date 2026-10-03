// SPDX-License-Identifier: GPL-2.0-or-later
// A BL-only initial guess on a separately constructed Euler grid, not a
// nested grid transfer or a converged fine-grid solution.
import { createStreamtubeBodySystem } from '../streamtube-body.js';
import { createStreamtubeBoundaryLayers } from '../streamtube-boundary-layers.js';
import { initializeCoupledStreamtubeBody } from '../streamtube-coupled-initializer.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';
import { initialStreamtubeDisplacement } from '../streamtube-geometry.js';
import { createContourTopology } from '../../geometry/contour-topology.js';

const require = (ok, message) => { if (!ok) throw new Error(message); };
const keys = ['lengthScale', 'massScale', 'center', 'mach', 'alpha', 'gamma', 'flowModel', 'streamwiseMode',
  'reynolds', 'ncrit', 'edgeMatching', 'transitionMode', 'blThermodynamics', 'pressureCorrectionFactor',
  'normalStencil', 'stagnationMotion', 'wakeGeometry', 'wakeOutlet', 'wakeDisplacementMotion', 'upwind', 'hybrid'];
const conditionValue = (conditions, key) => key === 'wakeDisplacementMotion' ? conditions[key] ?? 'fixed' : conditions[key];
const same = (a, b) => {
  if (a === b) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  if (Array.isArray(a) || ArrayBuffer.isView(a)) return (Array.isArray(b) || ArrayBuffer.isView(b))
    && a.length === b.length && Array.from(a).every((v, i) => same(v, b[i]));
  const names = Object.keys(a);
  return names.length === Object.keys(b).length && names.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};
const ordered = a => a.length >= 2 && a.every((v, i) => Number.isFinite(v) && (!i || v > a[i - 1]));
const mix = (a, b, t) => t === 0 ? a : t === 1 ? b : a + t * (b - a);
const interpolate = (a, b, t) => Object.fromEntries(['theta', 'deltaStar', 'ue', 'aux'].map(k => [k, mix(a[k], b[k], t)]));
function bracket(coordinates, x) {
  require(Number.isFinite(x) && x >= coordinates[0] && x <= coordinates.at(-1), 'BL transfer coordinate is outside its source branch.');
  if (x === coordinates.at(-1)) return { left: coordinates.length - 2, right: coordinates.length - 1, t: 1 };
  let left = 0, right = coordinates.length - 1;
  while (right - left > 1) { const mid = (left + right) >> 1; if (coordinates[mid] > x) right = mid; else left = mid; }
  return { left, right, t: (x - coordinates[left]) / (coordinates[right] - coordinates[left]) };
}
function validateDescriptor(d, source) {
  require(d?.conditions?.transitionMode === 'automatic' && Array.isArray(d.bodies) && d.bodies.length > 0
    && Array.isArray(d.surfaces) && d.surfaces.length === 2 * d.bodies.length
    && Array.isArray(d.wakes) && d.wakes.length === d.bodies.length, 'BL transfer requires complete automatic body profiles.');
  require(keys.filter(k => !['blThermodynamics', 'upwind', 'hybrid', 'wakeDisplacementMotion'].includes(k)).every(k => d.conditions[k] !== undefined)
    && ['fixed', 'te-center'].includes(conditionValue(d.conditions, 'wakeDisplacementMotion'))
    && ['lengthScale', 'massScale', 'mach', 'alpha', 'gamma', 'reynolds', 'ncrit', 'pressureCorrectionFactor']
      .every(k => Number.isFinite(d.conditions[k]))
    && Number.isFinite(d.conditions.center?.x) && Number.isFinite(d.conditions.center?.y)
    && d.conditions.massScale > 0 && d.conditions.gamma > 1 && d.conditions.ncrit > 0,
  'BL transfer requires complete finite physical/model conditions.');
  require(Number.isFinite(d.scale) && d.scale > 0 && d.conditions.lengthScale > 0
    && Number.isFinite(d.conditions.lengthScale) && d.conditions.reynolds > 0
    && Number.isFinite(d.conditions.reynolds) && d.scale === 1 / Math.sqrt(d.conditions.reynolds), 'Invalid BL transfer normalization.');
  require(Array.isArray(d.trips) && d.trips.length === d.bodies.length
    && d.trips.every(p => Array.isArray(p) && p.length === 2 && p.every(x => x === 1)), 'BL transfer currently requires terminal material trips.');
  const elements = new Set();
  for (const body of d.bodies) {
    require(Number.isInteger(body.element) && body.element >= 0 && !elements.has(body.element), 'BL transfer requires unique explicit element identities.');
    elements.add(body.element);
    require((body.trailingEdge === undefined || ['sharp', 'finite-base'].includes(body.trailingEdge?.kind)) && Array.isArray(body.points) && body.points.length >= 4
      && body.points.every(p => p && Number.isFinite(p.x) && Number.isFinite(p.y))
      && body.points[0].x === body.points.at(-1).x && body.points[0].y === body.points.at(-1).y,
    'BL transfer requires a closed solid contour with explicit sharp or finite-base topology.');
    if (body.trailingEdge?.kind === 'finite-base') {
      createContourTopology(body.points, body);
      require(d.conditions.wakeGeometry === 'independent-banks', 'Finite-base BL transfer requires independent wake banks.');
    }
  }
  const data = source ? d.states : d.coordinates;
  require(Array.isArray(data) && data.length > 0 && data.every(p => p && Number.isFinite(p.s) && p.s > 0
    && (p.wakeGap === undefined || Number.isFinite(p.wakeGap) && p.wakeGap >= 0)), 'BL transfer requires finite physical station arcs and nonnegative base gaps.');
  if (source) require(data.every(p => [p.theta, p.deltaStar, p.ue, p.aux].every(Number.isFinite)
    && p.theta > 0 && p.deltaStar - (p.wakeGap ?? 0) > p.theta && p.ue > 0 && p.aux >= 0), 'Invalid source BL physical profile.');
  const covered = new Set(), parts = new Set();
  for (const p of [...d.surfaces, ...d.wakes]) {
    require(Number.isInteger(p.body) && p.body >= 0 && p.body < d.bodies.length
      && Array.isArray(p.ids) && p.ids.length >= 2, 'BL transfer has incomplete body station metadata.');
    const surface = d.surfaces.includes(p), label = `${p.body}:${surface ? p.side : 'wake'}`;
    require(!parts.has(label) && (!surface || ['upper', 'lower'].includes(p.side)), 'BL transfer duplicates a body branch.');
    parts.add(label);
    require(p.ids.every(id => {
      if (!Number.isInteger(id) || id < 0 || id >= data.length || covered.has(id)) return false;
      covered.add(id); return true;
    }), 'BL transfer station ids overlap or leave their array.');
    require(ordered(p.ids.map(id => data[id].s)), 'BL transfer station arcs must increase downstream.');
    const finiteWake = !surface && d.bodies[p.body].trailingEdge?.kind === 'finite-base';
    require(p.ids.every(id => finiteWake ? Number.isFinite(data[id].wakeGap)
      : data[id].wakeGap === undefined || data[id].wakeGap === 0),
    finiteWake ? 'Finite-base BL transfer needs every decoded physical wake gap.' : 'Surface and sharp-wake stations require zero base gap.');
    if (finiteWake) require(data[p.ids[0]].wakeGap > 0, 'Finite-base wake start requires its positive solid-base gap.');
    if (source && surface) require(Number.isInteger(p.transition) && p.transition >= 0 && p.transition < p.ids.length
      && Number.isFinite(p.transitionS) && p.transitionS > 0 && p.transitionS <= data[p.ids.at(-1)].s,
    'BL transfer needs the source transition location and phase.');
    if (source) require(p.ids.every((id, i) => surface && i < p.transition
      ? data[id].aux < d.conditions.ncrit : data[id].aux > 0), 'Source BL auxiliary state disagrees with its phase.');
  }
  require(covered.size === data.length && parts.size === 3 * d.bodies.length, 'BL transfer omitted one or more stations or body branches.');
}

// Pure except for the optional pure scalar native closure used when a new
// leading station lies inside an already turbulent leading interval.
export function mapCoupledBLProfiles({ source, target, transitionShear } = {}) {
  validateDescriptor(source, true); validateDescriptor(target, false);
  for (const key of keys) require(same(conditionValue(source.conditions, key), conditionValue(target.conditions, key)),
    `BL transfer changed physical/model condition ${key}.`);
  require(source.bodies.length === target.bodies.length, 'BL transfer body count changed.');
  const bodyMap = target.bodies.map(body => {
    const b = source.bodies.findIndex(s => s.element === body.element);
    require(b >= 0 && same(source.bodies[b].points.map(p => [p.x, p.y]), body.points.map(p => [p.x, p.y])),
      'BL transfer element identity or physical solid contour changed.');
    const previous = source.bodies[b].trailingEdge, next = body.trailingEdge;
    require((previous?.kind ?? 'sharp') === (next?.kind ?? 'sharp')
      && (next?.kind !== 'finite-base' || previous.upperIndex === next.upperIndex && previous.lowerIndex === next.lowerIndex),
    'BL transfer changed material trailing-edge topology.');
    return b;
  });
  const values = new Float64Array(4 * target.coordinates.length), phase = [], mapping = { bodyMap, surfaces: [], wakes: [] };
  const finiteBase = target.bodies.some(b => b.trailingEdge?.kind === 'finite-base'), fluidWakeDisplacement = [];
  const put = (id, p) => values.set([p.aux, p.theta / target.scale, p.deltaStar / target.scale, p.ue], 4 * id);
  for (const surface of target.surfaces) {
    const old = source.surfaces.find(s => s.body === bodyMap[surface.body] && s.side === surface.side);
    const parents = old.ids.map(id => source.states[id]), oldLength = parents.at(-1).s;
    const newLength = target.coordinates[surface.ids.at(-1)].s, coordinates = parents.map(p => p.s / oldLength);
    const targetCoordinates = surface.ids.map(id => target.coordinates[id].s / newLength);
    const location = old.transitionS / oldLength, transitionIndex = targetCoordinates.findIndex(s => s >= location);
    require(transitionIndex >= 0, 'Mapped BL transition left its target branch.'); phase.push(transitionIndex);
    const rows = [];
    for (const [j, id] of surface.ids.entries()) {
      const coordinate = targetCoordinates[j]; let v, details;
      if (coordinate < coordinates[0]) {
        v = { ...parents[0], aux: 0, ue: parents[0].ue * (coordinate * oldLength) / parents[0].s };
        if (j >= transitionIndex) {
          require(typeof transitionShear === 'function', 'Leading turbulent BL transfer requires the native shear closure.');
          v.aux = transitionShear({ ...v, s: target.coordinates[id].s, aux: .03 });
          require(Number.isFinite(v.aux) && v.aux > 0, 'Invalid transferred leading turbulent shear.');
        }
        details = { leadingSimilarity: true };
      } else {
        const { left, right, t } = bracket(coordinates, coordinate); v = interpolate(parents[left], parents[right], t);
        const leftTurbulent = left >= old.transition, rightTurbulent = right >= old.transition;
        if (leftTurbulent !== rightTurbulent) v.aux = parents[leftTurbulent === (j >= transitionIndex) ? left : right].aux;
        details = { left, right, fraction: t, crossedAuxiliaryPhase: leftTurbulent !== rightTurbulent };
      }
      put(id, v); rows.push({ id, s: target.coordinates[id].s, sourceArcFraction: coordinate, ...details });
    }
    mapping.surfaces.push({ body: surface.body, sourceBody: old.body, element: target.bodies[surface.body].element,
      side: surface.side, sourceLength: oldLength, targetLength: newLength, sourceTransition: old.transition,
      mappedTransition: transitionIndex, transitionArcFraction: location, rows });
  }
  for (const wake of target.wakes) {
    const old = source.wakes.find(w => w.body === bodyMap[wake.body]), original = old.ids.map(id => source.states[id]);
    const finiteWake = target.bodies[wake.body].trailingEdge?.kind === 'finite-base';
    // The native unknown contains TOTAL displacement. XICALC's geometric
    // dead-air tail is evaluated on each grid's own physical mean-bank arc;
    // only the remaining fluid profile is an interpolated BL quantity.
    const parents = finiteWake ? original.map(p => ({ ...p, deltaStar: p.deltaStar - p.wakeGap })) : original;
    if (finiteWake) require(original[0].wakeGap === target.coordinates[wake.ids[0]].wakeGap,
      'BL transfer changed the physical trailing-edge base width.');
    const oldStart = parents[0].s, oldLength = parents.at(-1).s - oldStart;
    const newStart = target.coordinates[wake.ids[0]].s, newLength = target.coordinates[wake.ids.at(-1)].s - newStart;
    require(oldLength > 0 && newLength > 0, 'BL transfer wake has no positive length.');
    const coordinates = parents.map(p => (p.s - oldStart) / oldLength), rows = [];
    for (const id of wake.ids) {
      const coordinate = (target.coordinates[id].s - newStart) / newLength, b = bracket(coordinates, coordinate);
      const v = interpolate(parents[b.left], parents[b.right], b.t);
      if (finiteWake) {
        require(Number.isFinite(v.deltaStar) && v.deltaStar > v.theta, 'Transferred fluid wake displacement is outside the model domain.');
        fluidWakeDisplacement.push({ id, deltaStar: v.deltaStar });
        const gap = target.coordinates[id].wakeGap, endpoint = b.t === 0 ? original[b.left] : b.t === 1 ? original[b.right] : null;
        // A retained station with the same gap keeps its original total
        // exactly, including the floating-point representation.
        v.deltaStar = endpoint?.wakeGap === gap ? endpoint.deltaStar : v.deltaStar + gap;
      }
      put(id, v); rows.push({ id, sourceArcFraction: coordinate, ...b });
    }
    mapping.wakes.push({ body: wake.body, sourceBody: old.body, element: target.bodies[wake.body].element, sourceLength: oldLength, targetLength: newLength, rows });
  }
  return { initialBL: values, transitionState: phase, mapping, ...(finiteBase ? { fluidWakeDisplacement } : {}) };
}

const descriptor = (euler, bl, conditions, extra) => ({ conditions: { ...euler.conditions, ...conditions },
  bodies: euler.layout.bodies, scale: bl.scale, trips: bl.trips, surfaces: bl.surfaces, wakes: bl.wakes, ...extra });
const familiesPass = (families, tolerance) => ['euler', 'boundaryLayer', 'edgeMatching'].every(k =>
  Number.isFinite(families?.[k]) && families[k] >= 0 && families[k] <= tolerance);

export function prepareCoupledCoarseProfile({ sourceSystem, sourceResult, input, options, initialEuler } = {},
  { tolerance = 1e-10, onAttempt } = {}) {
  require(Number.isFinite(tolerance) && tolerance > 0 && (onAttempt === undefined || typeof onAttempt === 'function'), 'Invalid coarse BL preparation controls.');
  require(sourceSystem?.bl && sourceResult?.converged === true && sourceResult.mesh?.quality?.valid === true
    && familiesPass(sourceResult.families, tolerance), 'Coarse BL transfer requires a converged accepted source.');
  require(options?.transitionMode === 'automatic' && initialEuler?.x && initialEuler.nodes && initialEuler.undisplacedNodes,
    'Coarse BL preparation requires a complete target Euler state and automatic transition.');
  const before = sourceSystem.evaluate(sourceSystem.initial);
  require(same(before.residual, sourceResult.residual) && same(before.families, sourceResult.families)
    && familiesPass(before.families, tolerance) && sourceSystem.admissible(sourceSystem.initial), 'Accepted coarse BL source does not replay exactly.');
  require(streamtubeMeshSnapshot({ system: sourceSystem.euler, nodes: before.outer.nodes }).quality.valid, 'Replayed coarse source grid is invalid.');
  const nextInput = structuredClone(input), nextOptions = structuredClone(options), targetSeed = structuredClone(initialEuler);
  delete nextOptions.initialEuler; delete nextOptions.initialBL; delete nextOptions.transitionState;
  require(nextInput.displacement === undefined, 'Coarse BL transfer owns its prescribed displacement.');
  const displacement = { surfaces: nextInput.bodies.map(b => ({ upper: Array(b.trailingIndex - b.leadingIndex + 1).fill(0),
    lower: Array(b.trailingIndex - b.leadingIndex + 1).fill(0) })),
    wakes: nextInput.bodies.map(b => Array(nextInput.outerLower.length - 1 - b.trailingIndex).fill(0)) };
  // Only geometry and BL metadata are created here: no MRCHUE, density
  // inversion, missing-profile initialization, or fine Euler solve.
  const euler = createStreamtubeBodySystem({ ...nextInput, displacement });
  if (euler.baseGeometry?.some(Boolean)) euler.setDisplacement(initialStreamtubeDisplacement(euler.layout, euler.baseGeometry));
  const state = euler.adoptGeometry(targetSeed.x, targetSeed.undisplacedNodes);
  const historical = nextOptions.blThermodynamics === 'historical-common-isentrope';
  const bl = createStreamtubeBoundaryLayers(euler, state, { ...nextOptions, ...(historical ? { allowSupersonicEdge: true } : {}) });
  const geo = bl.geometry(state), targetFlow = euler.decode(state);
  const old = descriptor(sourceSystem.euler, sourceSystem.bl, sourceSystem.conditions, { states: before.layers.states,
    surfaces: sourceSystem.bl.surfaces.map(s => ({ ...s,
      transitionS: before.layers.transitions.find(t => t.body === s.body && t.side === s.side)?.s })) });
  const target = descriptor(euler, bl, { ...nextOptions, edgeMatching: nextOptions.edgeMatching ?? 'pressure', reynolds: bl.kernel.parameters.reynolds,
    ncrit: bl.kernel.parameters.ncrit }, { coordinates: geo.coordinates });
  const mapped = mapCoupledBLProfiles({ source: old, target,
    transitionShear: p => bl.kernel.station(p, 'turbulent').transitionShear });
  bl.restoreActive(mapped.transitionState);
  const phasePreparation = bl.updateActive(mapped.initialBL, state, { reinitializeAmplification: true });
  const supplied = { ...nextOptions, initialEuler: targetSeed, initialBL: mapped.initialBL, transitionState: bl.snapshotActive() };
  const prepared = initializeCoupledStreamtubeBody(nextInput, supplied, { initialThicknessFactor: 1, maximumBacktracks: 0, onAttempt });
  const system = prepared.system;
  let finiteBaseTransfer;
  if (mapped.fluidWakeDisplacement) {
    // The supplied-profile initializer may extend the TE center or adjust
    // wake correspondence. Restore the prescribed gap on that FINAL arc.
    // Independent banks are coordinates: changing total wake displacement
    // here cannot move them or change the surface BL/transition states.
    require(system.euler.layout.independentWakeBanks, 'Finite-base profile preparation lost its independent wake banks.');
    const finalGeo = system.bl.geometry(system.initial.subarray(0, system.ne));
    let maximumGapReconstructionChange = 0;
    for (const { id, deltaStar: fluid } of mapped.fluidWakeDisplacement) {
      const gap = finalGeo.coordinates[id].wakeGap, previous = geo.coordinates[id].wakeGap;
      require(Number.isFinite(gap) && gap >= 0, 'Final transferred wake gap is outside the geometric domain.');
      maximumGapReconstructionChange = Math.max(maximumGapReconstructionChange, Math.abs(gap - previous));
      if (gap !== previous) mapped.initialBL[4 * id + 2] = (fluid + gap) / system.bl.scale;
      system.initial[system.ne + 4 * id + 2] = mapped.initialBL[4 * id + 2];
    }
    system.euler.setDisplacement(system.bl.thicknesses(system.initial.subarray(system.ne)));
    const checkedGeo = system.bl.geometry(system.initial.subarray(0, system.ne));
    require(mapped.fluidWakeDisplacement.every(({ id }) => checkedGeo.coordinates[id].wakeGap === finalGeo.coordinates[id].wakeGap
      && checkedGeo.coordinates[id].s === finalGeo.coordinates[id].s), 'Restoring the transferred dead-air gap changed its physical wake arc.');
    finiteBaseTransfer = { model: 'XFOIL dead-air gap on final mean-bank arc from solid TE center',
      transferredWakeStations: mapped.fluidWakeDisplacement.length, maximumGapReconstructionChange,
      interpolation: 'fluid deltaStar = source total deltaStar minus source decoded wakeGap; target total = interpolated fluid plus final target wakeGap' };
  }
  const value = system.evaluate(system.initial);
  require(prepared.mesh.quality.valid && system.admissible(system.initial) && value.residual.every(Number.isFinite), 'Mapped fine BL seed failed its full physical/grid domain.');
  require(prepared.initialization.thicknessFactor === 1 && prepared.initialization.history.length === 1
    && system.initial.subarray(system.ne).every((v, k) => k % 4 === 0 || v === mapped.initialBL[k]), 'Coarse BL preparation changed mapped physical fields or thinned them.');
  if (finiteBaseTransfer) prepared.initialization.history.at(-1).families = { ...value.families };
  require(Array.from(system.initial.subarray(0, system.euler.layout.densityCount)).every((v, i) => v === targetSeed.x[i])
    && value.outer.sections.every((row, i) => row.every((group, g) => group.every((p, j) => p.rho ===
      (euler.layout.densityUnknowns ? Math.exp(targetSeed.x[euler.layout.densityIndex(i, g, j)]) : 1)))), 'Mapped seed changed fine Euler density.');
  for (const key of ['captured', 'stagnation']) require(same(value.outer[key], targetFlow[key]), `Mapped seed changed fine Euler ${key}.`);
  require(Array.from(system.initial.subarray(system.euler.layout.globalOffset, system.ne)).every((v, i) =>
    v === state[euler.layout.globalOffset + i]), 'Mapped seed changed fine Euler global unknowns.');
  require(value.outer.allocation.groups.every((group, g) => group.every((tube, j) =>
    tube.massFlow === targetFlow.allocation.groups[g][j].massFlow)), 'Mapped seed changed fine streamtube mass.');
  const transfer = { method: 'converged-coarse-bl-profile', initialGuessOnly: true, equationsChanged: false,
    sourceFamilies: { ...before.families }, targetFamilies: { ...value.families },
    sourceUnknowns: sourceSystem.n, targetUnknowns: system.n,
    sourceTubes: sourceSystem.euler.layout.tubes.slice(), targetTubes: system.euler.layout.tubes.slice(),
    sourceSegments: sourceSystem.euler.layout.nx, targetSegments: system.euler.layout.nx,
    initializationGeometryDomains: { source: sourceSystem.euler.conditions.geometryDomain, target: euler.conditions.geometryDomain },
    mapping: mapped.mapping, initialPhase: mapped.transitionState, preparedPhase: system.bl.snapshotActive(), phasePreparation,
    ...(finiteBaseTransfer ? { finiteBaseTransfer } : {}),
    physicalBLInterpolation: 'Normalized physical branch arc and TE-to-outlet wake arc; finite leading thickness and linear edge speed; phase-aware auxiliary interpolation followed by native N preparation.',
    targetCoordinateSource: 'Retained fine Euler undisplaced chart; no original failed BL profile is required.',
    maintenanceHistoryInherited: false, thicknessFactor: 1,
    operations: { sourceResidualReplays: 1, suppliedBLInitializations: 1, MRCHUE: 0, MRCHDU: 0, globalNewtonUpdates: 0, globalLinearSolves: 0 } };
  return { ...prepared, value, transfer, initialization: { ...prepared.initialization,
    method: 'converged coarse BL interpolation onto requested Euler topology', coarseProfileTransfer: structuredClone(transfer) } };
}
