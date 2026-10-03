// SPDX-License-Identifier: GPL-2.0-or-later
import { prolongStreamtubeDensities } from './streamtube-density-prolongation.js';
// Nested refinement for coupled flow convergence studies.
// Subdivide surface/cut intervals and selected captured streamtube masses.
// Interpolation supplies a new guess; the full coupled equations must be solved.
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { createStreamtubeBoundaryLayers } from './streamtube-boundary-layers.js';
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { initializeStreamtubeDensities } from './streamtube-initial-state.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { refineStreamtubeMassCoordinates } from '../geometry/streamtube-mass-interpolation.js';
import { nestedIntervalMap } from '../geometry/nested-intervals.js';
import { createStreamtubeRefinementCoordinates } from '../geometry/streamtube-refinement-coordinates.js';

const mix = (a, b, t) => t === 0 ? a : t === 1 ? b : a + t * (b - a);
const point = (a, b, t) => ({ x: mix(a.x, b.x, t), y: mix(a.y, b.y, t) });
const sample = (row, u, interpolate) => {
  const i = Math.min(row.length - 2, Math.floor(u));
  return interpolate(row[i], row[i + 1], u - i);
};
export function refineCoupledStreamtubeBody(input, source, { initial = source.initial, streamwiseFactor, streamwiseSubdivisions,
  normalFactor, normalSubdivisions, maxNodes = 50000, normalInterpolation = 'linear', streamwiseInterpolation = 'linear' } = {}) {
  if (!['linear', 'streamfunction-quadratic'].includes(normalInterpolation)) throw new Error('Unknown coupled refinement interpolation.');
  if (!['linear', 'surface-pchip'].includes(streamwiseInterpolation)) throw new Error('Unknown coupled streamwise refinement interpolation.');
  if (normalFactor !== undefined && normalSubdivisions !== undefined)
    throw new Error('Specify uniform or per-tube normal refinement, not both.');
  const factor = normalFactor ?? (normalSubdivisions === undefined ? 2 : null);
  const stations = nestedIntervalMap(source.euler.layout.nx, { factor: streamwiseFactor, subdivisions: streamwiseSubdivisions });
  streamwiseFactor = stations.uniformFactor;
  const stationMetadata = streamwiseSubdivisions === undefined ? {} : {
    streamwiseSubdivisions: stations.counts, retainedStreamwiseStations: stations.retained };
  const subdivide = (row, interpolate, origin = 0) => Array.from({
    length: stations.retained[origin + row.length - 1] - stations.retained[origin] + 1 },
    (_, i) => sample(row, stations.coordinate(i + stations.retained[origin], origin), interpolate));
  if ((factor !== null && (!Number.isInteger(factor) || factor < 1 || factor > 4))
    || !Number.isInteger(maxNodes) || maxNodes < 1)
    throw new Error('Invalid coupled refinement factors or node budget.');
  // Counts partition each parent's existing mass; they do not change total
  // capture or move a material boundary. A count of one leaves a tube intact.
  const counts = normalSubdivisions ?? source.euler.layout.tubes.map(n => Array(n).fill(factor));
  if (!Array.isArray(counts) || counts.length !== source.euler.layout.tubes.length
    || counts.some((row, g) => !Array.isArray(row) || row.length !== source.euler.layout.tubes[g]
      || row.some(n => !Number.isInteger(n) || n < 1 || n > 4))
    || (streamwiseFactor === 1 && counts.every(row => row.every(n => n === 1))))
    throw new Error('Invalid per-tube coupled refinement subdivisions.');
  if (!source.admissible(initial)) throw new Error('Coupled refinement requires an admissible parent state.');
  const normalCoordinates = counts.map(row => [0, ...row.flatMap((n, j) => Array.from({ length: n }, (_, k) => j + (k + 1) / n))]);
  const count = stations.coordinates.length
    * normalCoordinates.reduce((sum, row) => sum + row.length, 0);
  if (count > maxNodes) throw new Error('Coupled refinement exceeds its node budget.');
  const value = source.evaluate(initial), oldFlow = value.outer, oldBL = value.layers.states;
  const historicalTransonic = source.conditions.blThermodynamics === 'historical-common-isentrope';
  const finiteBase = source.bl.hasFiniteBase === true;
  const physicalEuler = source.euler.conditions.streamwiseMode === 'hybrid'
    || source.euler.conditions.streamwiseMode === 'momentum' && !!source.euler.conditions.upwind;
  const nextInput = structuredClone(input);
  nextInput.outerLower = subdivide(input.outerLower, point);
  nextInput.outerUpper = subdivide(input.outerUpper, point);
  nextInput.cutPaths = input.cutPaths.map(row => subdivide(row, point));
  nextInput.weights = input.weights.map((row, g) => row.flatMap((w, j) => Array(counts[g][j]).fill(w / counts[g][j])));
  const correspondence = streamwiseInterpolation === 'surface-pchip' ? createStreamtubeRefinementCoordinates({
    stations, bodies: source.euler.layout.bodies, fractions: source.euler.fractions, weights: nextInput.weights,
  }) : null;
  nextInput.bodies = input.bodies.map((body, b) => ({ ...structuredClone(body),
    leadingIndex: stations.retained[body.leadingIndex], trailingIndex: stations.retained[body.trailingIndex],
    surfaceFractions: Object.fromEntries(['upper', 'lower'].map(side =>
      [side, correspondence ? correspondence.surfaceFractions[b][side]
        : subdivide(source.euler.fractions[b][side], mix, body.leadingIndex)])) }));
  // Count metadata from a nonnested initializer must not describe this grid.
  nextInput.gridSpacing = { coordinate: 'nested parent intervals', streamwiseFactor, ...stationMetadata, normalFactor: factor,
    normalSubdivisions: structuredClone(counts), normalInterpolation,
    ...(correspondence ? { streamwiseInterpolation } : {}),
    surfaceIntervalsByElement: nextInput.bodies.map((b, element) => ({ element, intervals: b.trailingIndex - b.leadingIndex })) };
  const displacement = { surfaces: nextInput.bodies.map(b => Object.fromEntries(['upper', 'lower'].map(side =>
    [side, Array(b.trailingIndex - b.leadingIndex + 1).fill(0)]))),
    wakes: nextInput.bodies.map(b => Array(nextInput.outerLower.length - 1 - b.trailingIndex).fill(0)) };
  const euler = createStreamtubeBodySystem({ ...nextInput, displacement }), state = euler.initial.slice();
  const oldState = initial.subarray(0, source.ne);
  for (const [key, column] of Object.entries(euler.layout.globals)) {
    if (Array.isArray(column)) column.forEach((col, b) => { if (col !== null) state[col] = oldState[source.euler.layout.globals[key][b]]; });
    else state[column] = oldState[source.euler.layout.globals[key]];
  }
  for (const key of ['lengthScale', 'massScale']) if (euler.conditions[key] !== source.euler.conditions[key])
    throw new Error(`Refinement changed the ${key} normalization.`);
  const base = euler.decode(state).nodes;
  // Interpolate the actual fluid domain. Subdividing from a solid wall
  // through its displacement layer can put new fluid nodes inside that layer.
  const massInterpolation = normalInterpolation === 'streamfunction-quadratic' ? oldFlow.nodes.map((grid, g) =>
    refineStreamtubeMassCoordinates(grid, oldFlow.allocation.groups[g].map(t => t.massFlow), counts[g], {
      lowerStagnation: source.euler.layout.bodies[g - 1]?.leadingIndex ?? null,
      upperStagnation: source.euler.layout.bodies[g]?.leadingIndex ?? null
    })) : null;
  const interpolated = correspondence ? oldFlow.nodes.map((grid, g) => Array.from({ length: euler.layout.nx + 1 }, (_, i) =>
    Array.from({ length: euler.layout.tubes[g] + 1 }, (_, j) => {
      const u = correspondence.nodeCoordinatesByGroup[g][i][j];
      return massInterpolation ? sample(massInterpolation[g].nodes, u, (a, b, t) => point(a[j], b[j], t))
        : sample(grid, u, (a, b, t) => point(sample(a, normalCoordinates[g][j], point), sample(b, normalCoordinates[g][j], point), t));
    }))) : massInterpolation ? massInterpolation.map(({ nodes }) => Array.from({ length: euler.layout.nx + 1 }, (_, i) =>
    sample(nodes, stations.coordinates[i], (a, b, t) => a.map((p, j) => point(p, b[j], t)))))
    : oldFlow.nodes.map((grid, g) => Array.from({ length: euler.layout.nx + 1 }, (_, i) =>
    Array.from({ length: euler.layout.tubes[g] + 1 }, (_, j) => sample(grid, stations.coordinates[i],
      (a, b, t) => point(sample(a, normalCoordinates[g][j], point), sample(b, normalCoordinates[g][j], point), t)))));
  const raw = structuredClone(interpolated);
  for (let b = 0; b < euler.layout.elements; b++) for (let i = 0; i <= euler.layout.nx; i++) {
    if (euler.layout.active(b, i)) {
      raw[b][i][euler.layout.tubes[b]] = base[b][i].at(-1); raw[b + 1][i][0] = base[b + 1][i][0];
    } else if (!euler.layout.independentWakeBanks || i < euler.layout.bodies[b].leadingIndex) {
      // Only shared cuts have one underlying center coordinate. Independent
      // wake banks retain their interpolated physical positions; their gap
      // is a residual, not a coordinate reconstruction from BL thickness.
      const center = point(interpolated[b][i].at(-1), interpolated[b + 1][i][0], .5);
      raw[b][i][euler.layout.tubes[b]] = center; raw[b + 1][i][0] = center;
    }
  }
  let target = euler.adoptGeometry(state, raw);
  const options = { reynolds: source.conditions.reynolds, ncrit: source.conditions.ncrit,
    edgeMatching: source.conditions.edgeMatching, tripFractions: structuredClone(source.bl.trips),
    ...(source.conditions.hkFloorLinearization === undefined ? {} : { hkFloorLinearization: source.conditions.hkFloorLinearization }),
    ...(source.conditions.geometryReplay === undefined ? {} : { geometryReplay: source.conditions.geometryReplay }),
    ...(historicalTransonic ? { blThermodynamics: source.conditions.blThermodynamics } : {}),
    ...(source.bl.transitionMode === 'automatic' ? { transitionMode: 'automatic' } : {}) };
  const bl = createStreamtubeBoundaryLayers(euler, target, historicalTransonic ? { ...options, allowSupersonicEdge: true } : options), geo = bl.geometry(target), values = new Float64Array(4 * bl.stations.length);
  if (options.transitionMode === 'automatic') {
    const transitions = source.evaluate(initial).layers.transitions;
    bl.restoreActive(bl.surfaces.map(branch => {
      const location = transitions.find(t => t.body === branch.body && t.side === branch.side).s;
      const index = branch.ids.findIndex(id => geo.coordinates[id].s >= location);
      return index < 0 ? branch.ids.length - 1 : index;
    }));
  }
  const put = (id, v) => values.set([v.aux, v.theta / bl.scale, v.deltaStar / bl.scale, v.ue], 4 * id);
  const interpolate = (a, b, t) => Object.fromEntries(['aux', 'theta', 'deltaStar', 'ue'].map(key => [key, mix(a[key], b[key], t)]));
  for (const branch of bl.surfaces) {
    const parent = source.bl.surfaces.find(b => b.body === branch.body && b.side === branch.side);
    const parents = parent.ids.map(id => oldBL[id]);
    for (const id of branch.ids) {
      const station = bl.stations[id], le = source.euler.layout.bodies[branch.body].leadingIndex;
      const u = correspondence ? correspondence.surfaceCoordinates[branch.body][branch.side][station.i - stations.retained[le]] - le
        : stations.coordinate(station.i, le);
      let v;
      if (u < 1) {
        // Local stagnation similarity: finite thickness and linear edge speed.
        v = { ...parents[0], aux: 0, ue: parents[0].ue * geo.coordinates[id].s / parents[0].s };
        if (options.transitionMode === 'automatic' && ['leading-transition', 'transition', 'turbulent'].includes(station.regime))
          v.aux = bl.kernel.station({ ...v, s: geo.coordinates[id].s, aux: .03 }, 'turbulent').transitionShear;
      } else {
        v = sample(parents, u - 1, interpolate);
        // N and turbulent shear share a slot but cannot be interpolated
        // across a material trip as though they were the same quantity.
        if (!Number.isInteger(u)) {
          const left = Math.floor(u) - 1, right = left + 1, turbulent = ['leading-transition', 'transition', 'turbulent'].includes(station.regime);
          const isTurbulent = k => ['leading-transition', 'transition', 'turbulent'].includes(source.bl.stations[parent.ids[k]].regime);
          if (isTurbulent(left) !== isTurbulent(right)) v.aux = parents[isTurbulent(left) === turbulent ? left : right].aux;
        }
      }
      put(id, v);
    }
  }
  // The native wake unknown is total displacement, but the viscous profile
  // excludes XICALC's geometric dead-air gap. Interpolate only that fluid
  // displacement and restore the target gap at the target's physical arc.
  // Source decoded gaps are authoritative; never reconstruct them using a
  // different station distance or treat the total profile as a fluid shape.
  const fluidWakeDisplacement = finiteBase ? new Map() : null;
  for (const wake of bl.wakes) {
    const parent = source.bl.wakes.find(b => b.body === wake.body), parents = parent.ids.map(id => {
      const v = oldBL[id];
      return finiteBase ? { ...v, deltaStar: v.deltaStar - (v.wakeGap ?? 0) } : v;
    });
    for (const id of wake.ids) {
      const v = sample(parents,
        stations.coordinate(bl.stations[id].i, source.euler.layout.bodies[wake.body].trailingIndex), interpolate);
      if (finiteBase) {
        if (!(v.deltaStar > v.theta) || !Number.isFinite(v.deltaStar))
          throw new Error('Refined fluid wake displacement is outside the model domain.');
        fluidWakeDisplacement.set(id, v.deltaStar);
        v.deltaStar += geo.coordinates[id].wakeGap ?? 0;
      }
      put(id, v);
    }
  }
  if (options.transitionMode === 'automatic') {
    // Interpolated N does not satisfy the amplification equation on new
    // surface intervals, even when the active interval index is unchanged.
    // Normal-only and wake-only refinement retain the existing surface guess.
    const surfaceRefined = source.euler.layout.bodies.some(b =>
      stations.counts.slice(b.leadingIndex, b.trailingIndex).some(n => n > 1));
    bl.updateActive(values, target, { reinitializeAmplification: surfaceRefined });
    options.transitionState = bl.snapshotActive();
  }
  if (!bl.admissible(values)) throw new Error('Refined BL interpolation is outside the model domain.');
  euler.setDisplacement(bl.thicknesses(values));
  const boundary = euler.decode(target).nodes;
  const corrected = interpolated.map((grid, g) => {
    const weights = nextInput.weights[g], sum = weights.reduce((a, b) => a + b, 0), eta = [0], nj = euler.layout.tubes[g];
    for (const w of weights) eta.push(eta.at(-1) + w / sum); eta[nj] = 1;
    return grid.map((row, i) => {
      // Contour-exact wall offsets and curved wake normals give the boundary
      // correction. Coons interpolation extends it into the fluid seed.
      // Independent wake points above were interpolated in physical space
      // and temporarily adopted with zero displacement. Keep those physical
      // banks here; final adoption removes the new TE chart translation.
      // Applying the decoded translated bank would count that shift twice.
      const translating = euler.layout.wakeDisplacementMotion === 'te-center';
      const lower = translating && g > 0 && i > euler.layout.bodies[g - 1].trailingIndex ? row[0] : boundary[g][i][0];
      const upper = translating && g < euler.layout.elements && i > euler.layout.bodies[g].trailingIndex ? row[nj] : boundary[g][i][nj];
      const dl = { x: lower.x - row[0].x, y: lower.y - row[0].y }, du = { x: upper.x - row[nj].x, y: upper.y - row[nj].y };
      return row.map((p, j) => j === 0 ? { ...lower } : j === nj ? { ...upper } : {
        x: p.x + (1 - eta[j]) * dl.x + eta[j] * du.x, y: p.y + (1 - eta[j]) * dl.y + eta[j] * du.y });
    });
  });
  target = euler.adoptGeometry(target, corrected);
  let finiteBaseTransfer;
  if (finiteBase) {
    // Adoption is the last coordinate operation. In independent-bank mode
    // thickness does not move wakes; in centerline mode its two offsets
    // are symmetric. Both retain the mean-wake arc used by BL geometry.
    const finalGeo = bl.geometry(target);
    let maximumGapReconstructionChange = 0;
    for (const [id, fluid] of fluidWakeDisplacement) {
      const gap = finalGeo.coordinates[id].wakeGap ?? 0;
      maximumGapReconstructionChange = Math.max(maximumGapReconstructionChange,
        Math.abs(gap - (geo.coordinates[id].wakeGap ?? 0)));
      values[4 * id + 2] = (fluid + gap) / bl.scale;
    }
    euler.setDisplacement(bl.thicknesses(values));
    const checkedGeo = bl.geometry(target);
    for (const [id] of fluidWakeDisplacement) {
      const expected = finalGeo.coordinates[id].wakeGap ?? 0;
      const actual = checkedGeo.coordinates[id].wakeGap ?? 0;
      if (Math.abs(expected - actual) > 128 * Number.EPSILON * Math.max(1, Math.abs(expected)))
        throw new Error('Restoring the refined dead-air gap changed its physical wake arc.');
    }
    finiteBaseTransfer = { model: 'XFOIL dead-air gap on final mean-bank arc from solid TE center',
      transferredWakeStations: fluidWakeDisplacement.size, maximumGapReconstructionChange,
      interpolation: 'fluid deltaStar = source total deltaStar minus source decoded wakeGap; target total = interpolated fluid plus final target wakeGap' };
  }
  const geometry = euler.decode(target), mesh = streamtubeMeshSnapshot({ system: euler, nodes: geometry.nodes });
  if (!mesh.quality.valid) throw new Error('Refined displacement grid is not positive.', { cause: mesh.quality });
  let physicalDensityTransfer;
  if (physicalEuler) {
    if (euler.conditions.streamwiseMode !== source.euler.conditions.streamwiseMode) throw new Error('Refinement changed the physical Euler formulation.');
    const transferred = prolongStreamtubeDensities(source.euler, oldState, oldFlow, euler, target,
      { nodeCoordinates: stations.coordinates, subdivisions: counts,
        ...(correspondence ? { nodeCoordinatesByGroup: correspondence.nodeCoordinatesByGroup } : {}) });
    target = transferred.state; physicalDensityTransfer = transferred.diagnostics;
  } else target = initializeStreamtubeDensities(euler, target);
  const initialEuler = { x: target, undisplacedNodes: geometry.undisplacedNodes, nodes: geometry.nodes };
  const system = createCoupledStreamtubeBody(nextInput, { ...options, initialEuler, initialBL: values });
  system.evaluate(system.initial);
  if (!system.admissible(system.initial)) throw new Error('Refined coupled guess is inadmissible.');
  return { input: nextInput, options, initialEuler, initialBL: values, system,
    diagnostics: { streamwiseFactor, ...stationMetadata, normalFactor: factor, normalSubdivisions: structuredClone(counts), normalInterpolation,
      ...(correspondence ? { streamwiseInterpolation, stationCorrespondence: {
        version: 1, surface: 'PCHIP contour fraction; linear parent-interval inverse',
        interior: 'Mass-fraction blend of bounding surface source coordinates',
        density: 'Four-corner mean source coordinate at each child section',
        inactiveBanksAndWakeBL: 'Original linear nested coordinates; passage interiors blend their bounding banks',
      } } : {}),
      massInterpolation: massInterpolation?.map(r => r.diagnostics) ?? null,
      parentUnknowns: source.n, unknowns: system.n, quality: mesh.quality,
      ...(physicalDensityTransfer ? { physicalDensityTransfer } : {}),
      ...(finiteBaseTransfer ? { finiteBaseTransfer } : {}),
      initialization: physicalDensityTransfer
        ? `Nested intervals and conserved child masses; contour-exact boundary interpolation; phase-aware ${historicalTransonic ? 'historical ' : ''}BL interpolation; physical log-density prolongation. Full coupled solve required.`
        : 'Nested intervals and conserved child masses; contour-exact boundary interpolation; phase-aware BL interpolation; subsonic isentropic density inversion. Full coupled solve required.' } };
}
