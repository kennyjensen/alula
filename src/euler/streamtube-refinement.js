// SPDX-License-Identifier: GPL-2.0-or-later
import { prolongStreamtubeDensities } from './streamtube-density-prolongation.js';
// Nested Euler refinement for coarse-to-fine solution and accuracy studies.
// Preserve material banks, parent capture and contour parameters. New states
// are only guesses: initialization does not satisfy the refined equations.
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { initializeStreamtubeDensities } from './streamtube-initial-state.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { refineStreamtubeMassCoordinates } from '../geometry/streamtube-mass-interpolation.js';
import { initialStreamtubeDisplacement } from './streamtube-geometry.js';

const mix = (a, b, t) => t === 0 ? a : t === 1 ? b : a + t * (b - a);
const point = (a, b, t) => ({ x: mix(a.x, b.x, t), y: mix(a.y, b.y, t) });
const sample = (row, u, interpolate) => {
  const i = Math.min(row.length - 2, Math.floor(u));
  return interpolate(row[i], row[i + 1], u - i);
};
const subdivide = (row, factor, interpolate) => Array.from({ length: factor * (row.length - 1) + 1 }, (_, i) => sample(row, i / factor, interpolate));

export function refineStreamtubeBody(input, source, { initial = source.initial,
  streamwiseFactor = 2, normalSubdivisions = source.layout.tubes.map(n => Array(n).fill(2)), maxNodes = 50000,
  initializeFlow = true, normalInterpolation = 'linear' } = {}) {
  if (!['linear', 'streamfunction-quadratic'].includes(normalInterpolation))
    throw new Error('Unknown Euler refinement interpolation.');
  const finiteBaseWake = source.inviscidBaseWake === true;
  if (source.layout.displacedBoundaries && !finiteBaseWake || !source.layout.densityCount)
    throw new Error('Euler refinement requires a compressible system without displacement layers. Use coupled refinement for BL states.');
  if (!Number.isInteger(streamwiseFactor) || streamwiseFactor < 1 || streamwiseFactor > 4
    || !Number.isInteger(maxNodes) || maxNodes < 1
    || !Array.isArray(normalSubdivisions) || normalSubdivisions.length !== source.layout.tubes.length
    || normalSubdivisions.some((row, g) => !Array.isArray(row) || row.length !== source.layout.tubes[g]
      || row.some(n => !Number.isInteger(n) || n < 1 || n > 4))
    || streamwiseFactor === 1 && normalSubdivisions.every(row => row.every(n => n === 1)))
    throw new Error('Invalid Euler refinement subdivisions or node budget.');
  if (finiteBaseWake) {
    // The finite solid base enables the displacement chart even in an
    // inviscid system. Admit only that original zero-wall, constant-width
    // profile: setDisplacement can otherwise change it after construction.
    const expected = initialStreamtubeDisplacement(source.layout, source.baseGeometry), actual = source.displacement;
    const same = (a, b) => Array.isArray(a) && a.length === b.length && a.every((v, i) => v === b[i]);
    if (input.displacement !== undefined || !source.layout.displacedBoundaries
      || !Array.isArray(actual?.surfaces) || actual.surfaces.length !== expected.surfaces.length
      || actual.surfaces.some((surface, b) => ['upper', 'lower'].some(side => !same(surface?.[side], expected.surfaces[b][side])))
      || !Array.isArray(actual?.wakes) || actual.wakes.length !== expected.wakes.length
      || actual.wakes.some((row, b) => !same(row, expected.wakes[b])))
      throw new Error('Finite-base Euler refinement requires unchanged zero-wall, constant-width inviscid wake data. Use coupled refinement for BL states.');
    // Adding only crossline points retains every center and wake-bank
    // secant, hence its prescribed normal gap. New streamwise stations
    // would require a separate gap-preserving correspondence construction.
    if (streamwiseFactor !== 1)
      throw new Error('Finite-base Euler refinement currently supports normal subdivisions only (streamwiseFactor: 1); streamwise wake-gap transfer is not implemented.');
  }
  const coordinates = normalSubdivisions.map(row => [0, ...row.flatMap((n, j) => Array.from({ length: n }, (_, k) => j + (k + 1) / n))]);
  const nodeCount = (streamwiseFactor * source.layout.nx + 1) * coordinates.reduce((sum, row) => sum + row.length, 0);
  if (nodeCount > maxNodes) throw new Error('Euler refinement exceeds its node budget.');
  // Mesh preview must not need an admissible gas state or enter Newton.
  const parent = initializeFlow ? source.evaluate(initial) : source.decode(initial), nextInput = structuredClone(input);
  nextInput.outerLower = subdivide(input.outerLower, streamwiseFactor, point);
  nextInput.outerUpper = subdivide(input.outerUpper, streamwiseFactor, point);
  nextInput.cutPaths = input.cutPaths.map(row => subdivide(row, streamwiseFactor, point));
  nextInput.weights = input.weights.map((row, g) => row.flatMap((w, j) => Array(normalSubdivisions[g][j]).fill(w / normalSubdivisions[g][j])));
  nextInput.bodies = input.bodies.map((body, b) => ({ ...structuredClone(body),
    leadingIndex: streamwiseFactor * body.leadingIndex, trailingIndex: streamwiseFactor * body.trailingIndex,
    surfaceFractions: Object.fromEntries(['upper', 'lower'].map(side => [side, subdivide(source.fractions[b][side], streamwiseFactor, mix)])) }));
  nextInput.gridSpacing = { coordinate: 'nested parent intervals', streamwiseFactor,
    normalSubdivisions: structuredClone(normalSubdivisions), normalInterpolation,
    surfaceIntervalsByElement: nextInput.bodies.map((b, body) => ({ element: b.element ?? body, intervals: b.trailingIndex - b.leadingIndex })) };
  const system = createStreamtubeBodySystem(nextInput), state = system.initial.slice();
  for (const key of ['lengthScale', 'massScale']) if (system.conditions[key] !== source.conditions[key])
    throw new Error(`Euler refinement changed the ${key} normalization.`);
  for (const [key, column] of Object.entries(system.layout.globals)) {
    if (Array.isArray(column)) column.forEach((col, b) => { if (col !== null) state[col] = initial[source.layout.globals[key][b]]; });
    else state[column] = initial[source.layout.globals[key]];
  }
  const boundary = system.decode(state).nodes;
  // A geometric midpoint need not carry half the parent mass, particularly
  // at stagnation where streamfunction varies quadratically with distance.
  // Locate child mass levels on the existing edges before adding streamwise
  // stations. This is the same independently checked transfer used for BL
  // refinement; it changes neither parent nodes nor captured/tube masses.
  const massInterpolation = normalInterpolation === 'streamfunction-quadratic' ? parent.nodes.map((grid, g) =>
    refineStreamtubeMassCoordinates(grid, parent.allocation.groups[g].map(t => t.massFlow), normalSubdivisions[g], {
      lowerStagnation: source.layout.bodies[g - 1]?.leadingIndex ?? null,
      upperStagnation: source.layout.bodies[g]?.leadingIndex ?? null,
    })) : null;
  const nodes = parent.nodes.map((grid, g) => Array.from({ length: system.layout.nx + 1 }, (_, i) => {
    const row = massInterpolation
      ? sample(massInterpolation[g].nodes, i / streamwiseFactor, (a, b, t) => a.map((p, j) => point(p, b[j], t)))
      : coordinates[g].map(u => sample(grid, i / streamwiseFactor, (a, b, t) => point(sample(a, u, point), sample(b, u, point), t)));
    // Restore the exact contour at new surface stations. Spread only this
    // interpolation correction into the passage; original stations stay put.
    const lower = g > 0 && system.layout.active(g - 1, i) ? boundary[g][i][0] : row[0];
    const upper = g < system.layout.elements && system.layout.active(g, i) ? boundary[g][i].at(-1) : row.at(-1);
    const dl = { x: lower.x - row[0].x, y: lower.y - row[0].y }, du = { x: upper.x - row.at(-1).x, y: upper.y - row.at(-1).y };
    return row.map((p, j) => {
      const t = coordinates[g][j] / source.layout.tubes[g];
      return j === 0 ? { ...lower } : j === row.length - 1 ? { ...upper }
        : { x: p.x + (1 - t) * dl.x + t * du.x, y: p.y + (1 - t) * dl.y + t * du.y };
    });
  }));
  let target = system.adoptGeometry(state, nodes);
  const mesh = streamtubeMeshSnapshot({ system, nodes: system.decode(target).nodes });
  if (!mesh.quality.valid) throw new Error('Refined Euler seed is not a positive convex grid.', { cause: mesh.quality });
  // Hybrid and explicitly upwind momentum flow retain physical density.
  // Conserved child masses and final normal areas define q. The existing
  // nonhybrid cold inversion remains unchanged; full refinement equations
  // must still be solved in either case.
  let physicalDensityTransfer;
  if (initializeFlow && (source.conditions.streamwiseMode === 'hybrid'
    || source.conditions.streamwiseMode === 'momentum' && !!source.conditions.upwind)) {
    if (system.conditions.streamwiseMode !== source.conditions.streamwiseMode) throw new Error('Refinement changed the physical Euler formulation.');
    const transferred = prolongStreamtubeDensities(source, initial, parent, system, target,
      { nodeCoordinates: Array.from({ length: system.layout.nx + 1 }, (_, i) => i / streamwiseFactor),
        subdivisions: normalSubdivisions });
    target = transferred.state; physicalDensityTransfer = transferred.diagnostics;
  } else if (initializeFlow) target = initializeStreamtubeDensities(system, target);
  const value = initializeFlow ? system.evaluate(target) : system.decode(target);
  return { input: nextInput, system, initial: target, initialEuler: { x: target, nodes: value.nodes },
    diagnostics: { streamwiseFactor, normalSubdivisions: structuredClone(normalSubdivisions), normalInterpolation,
      massInterpolation: massInterpolation?.map(r => r.diagnostics) ?? null,
      parentUnknowns: source.layout.n, unknowns: system.layout.n, nodeCount, quality: mesh.quality,
      flowInitialized: initializeFlow,
      ...(physicalDensityTransfer ? { physicalDensityTransfer } : {}),
      ...(finiteBaseWake ? { finiteBaseWake: { model: 'constant-width inviscid', retainedStreamwiseStations: true,
        widths: system.baseGeometry.map(base => base?.width ?? 0) } } : {}),
      ...(initializeFlow ? { residual: value.diagnostics.residual, maxMach: value.diagnostics.maxMach } : {}),
      initialization: 'Nested parent intervals; conserved child masses; exact contour banks. '
        + (physicalDensityTransfer ? 'Physical log-density prolongation; refined flow equations still require solving.'
          : initializeFlow ? 'Subsonic isentropic density inversion; refined flow equations still require solving.' : 'Geometry only; gas initialization and flow solution have not run.') } };
}
