// SPDX-License-Identifier: GPL-2.0-or-later
// Connect the independently reconstructed isolated ISET geometry to the
// moving-quadrilateral chart. No panel, Euler, BL or density solve occurs here.
import { createIsetAirfoilInitialGrid } from '../../geometry/tests/iset-airfoil-initializer.js';
import { createStreamtubeBodySystem } from '../streamtube-body.js';
import { streamtubeMeshSnapshot, prepareStreamtubeMesh } from '../streamtube-mesh-preview.js';
import { pointInside } from '../../geometry/airfoil.js';

export function createIsetStreamtubeGrid(source, { alpha = 0, mach = .2, gamma = 1.4,
  flowModel = 'compressible', smoothing = true, onMesh } = {}) {
  if (typeof smoothing !== 'boolean' && (!smoothing || typeof smoothing !== 'object' || Array.isArray(smoothing)))
    throw new Error('ISET smoothing must be a boolean or an object of SLOR controls.');
  const grid = createIsetAirfoilInitialGrid(source);
  const { resampled, leadingIndex, trailingIndex } = grid;
  if (resampled.points.length < 9)
    throw new Error('The quadrilateral body chart requires at least five ISET stations per surface.');
  const input = { bodies: [{ points: resampled.points, leadingIndex, trailingIndex,
    stagnationParameter: resampled.leadingParameter, surfaceFractions: resampled.surfaceFractions, element: 0 }],
    primaryBody: 0, cutPaths: [grid.outline.upper.map((p, i) => ({
      x: .5 * (p.x + grid.outline.lower[i].x), y: .5 * (p.y + grid.outline.lower[i].y) }))],
    outerLower: grid.farfield.lower, outerUpper: grid.farfield.upper, weights: grid.massFlows,
    alpha, mach, gamma, flowModel };
  const system = createStreamtubeBodySystem(input);
  // Adopt the actual RESPLI geometry and explicit mass-coordinate seed.
  // Chart construction may differ by roundoff at spline knots; decode once
  // to use precisely identical shared endpoints in snapshots and restarts.
  const initial = system.adoptGeometry(system.initial, grid.nodes), nodes = system.decode(initial).nodes;
  const diagnostics = { initializer: 'ISET isolated-airfoil reconstruction', iset: grid.diagnostics,
    gridSpacing: { coordinate: 'supplied physical inlet/outlet and contour-parameter stations',
      inlet: { intervals: leadingIndex }, outlet: { intervals: grid.nodes[0].length - 1 - trailingIndex } },
    flowSolved: false };
  const guideField = {
    admissibleNode: p => !pointInside(p, resampled.points),
    diagnosticsForNodes: () => structuredClone(diagnostics),
  };
  const prepared = { input, system, initial, nodes, diagnostics, guideField,
    sourceInitialization: grid, status: 'ISET geometry only; Euler/BL equations not initialized or solved' };
  if (smoothing) return prepareStreamtubeMesh(prepared, { onMesh, ellipticSmoothing: {
    ...(typeof smoothing === 'object' ? smoothing : {}), discretization: 'giles-1985',
    streamwiseCoordinates: grid.xi, farfieldBoundary: 'giles-indexed-y',
  } });
  const mesh = streamtubeMeshSnapshot(prepared);
  onMesh?.(mesh, 'initial');
  return { ...prepared, mesh };
}
