// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour } from '../geometry/airfoil.js';
import { createDisplacementOperator } from '../inviscid/displacement.js';
import { triangularAirfoilMesh } from '../euler/triangular-mesh.js';

export const subcriticalMeshOptions = options => ({ type: 'triangular', surfaceScale: .3, growth: .35, padding: 8, ...options });

// Copy display geometry only. Never pass operators or mutable solver state
// across the worker boundary, or normalize the UI's coordinates twice.
export function potentialMeshSnapshot(mesh, referenceChord = 1) {
  return { topology: mesh.topology,
    vertices: mesh.vertices.map(p => ({ x: p.x * referenceChord, y: p.y * referenceChord })),
    cells: mesh.cells.map(cell => cell.vertices.slice()) };
}

// An initial computational mesh with inviscid wake guides. No boundary-layer
// or nonlinear flow solve is needed. The analysis can subsequently deform or
// rebuild this geometry as its viscous wakes change.
export function buildSubcriticalMeshPreview({ elements, alpha = 0, referenceChord = 1,
  wakeLength = 2, wakeLengths, wakeCount = 48, wakePaths } = {}, { mesh: meshOptions } = {}) {
  if (!Array.isArray(elements) || elements.length < 1 || elements.length > 6
    || !Number.isFinite(alpha) || Math.abs(alpha) > 20 || !Number.isFinite(referenceChord) || referenceChord <= 0
    || !Number.isFinite(wakeLength) || wakeLength <= 0 || !Number.isInteger(wakeCount) || wakeCount < 4 || wakeCount > 96)
    throw new Error('Invalid subcritical mesh geometry or wake controls.');
  const contours = elements.map(e => prepareContour(e.points));
  if (contours.some(p => p.length < 41) || contours.reduce((n, p) => n + p.length - 1, 0) > 700)
    throw new Error('Use at least 40 panels per element and at most 700 total.');
  if (wakeLengths && (!Array.isArray(wakeLengths) || wakeLengths.length !== elements.length || wakeLengths.some(v => !Number.isFinite(v) || v <= 0)))
    throw new Error('Supply one positive absolute wake length per element.');
  const normalized = contours.map(points => ({ points: points.map(p => ({ x: p.x / referenceChord, y: p.y / referenceChord })) }));
  const target = Math.max(...normalized.flatMap(e => e.points.map(p => p.x))) + wakeLength;
  const outer = createDisplacementOperator({ elements: normalized, alpha, wakeCount, computeInfluence: false, wakeInitialization: 'inviscid',
    wakeLengths: wakeLengths ? wakeLengths.map(v => v / referenceChord) : normalized.map(e => target - e.points[0].x),
    wakePaths: wakePaths?.map(path => path.map(p => ({ x: p.x / referenceChord, y: p.y / referenceChord }))) });
  const { type, ...options } = subcriticalMeshOptions(meshOptions);
  if (type !== 'triangular') throw new Error('The standalone subcritical mesh preview requires the triangular mesher.');
  const mesh = triangularAirfoilMesh(outer.bodies.map(b => b.points), { ...options, wakePaths: outer.wakes.map(w => w.points) });
  return potentialMeshSnapshot(mesh, referenceChord);
}
