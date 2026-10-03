// SPDX-License-Identifier: GPL-2.0-or-later
import { auditHarmonicGrid } from './harmonic-grid-audit.js';

export function streamtubeAuditRegions(mesh) {
  const { streamwiseSegments: nx, tubes, massFlows } = mesh?.initialization ?? {};
  if (mesh?.topology !== 'intrinsic-quadrilateral-streamtubes' || !Number.isInteger(nx) || nx < 2
    || !Array.isArray(tubes) || !Array.isArray(massFlows) || tubes.length !== massFlows.length
    || tubes.some((n, g) => !Number.isInteger(n) || n < 2 || massFlows[g].length !== n)
    || mesh.cells.length !== nx * tubes.reduce((s, n) => s + n, 0)) throw new Error('Build an initial quad grid with recorded tube masses before checking it.');
  const recordedXi = mesh.initialization.gridSmoothing?.stationCoordinate?.xi;
  if (recordedXi !== undefined && (!Array.isArray(recordedXi) || recordedXi.length !== nx + 1
    || recordedXi[0] !== 0 || recordedXi[nx] !== 1
    || Array.from(recordedXi).some((v, i) => !Number.isFinite(v) || i && !(v > recordedXi[i - 1]))))
    throw new Error('Recorded streamwise coordinates must increase strictly from zero to one and match the grid.');
  const streamwiseCoordinates = recordedXi === undefined ? undefined : Object.freeze([...recordedXi]);
  const recordedRegions = mesh.initialization.gridSmoothing?.regions;
  if (recordedRegions !== undefined && (!Array.isArray(recordedRegions) || recordedRegions.length !== tubes.length))
    throw new Error('Recorded boundary conditions must match the grid regions.');
  let offset = 0;
  return tubes.map((nt, g) => {
    const recorded = recordedRegions?.[g]?.coordinateEquations?.boundaryConditions;
    if (recorded !== undefined && (!recorded || typeof recorded !== 'object' || Array.isArray(recorded)
      || Object.keys(recorded).some(key => !['lower', 'upper'].includes(key))))
      throw new Error('Invalid recorded boundary conditions.');
    const boundaryConditions = Object.freeze(Object.fromEntries(['lower', 'upper'].map(side => {
      const mode = recorded?.[side] === undefined ? 'fixed' : recorded[side];
      if (!['fixed', 'giles-vertical', 'normal-curve'].includes(mode)) throw new Error('Unknown recorded boundary condition.');
      return [side, mode];
    })));
    const boundaryCurves = recordedRegions?.[g]?.coordinateEquations?.boundaryCurves;
    const indices = Array.from({ length: nx + 1 }, () => Array(nt + 1).fill(null)), cellIds = [];
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
      const id = offset + i * nt + j, cell = mesh.cells[id];
      if (!Array.isArray(cell) || cell.length !== 4) throw new Error('Invalid quad grid connectivity.');
      [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]].forEach(([a, b], k) => {
        if (indices[a][b] !== null && indices[a][b] !== cell[k]) throw new Error('Disconnected quad grid audit region.');
        indices[a][b] = cell[k];
      }); cellIds.push(id);
    }
    offset += nx * nt;
    return { nodes: indices.map(row => row.map(id => mesh.vertices[id])), indices, cellIds, massFlows: massFlows[g], boundaryConditions,
      ...(boundaryCurves === undefined ? {} : { boundaryCurves: structuredClone(boundaryCurves) }),
      ...(streamwiseCoordinates === undefined ? {} : { streamwiseCoordinates }) };
  });
}

export function measureStreamtubeGridSpacing(regions, { adjacentRatioLimit = 1.5, cornerSineLimit = .1 } = {}) {
  if (!Number.isFinite(adjacentRatioLimit) || adjacentRatioLimit < 1 || !Number.isFinite(cornerSineLimit) || !(cornerSineLimit > 0 && cornerSineLimit <= 1))
    throw new Error('Invalid grid spacing screening limits.');
  const flagged = new Map(); let maximumAdjacentRatio = 1, minimumCornerSine = 1;
  const worstSpacing = [], worstCorners = [];
  const flag = (id, kind, score) => {
    const old = flagged.get(id) ?? { cell: id, spacing: 1, shear: 0 }; old[kind] = Math.max(old[kind], score); flagged.set(id, old);
  };
  regions.forEach(({ nodes, cellIds, boundaryConditions }, g) => {
    const nx = nodes.length - 1, nt = nodes[0].length - 1;
    for (let j = 0; j <= nt; j++) {
      const length = nodes.slice(1).map((row, i) => Math.hypot(row[j].x - nodes[i][j].x, row[j].y - nodes[i][j].y));
      for (let i = 1; i < nx; i++) {
        const ratio = Math.max(length[i] / length[i - 1], length[i - 1] / length[i]);
        maximumAdjacentRatio = Math.max(maximumAdjacentRatio, ratio);
        if (ratio > adjacentRatioLimit) {
          const fixedBoundary = (j === 0 && (boundaryConditions?.lower ?? 'fixed') === 'fixed')
            || (j === nt && (boundaryConditions?.upper ?? 'fixed') === 'fixed');
          worstSpacing.push({ region: g, i, j, ratio, point: nodes[i][j], fixedBoundary });
          for (const a of [i - 1, i]) for (const b of [j - 1, j]) if (b >= 0 && b < nt) flag(cellIds[a * nt + b], 'spacing', ratio);
        }
      }
    }
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
      const p = [nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1]], id = cellIds[i * nt + j];
      let minimum = 1;
      for (let k = 0; k < 4; k++) {
        const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4], x = b.x - a.x, y = b.y - a.y, u = c.x - b.x, v = c.y - b.y;
        minimum = Math.min(minimum, (x * v - y * u) / (Math.hypot(x, y) * Math.hypot(u, v)));
      }
      minimumCornerSine = Math.min(minimumCornerSine, minimum);
      if (!(minimum >= cornerSineLimit)) { worstCorners.push({ region: g, i, j, cell: id, sine: minimum }); flag(id, 'shear', 1 - minimum); }
    }
  });
  worstSpacing.sort((a, b) => b.ratio - a.ratio); worstCorners.sort((a, b) => a.sine - b.sine);
  return { limits: { adjacentRatioLimit, cornerSineLimit }, maximumAdjacentRatio, minimumCornerSine,
    flaggedCells: [...flagged.values()], spacingLocations: worstSpacing, shearLocations: worstCorners,
    scope: 'Spacing/shear screening limits are explicit engineering targets, not universal physical laws. Boundary-imposed clustering is included.' };
}

export function auditStreamtubeMesh(mesh, { onProgress, ...controls } = {}) {
  if (mesh.initialization?.flowSolved) throw new Error('This audit checks incompressible initial grids, not a converged compressible flow grid.');
  const regions = streamtubeAuditRegions(mesh), spacing = measureStreamtubeGridSpacing(regions), harmonic = [];
  regions.forEach((region, g) => {
    onProgress?.({ region: g, groups: regions.length });
    const r = auditHarmonicGrid(region, controls); harmonic.push({ region: g, ...r });
  });
  const fieldErrors = mesh.cells.map(() => ({ crosslineIntervals: 0, tubeIntervals: 0 }));
  regions.forEach((r, g) => {
    const nx = r.nodes.length - 1, nt = r.nodes[0].length - 1;
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) for (const key of ['crosslineIntervals', 'tubeIntervals'])
      fieldErrors[r.cellIds[i * nt + j]][key] = Math.max(...[harmonic[g].errors[i][j], harmonic[g].errors[i + 1][j], harmonic[g].errors[i + 1][j + 1], harmonic[g].errors[i][j + 1]].map(p => Math.abs(p[key])));
  });
  return { status: harmonic.every(r => r.passed) && !spacing.flaggedCells.length ? 'selected interior and spacing checks passed; boundary/flow validation pending' : 'grid checks need attention',
    accepted: false, scope: 'Independent initial-grid diagnostic. General boundary placement, cut velocity continuity and physical Euler/BL validation remain required.',
    spacing, harmonic: harmonic.map(({ values, errors, ...r }) => r), fieldErrors };
}
