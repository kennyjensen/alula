// SPDX-License-Identifier: GPL-2.0-or-later
import { streamtubeMeshConnectivity } from '../geometry/streamtube-mesh-connectivity.js';
import { untangleQuadrilaterals } from '../geometry/untangle-quadrilaterals.js';
import { smoothQuadrilaterals } from '../geometry/smooth-quadrilaterals.js';
import { isentropicSonicMassFlux } from './streamtube-initial-state.js';
import { streamtubeCellGeometry } from './streamtube-cell.js';
import { streamtubeGridConvexity } from '../geometry/streamtube-convex-step.js';
import { segmentsTouch } from '../geometry/airfoil.js';

export function repairStreamtubeGrid({ system, nodes }, controls) {
  const mesh = streamtubeMeshConnectivity(system.layout, nodes), repaired = untangleQuadrilaterals(mesh, controls);
  return { ...repaired, nodes: mesh.indices.map(group => group.map(row => row.map(id => ({ ...repaired.vertices[id] })))) };
}

// Initializer geometry only: require enough normal width for every assigned
// tube mass to admit the subsonic density inversion. The fraction supplies
// explicit headroom below choking; it changes neither Euler equations nor
// their Mach/admissibility checks. Pressure positivity must still be checked
// by the full system after gas initialization.
export function repairSubsonicStreamtubeGrid({ system, nodes, initial = system.initial }, {
  massFluxFraction = .99, ...controls
} = {}) {
  if (system.conditions.flowModel !== 'compressible' || !Number.isFinite(massFluxFraction)
    || massFluxFraction <= 0 || massFluxFraction >= 1)
    throw new Error('Invalid subsonic grid-repair conditions.');
  const massFluxLimit = massFluxFraction * isentropicSonicMassFlux(system.conditions);
  const { allocation } = system.decode(initial), { nx, tubes } = system.layout;
  const minimumNormalAreas = tubes.flatMap((count, g) => Array.from({ length: nx }, () =>
    Array.from({ length: count }, (_, j) => {
      const mass = allocation.groups[g][j].massFlow;
      if (!(mass > 0) || !Number.isFinite(mass)) throw new Error('Invalid assigned mass in subsonic grid repair.');
      return mass / massFluxLimit;
    })).flat());
  const repaired = repairStreamtubeGrid({ system, nodes }, { ...controls, minimumNormalAreas });
  return { ...repaired, massFluxLimit, massFluxFraction, minimumNormalAreas,
    status: 'corner/section-width feasibility only; pressure and flow equations have not been solved' };
}

export function smoothStreamtubeGrid({ system, nodes }, controls) {
  const mesh = streamtubeMeshConnectivity(system.layout, nodes), smoothed = smoothQuadrilaterals(mesh, controls);
  return { ...smoothed, nodes: mesh.indices.map(group => group.map(row => row.map(id => ({ ...smoothed.vertices[id] })))) };
}

// Additional geometry-only gate for a retained partial SLOR initial guess.
// No gas state is evaluated and no coordinates are repaired here.
export function requireStreamtubeInitialGridDomain(nodes) {
  const quality = streamtubeGridConvexity(nodes);
  if (!quality.valid) throw new Error('Retained initial grid is folded or degenerate.');
  let localStencils = 0;
  for (let g = 0; g < nodes.length; g++) {
    if (nodes[g].length < 3) throw new Error('Retained initial grid needs at least two streamwise intervals.');
    for (let i = 1; i < nodes[g].length - 1; i++) for (let j = 0; j < nodes[g][i].length - 1; j++) {
      streamtubeCellGeometry([nodes[g][i - 1][j], nodes[g][i][j], nodes[g][i + 1][j]],
        [nodes[g][i - 1][j + 1], nodes[g][i][j + 1], nodes[g][i + 1][j + 1]], { geometryDomain: 'positive-simple' });
      localStencils++;
    }
  }
  // Coordinates, rather than solver ids, identify coincident cuts/banks.
  // Two cells may traverse a shared edge only in opposite directions.
  const key = p => `${p.x},${p.y}`, edges = new Map();
  for (const region of nodes) for (let i = 0; i < region.length - 1; i++) for (let j = 0; j < region[i].length - 1; j++) {
    const p = [region[i][j], region[i + 1][j], region[i + 1][j + 1], region[i][j + 1]];
    for (let k = 0; k < 4; k++) {
      const a = p[k], b = p[(k + 1) % 4], ids = [key(a), key(b)].sort(), id = ids.join('|'), sign = key(a) === ids[0] ? 1 : -1;
      const old = edges.get(id);
      if (old) {
        if (old.count !== 1 || old.sign === sign) throw new Error('Retained initial grid has nonmanifold or overlapping cells.');
        old.count++; continue;
      }
      edges.set(id, { a, b, ids, sign, count: 1, minX: Math.min(a.x, b.x), maxX: Math.max(a.x, b.x),
        minY: Math.min(a.y, b.y), maxY: Math.max(a.y, b.y) });
    }
  }
  let active = [], testedPairs = 0;
  for (const current of [...edges.values()].sort((a, b) => a.minX - b.minX)) {
    active = active.filter(e => e.maxX >= current.minX);
    for (const e of active) {
      if (e.maxY < current.minY || current.maxY < e.minY) continue;
      testedPairs++;
      const shared = e.ids.find(v => current.ids.includes(v));
      if (shared !== undefined) {
        const a = key(e.a) === shared ? e.a : e.b, b = key(e.a) === shared ? e.b : e.a;
        const c = key(current.a) === shared ? current.b : current.a;
        if ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) === 0
          && (b.x - a.x) * (c.x - a.x) + (b.y - a.y) * (c.y - a.y) > 0)
          throw new Error('Retained initial grid has overlapping edges.');
      } else if (segmentsTouch(e.a, e.b, current.a, current.b, 0))
        throw new Error('Retained initial grid has crossed or touching nonadjacent edges.');
    }
    active.push(current);
  }
  return { valid: true, minCornerSine: quality.minCornerSine, localStencils,
    edges: edges.size, testedPairs, primalHalfDualVolumes: 'positive simple', embedding: 'no edge crossings or nonmanifold edges' };
}
