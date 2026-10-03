// SPDX-License-Identifier: GPL-2.0-or-later
// Nested reference refinement with actual boundary-curve sampling. Coons
// blending extends boundary changes into the seed; a later elliptic solve
// determines the interior. This neither repairs cells nor relaxes a field.
import { potentialGridQuality } from '../potential-plane-grid.js';

export function refinePotentialGrid({ nodes, coordinates }, { boundaryAt, factor = 2, maxNodes = 100000 } = {}) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  const shape = grid => Array.isArray(grid) && grid.length === nx + 1 && grid.every(row => Array.isArray(row)
    && row.length === nt + 1 && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)));
  if (!(nx >= 2 && nt >= 2) || !shape(nodes) || !shape(coordinates) || typeof boundaryAt !== 'function'
    || !Number.isInteger(factor) || factor < 2 || factor > 8 || !Number.isInteger(maxNodes) || maxNodes < 1)
    throw new Error('Invalid potential-grid refinement controls.');
  const ni = factor * nx, nj = factor * nt;
  if ((ni + 1) * (nj + 1) > maxNodes) throw new Error('Potential-grid refinement exceeds its node budget.');
  if (!potentialGridQuality(nodes).valid || !potentialGridQuality(coordinates).valid)
    throw new Error('Potential-grid refinement requires positive parent cells.');
  const levels = coordinates[0].map(p => p.y), span = levels.at(-1) - levels[0];
  const psiScale = Math.max(Math.abs(levels[0]), Math.abs(levels.at(-1)), span);
  if (!(span > 0) || !Number.isFinite(span) || levels.some((v, j) => j && !(v > levels[j - 1]))
    || coordinates.some(row => row.some((p, j) => Math.abs(p.y - levels[j]) > 64 * Number.EPSILON * psiScale)))
    throw new Error('Reference refinement requires ordered constant-streamfunction rows.');
  const sample = (grid, i, j) => {
    const a = Math.min(nx - 1, Math.floor(i / factor)), b = Math.min(nt - 1, Math.floor(j / factor));
    const s = i / factor - a, t = j / factor - b, p = grid[a][b], q = grid[a + 1][b], r = grid[a + 1][b + 1], v = grid[a][b + 1];
    return Object.fromEntries(['x', 'y'].map(key => [key, p[key] + s * (q[key] - p[key]) + t * (v[key] - p[key])
      + s * t * (r[key] - q[key] - v[key] + p[key])]));
  };
  const physical = [], plane = [], changes = [];
  const finite = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
  for (let i = 0; i <= ni; i++) {
    physical.push([]); plane.push([]); changes.push([]);
    for (let j = 0; j <= nj; j++) {
      const p = sample(nodes, i, j), w = sample(coordinates, i, j);
      // Exact coarse-node copies make the nested correspondence explicit,
      // including far-boundary endpoint roundoff in the traced data.
      if (i % factor === 0 && j % factor === 0) {
        physical[i].push({ ...nodes[i / factor][j / factor] }); plane[i].push({ ...coordinates[i / factor][j / factor] });
        changes[i].push({ x: 0, y: 0, phi: 0 }); continue;
      }
      if (!i || i === ni || !j || j === nj) {
        const value = boundaryAt(i / factor, j / factor);
        if (!finite(value?.point) || !finite(value?.coordinate) || Math.abs(value.coordinate.y - w.y) > 64 * Number.EPSILON * psiScale)
          throw new Error('Refined boundary must preserve the prescribed streamfunction labels.');
        physical[i].push({ ...value.point }); plane[i].push({ x: value.coordinate.x, y: w.y });
        changes[i].push({ x: value.point.x - p.x, y: value.point.y - p.y, phi: value.coordinate.x - w.x });
      } else { physical[i].push(p); plane[i].push(w); changes[i].push(null); }
    }
  }
  for (let i = 1; i < ni; i++) for (let j = 1; j < nj; j++) {
    if (i % factor === 0 && j % factor === 0) continue;
    const u = i / ni, v = (plane[i][j].y - levels[0]) / span;
    const correction = key => (1 - v) * changes[i][0][key] + v * changes[i][nj][key]
      + (1 - u) * changes[0][j][key] + u * changes[ni][j][key];
    physical[i][j].x += correction('x'); physical[i][j].y += correction('y'); plane[i][j].x += correction('phi');
  }
  return { nodes: physical, coordinates: plane, factor, quality: potentialGridQuality(physical), computationalQuality: potentialGridQuality(plane),
    method: 'Nested physical-boundary sampling and Coons seed interpolation; old physical nodes, potential labels and streamfunction levels preserved. No cell repair.' };
}
