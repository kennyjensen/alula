// SPDX-License-Identifier: GPL-2.0-or-later
// Independent physical-space Q1 finite-element Laplace audit. No SLOR
// metrics, stencil or residuals are reused. Polygonal boundaries are fixed
// under refinement unless a prescribed normal-curve graph is supplied. Those
// graphs are approximated by refined Q1 polygons, not curved-cell integration.
import { sparseMatrix, sparseAdd } from '../numerics/sparse.js';
import { solveSparseDirect } from '../numerics/klu.js';
import { createNormalGraphBoundary } from './normal-graph-boundary.js';

import { quadLaplaceMatrix } from '../numerics/quad-laplace.js';
export { quadLaplaceMatrix } from '../numerics/quad-laplace.js';
const maximum = a => Math.max(0, ...a.map(Math.abs));

export function solveHarmonicGridReference({ nodes, massFlows, streamwiseCoordinates, boundaryConditions = {}, boundaryCurves = {} }, { refinement = 1, maxUnknowns = 50000, includeRefinedField = false } = {}) {
  const nx = nodes?.length - 1, nt = massFlows?.length;
  if (!Number.isInteger(refinement) || refinement < 1 || refinement > 8 || !Number.isInteger(maxUnknowns) || maxUnknowns < 1
    || !Array.isArray(nodes) || nx < 2 || !Array.isArray(massFlows) || nt < 2
    || massFlows.some(m => !(m > 0) || !Number.isFinite(m))
    || nodes.some(row => !Array.isArray(row) || row.length !== nt + 1 || row.some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y))))
    throw new Error('Invalid harmonic grid audit input.');
  if (streamwiseCoordinates !== undefined && (!Array.isArray(streamwiseCoordinates) || streamwiseCoordinates.length !== nx + 1
    || streamwiseCoordinates[0] !== 0 || streamwiseCoordinates[nx] !== 1
    || Array.from(streamwiseCoordinates).some((v, i) => !Number.isFinite(v) || i && !(v > streamwiseCoordinates[i - 1]))))
    throw new Error('Harmonic audit streamwise coordinates must increase strictly from zero to one and match the grid.');
  if (!boundaryConditions || typeof boundaryConditions !== 'object' || Array.isArray(boundaryConditions)
    || Object.keys(boundaryConditions).some(key => !['lower', 'upper'].includes(key))
    || Object.values(boundaryConditions).some(value => value !== undefined && !['fixed', 'giles-vertical', 'normal-curve'].includes(value)))
    throw new Error('Invalid harmonic audit boundary conditions.');
  const boundaries = Object.freeze(Object.fromEntries(['lower', 'upper'].map(side =>
    [side, boundaryConditions[side] === undefined ? 'fixed' : boundaryConditions[side]])));
  const lowerFree = boundaries.lower !== 'fixed', upperFree = boundaries.upper !== 'fixed';
  const mixed = lowerFree || upperFree;
  for (const [side, j] of [['lower', 0], ['upper', nt]]) if (boundaries[side] === 'giles-vertical'
    && nodes.some(row => row[j].y !== nodes[0][j].y))
    throw new Error(`Harmonic Giles ${side} boundary must be exactly horizontal.`);
  if (!boundaryCurves || typeof boundaryCurves !== 'object' || Array.isArray(boundaryCurves)
    || Object.keys(boundaryCurves).some(side => !['lower', 'upper'].includes(side) || boundaries[side] !== 'normal-curve'))
    throw new Error('Invalid harmonic audit boundary curves.');
  const curves = Object.fromEntries(['lower', 'upper'].filter(side => boundaries[side] === 'normal-curve')
    .map(side => {
      const spec = boundaryCurves[side];
      if (!spec || typeof spec !== 'object' || Array.isArray(spec))
        throw new Error(`Harmonic ${side} normal-curve boundary requires a graph descriptor.`);
      return [side, createNormalGraphBoundary(spec)];
    }));
  const curved = Object.keys(curves).length > 0;
  if (curved) {
    const origin = nodes[0][0]; let scale = 0;
    for (const row of nodes) for (const p of row) scale = Math.max(scale, Math.hypot(p.x - origin.x, p.y - origin.y));
    for (const [side, curve] of Object.entries(curves)) {
      const j = side === 'lower' ? 0 : nt;
      for (let i = 0; i <= nx; i++) {
        const p = nodes[i][j], expected = curve.evaluate(p.x).point;
        const tolerance = 64 * Number.EPSILON * Math.max(scale, Math.abs(p.y), Math.abs(expected.y));
        if (Math.abs(p.y - expected.y) > tolerance)
          throw new Error(`Harmonic ${side} boundary node does not lie on its prescribed curve.`);
        if (i && !(p.x > nodes[i - 1][j].x))
          throw new Error(`Harmonic ${side} curve stations must be strictly ordered.`);
      }
      for (const [i, p] of [[0, curve.points[0]], [nx, curve.points.at(-1)]]) {
        const x = nodes[i][j].x, tolerance = 64 * Number.EPSILON * Math.max(scale, Math.abs(x), Math.abs(p.x));
        if (Math.abs(x - p.x) > tolerance)
          throw new Error(`Harmonic ${side} curve endpoints must match the fixed end nodes.`);
      }
    }
  }
  const xi = Object.freeze(streamwiseCoordinates === undefined ? Array.from({ length: nx + 1 }, (_, i) => i / nx) : [...streamwiseCoordinates]);
  const ni = nx * refinement, nj = nt * refinement, n = (ni - 1) * (nj - 1);
  const xiRows = nj - 1 + Number(lowerFree) + Number(upperFree), nXi = (ni - 1) * xiRows;
  if (n > maxUnknowns) throw new Error(`Harmonic audit needs ${n} unknowns, exceeding its ${maxUnknowns} budget.`);
  if (nXi > maxUnknowns) throw new Error(`Harmonic xi audit needs ${nXi} unknowns, exceeding its ${maxUnknowns} budget.`);
  const total = massFlows.reduce((s, v) => s + v, 0), eta = [0];
  for (const m of massFlows) eta.push(eta.at(-1) + m / total);
  eta[nt] = 1;
  if (!Number.isFinite(total) || eta.some((e, j) => j && e <= eta[j - 1])) throw new Error('Unresolved harmonic audit mass coordinate.');
  const index = (i, j) => i > 0 && i < ni && j > 0 && j < nj ? (i - 1) * (nj - 1) + j - 1 : -1;
  // Natural homogeneous normal-Neumann data for xi on selected farfields.
  // Inlet/outlet and corners remain Dirichlet; eta retains all its boundary
  // mass labels. These are independent physical-space weak equations.
  const xiIndex = (i, j) => i > 0 && i < ni && (j > 0 || lowerFree) && (j < nj || upperFree)
    ? (i - 1) * xiRows + j - (lowerFree ? 0 : 1) : -1;
  const refined = Array.from({ length: ni + 1 }, (_, i) => Array.from({ length: nj + 1 }, (_, j) => {
    const a = Math.min(nx - 1, Math.floor(i / refinement)), b = Math.min(nt - 1, Math.floor(j / refinement));
    const s = i / refinement - a, t = j / refinement - b;
    const weights = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
    const corners = [nodes[a][b], nodes[a + 1][b], nodes[a + 1][b + 1], nodes[a][b + 1]];
    const point = { ...corners[0] };
    for (let k = 1; k < 4; k++) { point.x += weights[k] * (corners[k].x - corners[0].x); point.y += weights[k] * (corners[k].y - corners[0].y); }
    // Refine the supplied geometry and interpolate its seed labels. These
    // labels remain Dirichlet only on fixed sides; selected xi boundary
    // values are solved below. Uniform index labels would change stretched xi.
    const refinedXi = s === 0 ? xi[a] : s === 1 ? xi[a + 1] : xi[a] + s * (xi[a + 1] - xi[a]);
    return { ...point, xi: refinedXi, eta: (1 - t) * eta[b] + t * eta[b + 1], id: index(i, j),
      ...(mixed ? { xiId: xiIndex(i, j) } : {}) };
  }));
  if (curved) {
    for (let i = 0; i <= ni; i++) {
      // The graph agrees at every original station. Preserve these columns
      // exactly, including both fixed end lines and all observation nodes.
      // Between stations, spread each graph-versus-chord y discrepancy by
      // linear Coons weights. This refines geometry, not the prescribed data.
      if (i % refinement === 0) {
        for (let j = 0; j <= nj; j += refinement) {
          const p = nodes[i / refinement][j / refinement];
          refined[i][j].x = p.x; refined[i][j].y = p.y;
        }
        continue;
      }
      const lower = curves.lower ? curves.lower.evaluate(refined[i][0].x).point.y - refined[i][0].y : 0;
      const upper = curves.upper ? curves.upper.evaluate(refined[i][nj].x).point.y - refined[i][nj].y : 0;
      for (let j = 0; j <= nj; j++) refined[i][j].y += (1 - j / nj) * lower + (j / nj) * upper;
    }
  }
  const linear = {}, correction = {};
  if (!mixed) {
    // Keep the existing shared-matrix numerical path for fully prescribed
    // boundaries, including assembly and solve ordering.
    const rows = Array.from({ length: n }, () => new Map()), rhs = { xi: new Float64Array(n), eta: new Float64Array(n) };
    for (let i = 0; i < ni; i++) for (let j = 0; j < nj; j++) {
      const p = [refined[i][j], refined[i + 1][j], refined[i + 1][j + 1], refined[i][j + 1]], k = quadLaplaceMatrix(p);
      for (let a = 0; a < 4; a++) if (p[a].id >= 0) for (let b = 0; b < 4; b++) {
        const row = rows[p[a].id], coefficient = k[4 * a + b];
        if (p[b].id >= 0) row.set(p[b].id, (row.get(p[b].id) ?? 0) + coefficient);
        // Solve for the correction to the supplied labels. Difference form
        // preserves an exactly constant field without summing large terms.
        for (const key of ['xi', 'eta']) rhs[key][p[a].id] -= coefficient * (p[b][key] - p[a][key]);
      }
    }
    const matrix = sparseMatrix(rows.map(r => r.keys()));
    rows.forEach((row, i) => { for (const [j, value] of row) sparseAdd(matrix, i, j, value); });
    for (const key of ['xi', 'eta']) {
      const solved = solveSparseDirect(matrix, rhs[key]); correction[key] = solved.x;
      linear[key] = { relativeResidual: solved.relativeResidual, ordering: solved.ordering, refinements: solved.refinements };
    }
  } else {
    const rows = { xi: Array.from({ length: nXi }, () => new Map()), eta: Array.from({ length: n }, () => new Map()) };
    const rhs = { xi: new Float64Array(nXi), eta: new Float64Array(n) };
    for (let i = 0; i < ni; i++) for (let j = 0; j < nj; j++) {
      const p = [refined[i][j], refined[i + 1][j], refined[i + 1][j + 1], refined[i][j + 1]], k = quadLaplaceMatrix(p);
      for (const [key, id] of [['xi', 'xiId'], ['eta', 'id']]) for (let a = 0; a < 4; a++) if (p[a][id] >= 0) {
        const row = rows[key][p[a][id]];
        for (let b = 0; b < 4; b++) {
          const coefficient = k[4 * a + b];
          if (p[b][id] >= 0) row.set(p[b][id], (row.get(p[b][id]) ?? 0) + coefficient);
          rhs[key][p[a][id]] -= coefficient * (p[b][key] - p[a][key]);
        }
      }
      // The selected xi boundary has zero natural load. No boundary copy,
      // finite-difference derivative, or SLOR equation enters this assembly.
    }
    for (const key of ['xi', 'eta']) {
      const matrix = sparseMatrix(rows[key].map(r => r.keys()));
      rows[key].forEach((row, i) => { for (const [j, value] of row) sparseAdd(matrix, i, j, value); });
      const solved = solveSparseDirect(matrix, rhs[key]); correction[key] = solved.x;
      linear[key] = { relativeResidual: solved.relativeResidual, ordering: solved.ordering, refinements: solved.refinements };
    }
  }
  const valueAt = q => {
    const xiId = mixed ? q.xiId : q.id;
    return { xi: q.xi + (xiId < 0 ? 0 : correction.xi[xiId]), eta: q.eta + (q.id < 0 ? 0 : correction.eta[q.id]) };
  };
  const values = nodes.map((row, i) => row.map((p, j) => {
    const q = refined[i * refinement][j * refinement];
    return valueAt(q);
  }));
  const errors = values.map((row, i) => row.map((v, j) => {
    const width = Math.min(j ? eta[j] - eta[j - 1] : Infinity, j < nt ? eta[j + 1] - eta[j] : Infinity);
    const crosslineWidth = Math.min(i ? xi[i] - xi[i - 1] : Infinity, i < nx ? xi[i + 1] - xi[i] : Infinity);
    return { xi: v.xi - xi[i], eta: v.eta - eta[j], crosslineIntervals: (v.xi - xi[i]) / crosslineWidth, tubeIntervals: (v.eta - eta[j]) / width };
  }));
  return { refinement, unknowns: Math.max(n, nXi), unknownsByCoordinate: { xi: nXi, eta: n },
    boundaryConditions: boundaries,
    ...(curved ? { boundaryCurves: Object.fromEntries(Object.entries(curves).map(([side, curve]) => [side, curve.descriptor])),
      geometryRefinement: 'Q1 polygon approximations of fixed cubic Hermite graphs with transverse Coons corrections; fixed end lines and original nodes preserved. No exact curved-cell integration.' } : {}),
    streamwiseCoordinates: xi, linear, values, errors,
    ...(includeRefinedField ? { refinedField: refined.map(row => row.map(p => ({ x: p.x, y: p.y,
      ...valueAt(p) }))) } : {}),
    maximum: { crosslineIntervals: maximum(errors.flat().map(p => p.crosslineIntervals)), tubeIntervals: maximum(errors.flat().map(p => p.tubeIntervals)) } };
}

export function auditHarmonicGrid(region, { refinements = [1, 2, 4], coordinateLimit = .05, referenceLimit = .01, maxUnknowns = 50000 } = {}) {
  if (!Array.isArray(refinements) || refinements.length < 2 || refinements.some((n, k) => !Number.isInteger(n) || n < 1 || n > 8 || k && n !== 2 * refinements[k - 1])
    || !Number.isFinite(coordinateLimit) || !Number.isFinite(referenceLimit) || !(coordinateLimit > 0) || !(referenceLimit > 0)) throw new Error('Invalid harmonic audit refinement controls.');
  const levels = refinements.map(refinement => solveHarmonicGridReference(region, { refinement, maxUnknowns }));
  const changes = levels.slice(1).map((level, k) => {
    const previous = levels[k].errors.flat(), next = level.errors.flat();
    return Object.fromEntries(['crosslineIntervals', 'tubeIntervals'].map(key => [key, maximum(next.map((p, j) => p[key] - previous[j][key]))]));
  });
  const last = levels.at(-1), referenceChange = changes.at(-1);
  const converging = changes.length > 1 && Object.keys(referenceChange).every(key => referenceChange[key] < referenceLimit || referenceChange[key] < .7 * changes.at(-2)[key]);
  const resolved = converging && Object.values(referenceChange).every(e => e <= referenceLimit);
  const coordinatePass = Object.values(last.maximum).every(e => e <= coordinateLimit);
  const mixed = Object.values(last.boundaryConditions).some(mode => mode !== 'fixed');
  const curved = Object.values(last.boundaryConditions).includes('normal-curve');
  return { method: `Independent physical-space bilinear quadrilateral finite elements, 3x3 Gauss integration, ${curved ? 'Q1 polygon approximations of fixed cubic Hermite graph boundaries, transverse Coons geometry refinement; no exact curved-cell integration' : 'fixed polygonal boundary, nested reference refinement'}.${mixed ? ` Homogeneous normal-Neumann xi on selected ${curved ? 'farfields' : 'horizontal farfields'}; eta boundary mass labels remain prescribed.` : ''}`,
    scope: `${curved ? 'Harmonic coordinates on polygonal approximations to the prescribed fixed graphs, with mixed coordinate boundary conditions. Refinement reduces chord and FE error, not the graph approximation or fixed end-data error.' : mixed ? 'Harmonic coordinates on each fixed physical domain with the stated mixed coordinate boundary conditions.' : 'Harmonic coordinates inside each fixed-boundary region only.'} Does not certify dividing-cut placement, cross-cut velocity continuity, surface spacing or Euler/BL accuracy.`,
    limits: { coordinateIntervals: coordinateLimit, referenceChangeIntervals: referenceLimit },
    status: !coordinatePass ? 'harmonic coordinate discrepancy' : !resolved ? 'reference refinement unresolved' : 'sampled harmonic coordinate checks passed',
    passed: coordinatePass && resolved, referenceResolved: resolved, referenceChange, changes,
    levels: levels.map(({ values, errors, ...r }) => r), values: last.values, errors: last.errors };
}
