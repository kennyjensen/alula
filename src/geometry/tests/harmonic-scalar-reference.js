// SPDX-License-Identifier: GPL-2.0-or-later
// Physical-space Q1/Q2 Dirichlet Laplace reference for arbitrary scalar labels.
// No inverse-grid metric, SLOR stencil or control-map derivative is reused.
import { quadLaplaceMatrix } from '../../numerics/quad-laplace.js';
import { quadraticQuadLaplaceMatrix } from '../../numerics/tests/quadratic-quad-laplace.js';
import { sparseMatrix, sparseAdd } from '../../numerics/sparse.js';
import { solveSparseDirect } from '../../numerics/klu.js';

export function solveHarmonicScalarReference({ nodes, labels }, { refinement = 1, degree = 1,
  quadratureOrder = degree === 2 ? 5 : 3, maxUnknowns = 50000 } = {}) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  if (!Array.isArray(nodes) || !(nx >= 1 && nt >= 1) || !Array.isArray(labels) || labels.length !== nx + 1
    || Array.from(nodes).some(row => !Array.isArray(row) || row.length !== nt + 1 || Array.from(row).some(p => !Number.isFinite(p?.x) || !Number.isFinite(p?.y)))
    || Array.from(labels).some(row => !Array.isArray(row) || row.length !== nt + 1 || !Array.from(row).every(Number.isFinite))
    || !Number.isInteger(refinement) || refinement < 1 || refinement > 8 || ![1, 2].includes(degree)
    || !(degree === 1 ? quadratureOrder === 3 : [3, 5].includes(quadratureOrder)) || !Number.isInteger(maxUnknowns) || maxUnknowns < 1)
    throw new Error('Invalid scalar harmonic reference input.');
  const stride = refinement * degree, ni = nx * stride, nj = nt * stride, n = (ni - 1) * (nj - 1);
  if (n > maxUnknowns) throw new Error(`Scalar reference needs ${n} unknowns, exceeding its ${maxUnknowns} budget.`);
  const index = (i, j) => !i || !j || i === ni || j === nj ? -1 : (i - 1) * (nj - 1) + j - 1;
  const refined = Array.from({ length: ni + 1 }, (_, i) => Array.from({ length: nj + 1 }, (_, j) => {
    const a = Math.min(nx - 1, Math.floor(i / stride)), b = Math.min(nt - 1, Math.floor(j / stride));
    const s = i / stride - a, t = j / stride - b, weights = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
    const corners = [[a, b], [a + 1, b], [a + 1, b + 1], [a, b + 1]];
    const point = { x: nodes[a][b].x, y: nodes[a][b].y, label: labels[a][b], id: index(i, j) };
    for (let k = 1; k < 4; k++) {
      const [u, v] = corners[k];
      point.x += weights[k] * (nodes[u][v].x - nodes[a][b].x);
      point.y += weights[k] * (nodes[u][v].y - nodes[a][b].y);
      point.label += weights[k] * (labels[u][v] - labels[a][b]);
    }
    return point;
  }));
  const rows = Array.from({ length: n }, () => new Map()), rhs = new Float64Array(n);
  for (let i = 0; i < ni; i += degree) for (let j = 0; j < nj; j += degree) {
    const corners = [refined[i][j], refined[i + degree][j], refined[i + degree][j + degree], refined[i][j + degree]];
    const p = degree === 1 ? corners : Array.from({ length: 9 }, (_, k) => refined[i + Math.floor(k / 3)][j + k % 3]);
    const k = degree === 1 ? quadLaplaceMatrix(corners) : quadraticQuadLaplaceMatrix(corners, { quadratureOrder });
    for (let a = 0; a < p.length; a++) if (p[a].id >= 0) for (let b = 0; b < p.length; b++) {
      const row = rows[p[a].id], value = k[p.length * a + b];
      if (p[b].id >= 0) row.set(p[b].id, (row.get(p[b].id) ?? 0) + value);
      // Interior labels are only a seed, not constraints. Solve A*delta=-A*seed
      // with zero boundary correction, in difference form for constant fields.
      rhs[p[a].id] -= value * (p[b].label - p[a].label);
    }
  }
  let solved = { x: new Float64Array(0), relativeResidual: 0 };
  if (n) {
    const matrix = sparseMatrix(rows.map(row => row.keys()));
    rows.forEach((row, i) => { for (const [j, value] of row) sparseAdd(matrix, i, j, value); });
    solved = solveSparseDirect(matrix, rhs);
  }
  let maximumError = 0;
  const values = nodes.map((row, i) => row.map((_, j) => {
    const id = index(i * stride, j * stride), correction = id < 0 ? 0 : solved.x[id];
    maximumError = Math.max(maximumError, Math.abs(correction));
    return labels[i][j] + correction;
  }));
  return { refinement, degree, quadratureOrder, unknowns: n, values, maximumError,
    linear: { relativeResidual: solved.relativeResidual, ordering: solved.ordering, refinements: solved.refinements },
    scope: `Scalar physical-space Q${degree} Laplace solution on the fixed boundary polygon. Boundary labels interpolate linearly along its edges for both degrees; interior labels are comparison targets only. No cross-cut flux continuity or curved-boundary accuracy certification.`,
    physicalAcceptance: false };
}
