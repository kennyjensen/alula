// SPDX-License-Identifier: GPL-2.0-or-later
// Independent scalar Q1 field reference integrated on a supplied fixed curved
// geometry map. Reference subdivision retains that map and its boundary data.
import { sparseMatrix, sparseAdd } from '../../numerics/sparse.js';
import { solveSparseDirect } from '../../numerics/klu.js';
import { gaussUnitRule } from '../../numerics/gauss-unit.js';

const finite = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);

export function solveCurvedHarmonicReference({ nodes, massFlows, geometry }, { refinement = 1, maxUnknowns = 50000, quadratureOrder = 3 } = {}) {
  const { nodes: gauss, weights } = gaussUnitRule(quadratureOrder);
  const nx = nodes?.length - 1, nt = massFlows?.length;
  if (!(nx >= 1 && nt >= 1) || !Array.isArray(nodes) || !Array.isArray(massFlows)
    || nodes.some(row => !Array.isArray(row) || row.length !== nt + 1 || row.some(p => !finite(p)))
    || !massFlows.every(m => Number.isFinite(m) && m > 0) || typeof geometry?.at !== 'function'
    || !Number.isInteger(refinement) || refinement < 1 || refinement > 8
    || !Number.isInteger(maxUnknowns) || maxUnknowns < 1) throw new Error('Invalid curved harmonic reference.');
  const ni = nx * refinement, nj = nt * refinement, n = (ni - 1) * (nj - 1);
  if (n > maxUnknowns) throw new Error(`Curved reference needs ${n} unknowns, exceeding its ${maxUnknowns} budget.`);
  const totalMass = massFlows.reduce((s, m) => s + m, 0), eta = [0];
  for (const m of massFlows) eta.push(eta.at(-1) + m / totalMass);
  eta[nt] = 1;
  if (!Number.isFinite(totalMass) || eta.some((v, j) => j && !(v > eta[j - 1]))) throw new Error('Unresolved curved-reference mass coordinates.');
  let lengthScale = 0;
  for (const row of nodes) for (const p of row) lengthScale = Math.max(lengthScale, Math.hypot(p.x - nodes[0][0].x, p.y - nodes[0][0].y));
  if (!(lengthScale > 0) || !Number.isFinite(lengthScale)) throw new Error('Degenerate curved-reference geometry.');
  // The returned values are claimed at these original nodes, so a supplied
  // geometry cannot silently move those observation points.
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) for (const [s, t] of [[0, 0], [1, 0], [1, 1], [0, 1]]) {
    const p = geometry.at(i, j, s, t).point, q = nodes[i + s][j + t];
    if (!finite(p) || Math.hypot(p.x - q.x, p.y - q.y) > 128 * Number.EPSILON * lengthScale)
      throw new Error('Curved geometry does not match its observation nodes.');
  }
  const index = (i, j) => !i || !j || i === ni || j === nj ? -1 : (i - 1) * (nj - 1) + j - 1;
  const rows = Array.from({ length: n }, () => new Map()), rhs = new Float64Array(n);
  let minimumQuadratureJacobian = Infinity;
  for (let i = 0; i < ni; i++) for (let j = 0; j < nj; j++) {
    const a = Math.floor(i / refinement), b = Math.floor(j / refinement);
    const localI = i - a * refinement, localJ = j - b * refinement, deltaEta = (eta[b + 1] - eta[b]) / refinement;
    // Only cell-local differences enter the RHS. Removing the constant
    // offset avoids subtracting nearly equal labels close to eta=1.
    const labels = [0, 0, deltaEta, deltaEta];
    const ids = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]].map(([x, y]) => index(x, y));
    const k = new Float64Array(16);
    for (let q = 0; q < gauss.length; q++) for (let r = 0; r < gauss.length; r++) {
      const s = gauss[q], t = gauss[r], value = geometry.at(a, b, (localI + s) / refinement, (localJ + t) / refinement);
      if (!finite(value?.point) || !finite(value?.ds) || !finite(value?.dt)) throw new Error('Nonfinite curved geometry sample.');
      // Derivatives of the parent map with respect to a refined subcell.
      const xs = value.ds.x / refinement, ys = value.ds.y / refinement;
      const xt = value.dt.x / refinement, yt = value.dt.y / refinement, jacobian = xs * yt - ys * xt;
      minimumQuadratureJacobian = Math.min(minimumQuadratureJacobian, jacobian * refinement ** 2);
      if (!(jacobian > 0) || !Number.isFinite(jacobian)) throw new Error('Nonpositive curved-reference quadrature Jacobian.');
      const ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
      const gx = ds.map((v, l) => (yt * v - ys * dt[l]) / jacobian);
      const gy = ds.map((v, l) => (-xt * v + xs * dt[l]) / jacobian);
      const weight = weights[q] * weights[r] * jacobian;
      for (let x = 0; x < 4; x++) for (let y = 0; y < 4; y++) k[4 * x + y] += weight * (gx[x] * gx[y] + gy[x] * gy[y]);
    }
    for (let x = 0; x < 4; x++) if (ids[x] >= 0) for (let y = 0; y < 4; y++) {
      const coefficient = k[4 * x + y], row = rows[ids[x]];
      if (ids[y] >= 0) row.set(ids[y], (row.get(ids[y]) ?? 0) + coefficient);
      rhs[ids[x]] -= coefficient * (labels[y] - labels[x]);
    }
  }
  let solved = { x: new Float64Array(0), relativeResidual: 0 };
  if (n) {
    const matrix = sparseMatrix(rows.map(row => row.keys()));
    rows.forEach((row, i) => { for (const [j, value] of row) sparseAdd(matrix, i, j, value); });
    solved = solveSparseDirect(matrix, rhs);
  }
  let maximumTubeIntervals = 0;
  const values = nodes.map((row, i) => row.map((p, j) => {
    const id = index(i * refinement, j * refinement), correction = id < 0 ? 0 : solved.x[id];
    const width = Math.min(j ? eta[j] - eta[j - 1] : Infinity, j < nt ? eta[j + 1] - eta[j] : Infinity);
    maximumTubeIntervals = Math.max(maximumTubeIntervals, Math.abs(correction) / width);
    return eta[j] + correction;
  }));
  return { refinement, quadratureOrder, unknowns: n, values, maximumTubeIntervals,
    linear: { relativeResidual: solved.relativeResidual, ordering: solved.ordering, refinements: solved.refinements },
    geometryCheck: { minimumQuadratureJacobian, globallyCertified: false },
    scope: 'Q1 scalar Laplace reference on the supplied curved map with prescribed boundary mass labels. Checks quadrature Jacobians and observation nodes; caller must establish global geometry validity and physical boundary-data accuracy.',
    physicalAcceptance: false };
}
