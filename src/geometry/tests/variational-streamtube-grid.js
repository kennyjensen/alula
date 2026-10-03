// SPDX-License-Identifier: GPL-2.0-or-later
// Research Winslow functional: harmonic FORWARD coordinates on a prescribed,
// nonuniform potential/streamfunction mesh. Unlike inverse Laplace smoothing,
// this does not assume that the boundary correspondence is conformal.
// The browser continues to use the existing, separately audited method.
import { potentialGridQuality } from '../potential-plane-grid.js';
import { sparseMatrix, sparseIndex, sparseProduct } from '../../numerics/sparse.js';
import { solveStreamwiseBlockLines } from '../../numerics/tests/streamwise-block-lines.js';

const gauss = [.5 - Math.sqrt(3 / 5) / 2, .5, .5 + Math.sqrt(3 / 5) / 2], weights = [5 / 18, 4 / 9, 5 / 18];
const copy = grid => grid.map(row => row.map(p => ({ x: p.x, y: p.y })));

export function createVariationalStreamtubeGrid({ nodes, coordinates }) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  const shape = grid => Array.isArray(grid) && grid.length === nx + 1 && grid.every(row => Array.isArray(row)
    && row.length === nt + 1 && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)));
  if (!(nx >= 2 && nt >= 2) || !shape(nodes) || !shape(coordinates)) throw new Error('Invalid variational grid dimensions.');
  if (!potentialGridQuality(nodes).valid || !potentialGridQuality(coordinates).valid)
    throw new Error('The variational grid requires positive physical and computational cells.');
  const initial = copy(nodes), plane = copy(coordinates), origin = initial[0][0];
  let lengthScale = 0;
  for (const row of initial) for (const p of row) lengthScale = Math.max(lengthScale, Math.hypot(p.x - origin.x, p.y - origin.y));
  const n = 2 * (nx - 1) * (nt - 1), index = (i, j) => !i || i === nx || !j || j === nt ? -1 : 2 * ((i - 1) * (nt - 1) + j - 1);
  const elements = [], pattern = Array.from({ length: n }, () => new Set());
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ij = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], p = ij.map(([a, b]) => plane[a][b]);
    const ids = ij.flatMap(([a, b]) => { const k = index(a, b); return [k, k < 0 ? -1 : k + 1]; }), points = [];
    for (let q = 0; q < 3; q++) for (let r = 0; r < 3; r++) {
      const s = gauss[q], t = gauss[r], ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
      let us = 0, ut = 0, vs = 0, vt = 0;
      for (let k = 1; k < 4; k++) {
        us += ds[k] * (p[k].x - p[0].x); ut += dt[k] * (p[k].x - p[0].x);
        vs += ds[k] * (p[k].y - p[0].y); vt += dt[k] * (p[k].y - p[0].y);
      }
      const det = us * vt - ut * vs;
      if (!(det > 0) || !Number.isFinite(det)) throw new Error('Unresolved variational computational cell.');
      points.push({ weight: weights[q] * weights[r] * det,
        du: ds.map((v, k) => (v * vt - dt[k] * vs) / det), dv: ds.map((v, k) => (dt[k] * us - v * ut) / det) });
    }
    for (const a of ids) if (a >= 0) for (const b of ids) if (b >= 0) pattern[a].add(b);
    elements.push({ ij, ids, points });
  }
  const stencil = sparseMatrix(pattern);
  elements.forEach(e => { e.entries = e.ids.flatMap(a => e.ids.map(b => a < 0 || b < 0 ? -1 : sparseIndex(stencil, a, b))); });
  const validate = grid => {
    if (!shape(grid)) throw new Error('Invalid variational grid state.');
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) if (index(i, j) < 0
      && (grid[i][j].x !== initial[i][j].x || grid[i][j].y !== initial[i][j].y)) throw new Error('Variational grid boundaries must remain fixed.');
  };
  const evaluate = (grid, { linearize = false } = {}) => {
    validate(grid);
    let energy = 0;
    const gradient = linearize ? new Float64Array(n) : null;
    const matrix = linearize ? { ...stencil, values: new Float64Array(stencil.values.length) } : null;
    for (const { ij, ids, entries, points } of elements) {
      const p = ij.map(([a, b]) => grid[a][b]);
      for (const { weight, du, dv } of points) {
        let a = 0, b = 0, c = 0, d = 0;
        for (let k = 1; k < 4; k++) {
          const x = p[k].x - p[0].x, y = p[k].y - p[0].y;
          a += du[k] * x; b += dv[k] * x; c += du[k] * y; d += dv[k] * y;
        }
        const jacobian = a * d - b * c;
        if (!(jacobian > 0)) return { energy: Infinity };
        // E/2 minus the fixed computational area = 1/2 integral
        // ((x_u-y_v)^2+(x_v+y_u)^2)/J du dv. Its stationary point,
        // not zero conformal defect, is the target for general boundaries.
        const s1 = a - d, s2 = b + c, t = Math.sqrt(weight / jacobian), f1 = t * s1, f2 = t * s2;
        energy += .5 * (f1 * f1 + f2 * f2);
        if (linearize) {
          const df1 = [], df2 = [];
          for (let k = 0; k < 4; k++) {
            const jx = du[k] * d - dv[k] * c, jy = a * dv[k] - b * du[k];
            df1.push(t * (du[k] - .5 * s1 / jacobian * jx), t * (-dv[k] - .5 * s1 / jacobian * jy));
            df2.push(t * (dv[k] - .5 * s2 / jacobian * jx), t * (du[k] - .5 * s2 / jacobian * jy));
          }
          for (let k = 0; k < 8; k++) if (ids[k] >= 0) {
            gradient[ids[k]] += df1[k] * f1 + df2[k] * f2;
            for (let l = 0; l < 8; l++) if (ids[l] >= 0) matrix.values[entries[8 * k + l]] += df1[k] * df1[l] + df2[k] * df2[l];
          }
        }
      }
    }
    let residual = 0;
    if (linearize) for (let i = 0; i < n; i++) {
      let rowScale = 0;
      for (let p = matrix.rowPtr[i]; p < matrix.rowPtr[i + 1]; p++) rowScale += Math.abs(matrix.values[p]);
      residual = Math.max(residual, Math.abs(gradient[i]) / (rowScale * lengthScale));
    }
    if (!Number.isFinite(energy) || !Number.isFinite(residual)) throw new Error('Nonfinite variational grid functional.');
    return { energy, ...(linearize ? { gradient, matrix, residual } : {}) };
  };
  const move = (grid, delta, step) => grid.map((row, i) => row.map((p, j) => {
    const k = index(i, j); return k < 0 ? { ...p } : { x: p.x + step * delta[k], y: p.y + step * delta[k + 1] };
  }));
  return { nx, nt, n, initial, coordinates: plane, lengthScale, evaluate, move };
}

export function smoothVariationalStreamtubeGrid(system, { maxIterations = 40, tolerance = 1e-9,
  linearSolve, linearSolverName, onIteration } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !(tolerance > 0) || !Number.isFinite(tolerance))
    throw new Error('Invalid variational smoothing controls.');
  if (linearSolve !== undefined && typeof linearSolve !== 'function') throw new Error('Invalid variational linear solver.');
  const solve = linearSolve ?? ((matrix, rhs) => solveStreamwiseBlockLines(matrix, rhs, system));
  const linearBackend = linearSolve ? (linearSolverName ?? 'custom certified linear solver') : 'streamwise block SLOR';
  let nodes = copy(system.initial), reason = 'iteration limit'; const history = [];
  for (let iteration = 0; iteration <= maxIterations; iteration++) {
    const e = system.evaluate(nodes, { linearize: true });
    const linear = solve(e.matrix, Float64Array.from(e.gradient, v => -v));
    if (linear.converged === false || linear.x?.length !== system.n || !linear.x.every(Number.isFinite)) {
      reason = 'linear solve failed'; history.push({ iteration, energy: e.energy, residual: e.residual, linearResidual: linear.relativeResidual }); break;
    }
    const norm = a => Math.sqrt(a.reduce((s, v) => s + v * v, 0)), rhsNorm = norm(e.gradient);
    const defect = norm(sparseProduct(e.matrix, linear.x).map((v, k) => v + e.gradient[k]));
    const linearResidual = rhsNorm ? defect / rhsNorm : defect === 0 ? 0 : Infinity;
    if (!(linearResidual <= 1e-10)) {
      reason = 'linear solve failed'; history.push({ iteration, energy: e.energy, residual: e.residual, linearResidual }); break;
    }
    const maxUpdate = linear.x.reduce((m, v) => Math.max(m, Math.abs(v) / system.lengthScale), 0);
    const h = { iteration, energy: e.energy, residual: e.residual, maxUpdate, linearResidual, linearSweeps: linear.sweeps };
    history.push(h); onIteration?.(h, nodes);
    if (maxUpdate <= tolerance && e.residual <= tolerance) { reason = 'stationary'; break; }
    if (iteration === maxIterations) break;
    const descent = e.gradient.reduce((s, v, i) => s + v * linear.x[i], 0);
    if (!(descent < 0)) { reason = 'direction is not descent'; break; }
    let next, step = 1;
    // Standard feasible Armijo search. It scales a single variational step;
    // no node clipping, cell repair, boundary redistribution or extra penalty.
    for (; step >= 2 ** -24; step *= .5) {
      next = system.move(nodes, linear.x, step);
      if (potentialGridQuality(next).valid && system.evaluate(next).energy <= e.energy + 1e-4 * step * descent) break;
    }
    h.step = step;
    if (step < 2 ** -24) { reason = 'line search stalled'; break; }
    nodes = next;
  }
  return { nodes, reason, converged: reason === 'stationary', history, quality: potentialGridQuality(nodes), physicsValidated: false, linearBackend,
    formulation: `Q1 Winslow functional with prescribed computational coordinates and fixed physical boundaries; Gauss-Newton directions using ${linearBackend}. Research only.` };
}
