// SPDX-License-Identifier: GPL-2.0-or-later
// Research one-coordinate harmonic grid. Only streamfunction is harmonic;
// prescribed transverse guide lines fix the independent station placement.
// E = 1/2 integral_physical |grad(eta)|^2 dA. Q1 geometry and fixed mass labels.
// This is not a recovered MSET stencil and is not connected to the GUI.
import { potentialGridQuality } from './potential-plane-grid.js';
import { sparseMatrix, sparseIndex, sparseProduct } from '../numerics/sparse.js';
import { solveStreamwiseScalarLines } from '../numerics/streamwise-scalar-lines.js';

const gauss = [.5 - Math.sqrt(3 / 5) / 2, .5, .5 + Math.sqrt(3 / 5) / 2];
const weights = [5 / 18, 4 / 9, 5 / 18];
const copy = grid => grid.map(row => row.map(p => ({ x: p.x, y: p.y })));
const cross = (a, b) => a.x * b.y - a.y * b.x;
const minus = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const norm = v => Math.hypot(v.x, v.y);

export function createTransverseStreamfunctionGrid({ nodes, massFlows, directions }) {
  const nx = nodes?.length - 1, nt = massFlows?.length;
  const shape = grid => Array.isArray(grid) && grid.length === nx + 1 && grid.every(row => Array.isArray(row)
    && row.length === nt + 1 && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)));
  if (!(nx >= 2 && nt >= 2) || !shape(nodes) || !shape(directions) || !Array.isArray(massFlows)
    || !massFlows.every(m => Number.isFinite(m) && m > 0)) throw new Error('Invalid transverse streamfunction grid.');
  const totalMass = massFlows.reduce((s, m) => s + m, 0), eta = [0];
  for (const m of massFlows) eta.push(eta.at(-1) + m / totalMass);
  eta[nt] = 1;
  if (!Number.isFinite(totalMass) || eta.some((v, j) => j && !(v > eta[j - 1])))
    throw new Error('Unresolved transverse mass coordinates.');
  const initial = copy(nodes), guide = directions.map(row => row.map(p => {
    const length = norm(p);
    if (!(length > 0) || !Number.isFinite(length)) throw new Error('A transverse guide must be nonzero.');
    return { x: p.x / length, y: p.y / length };
  }));
  const origin = initial[0][0]; let lengthScale = 0;
  for (const row of initial) for (const p of row) lengthScale = Math.max(lengthScale, norm(minus(p, origin)));
  if (!(lengthScale > 0)) throw new Error('Degenerate transverse grid.');
  const n = (nx - 1) * (nt - 1), index = (i, j) => !i || i === nx || !j || j === nt ? -1 : (i - 1) * (nt - 1) + j - 1;
  const pattern = Array.from({ length: n }, () => new Set()), elements = [];
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ij = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], ids = ij.map(([a, b]) => index(a, b));
    const dir = ij.map(([a, b]) => guide[a][b]), points = [];
    for (let q = 0; q < 3; q++) for (let r = 0; r < 3; r++) {
      const s = gauss[q], t = gauss[r], ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
      points.push({ ds, dt, weight: weights[q] * weights[r] * (eta[j + 1] - eta[j]) ** 2,
        da: dir.map((d, k) => ({ x: ds[k] * d.x, y: ds[k] * d.y })),
        db: dir.map((d, k) => ({ x: dt[k] * d.x, y: dt[k] * d.y })) });
    }
    for (const a of ids) if (a >= 0) for (const b of ids) if (b >= 0) pattern[a].add(b);
    elements.push({ ij, ids, dir, points });
  }
  const stencil = sparseMatrix(pattern);
  for (const e of elements) e.entries = e.ids.flatMap(a => e.ids.map(b => a < 0 || b < 0 ? -1 : sparseIndex(stencil, a, b)));
  const validate = grid => {
    if (!shape(grid)) throw new Error('Invalid transverse grid state.');
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) {
      const d = minus(grid[i][j], initial[i][j]);
      if (index(i, j) < 0 && (d.x !== 0 || d.y !== 0)) throw new Error('Transverse grid boundaries must remain fixed.');
      if (Math.abs(cross(d, guide[i][j])) > 64 * Number.EPSILON * lengthScale)
        throw new Error('A transverse grid node left its prescribed guide line.');
    }
  };
  const quality = grid => {
    validate(grid);
    const geometry = potentialGridQuality(grid); let minimumTransversality = Infinity;
    // r_s is affine in t; the interpolated guide D(s,t) is bilinear.
    // cross(r_s,D) is affine in s and quadratic in t. Its exact minimum
    // occurs on s=0 or 1, at a t endpoint or an interior quadratic vertex.
    // A positive lower bound excludes a tangent/vanishing guide everywhere,
    // not just at quadrature points. Corner norms bound the denominator.
    for (const { ij, dir } of elements) {
      const p = ij.map(([i, j]) => grid[i][j]), a0 = minus(p[1], p[0]), a1 = minus(minus(p[2], p[3]), a0);
      for (const [lo, hi] of [[0, 3], [1, 2]]) {
        const d0 = dir[lo], d1 = minus(dir[hi], d0);
        const c0 = cross(a0, d0), c1 = cross(a0, d1) + cross(a1, d0), c2 = cross(a1, d1);
        let minimum = Math.min(c0, c0 + c1 + c2);
        const t = -c1 / (2 * c2);
        if (c2 > 0 && t > 0 && t < 1) minimum = Math.min(minimum, c0 + t * (c1 + t * c2));
        const bound = Math.max(norm(a0), norm({ x: a0.x + a1.x, y: a0.y + a1.y })) * Math.max(norm(d0), norm(dir[hi]));
        minimumTransversality = Math.min(minimumTransversality, bound > 0 ? minimum / bound : -Infinity);
      }
    }
    const transverse = minimumTransversality > 1e-10;
    return { ...geometry, positiveCells: geometry.valid, transverse, minimumTransversality, valid: geometry.valid && transverse };
  };
  if (!quality(initial).valid) throw new Error('Transverse initialization requires positive cells and transverse guides.');
  const evaluate = (grid, { linearize = false, hessian = 'exact' } = {}) => {
    validate(grid);
    if (!['exact', 'gauss-newton'].includes(hessian)) throw new Error('Unknown transverse energy Hessian.');
    const gradient = linearize ? new Float64Array(n) : null;
    const matrix = linearize ? { ...stencil, values: new Float64Array(stencil.values.length) } : null;
    let energy = 0;
    for (const { ij, ids, entries, points } of elements) {
      const p = ij.map(([i, j]) => grid[i][j]);
      for (const { ds, dt, da, db, weight } of points) {
        const a = { x: 0, y: 0 }, b = { x: 0, y: 0 };
        for (let k = 1; k < 4; k++) {
          const d = minus(p[k], p[0]);
          a.x += ds[k] * d.x; a.y += ds[k] * d.y;
          b.x += dt[k] * d.x; b.y += dt[k] * d.y;
        }
        const jacobian = cross(a, b), numerator = a.x * a.x + a.y * a.y;
        if (!(jacobian > 0)) return { energy: Infinity };
        energy += .5 * weight * numerator / jacobian;
        if (!linearize) continue;
        const dJ = da.map((d, k) => cross(d, b) + cross(a, db[k]));
        const halfDQ = da.map(d => a.x * d.x + a.y * d.y);
        for (let k = 0; k < 4; k++) if (ids[k] >= 0) {
          gradient[ids[k]] += weight * (halfDQ[k] / jacobian - .5 * numerator * dJ[k] / jacobian ** 2);
          for (let l = 0; l < 4; l++) if (ids[l] >= 0) {
            const ddJ = cross(da[k], db[l]) + cross(da[l], db[k]);
            const daDot = da[k].x * da[l].x + da[k].y * da[l].y;
            const h = hessian === 'exact'
              ? daDot / jacobian - (halfDQ[k] * dJ[l] + halfDQ[l] * dJ[k] + .5 * numerator * ddJ) / jacobian ** 2
                + numerator * dJ[k] * dJ[l] / jacobian ** 3
              : daDot / jacobian - .5 * (halfDQ[k] * dJ[l] + halfDQ[l] * dJ[k]) / jacobian ** 2
                + .25 * numerator * dJ[k] * dJ[l] / jacobian ** 3;
            matrix.values[entries[4 * k + l]] += weight * h;
          }
        }
      }
    }
    let residual = 0;
    if (linearize) for (let i = 0; i < n; i++) {
      let rowScale = 0;
      for (let k = matrix.rowPtr[i]; k < matrix.rowPtr[i + 1]; k++) rowScale += Math.abs(matrix.values[k]);
      residual = Math.max(residual, Math.abs(gradient[i]) / (rowScale * lengthScale));
    }
    if (!Number.isFinite(energy) || !Number.isFinite(residual)) throw new Error('Nonfinite transverse streamfunction energy.');
    return { energy, ...(linearize ? { gradient, matrix, residual } : {}) };
  };
  const move = (grid, delta, step = 1) => {
    validate(grid);
    if (delta?.length !== n || !Array.from(delta).every(Number.isFinite) || !Number.isFinite(step)) throw new Error('Invalid transverse grid direction.');
    return grid.map((row, i) => row.map((p, j) => {
      const k = index(i, j), dir = guide[i][j];
      return k < 0 ? { ...p } : { x: p.x + step * delta[k] * dir.x, y: p.y + step * delta[k] * dir.y };
    }));
  };
  return { nx, nt, n, initial, eta, directions: guide, lengthScale, evaluate, quality, move, validate };
}

export function smoothTransverseStreamfunctionGrid(system, { maxIterations = 40, tolerance = 1e-9,
  linearSolve, linearSolverName, onIteration } = {}) {
  if (!Number.isInteger(maxIterations) || maxIterations < 0 || !Number.isFinite(tolerance) || !(tolerance > 0)
    || (linearSolve !== undefined && typeof linearSolve !== 'function')) throw new Error('Invalid transverse smoothing controls.');
  const solve = linearSolve ?? ((matrix, rhs) => solveStreamwiseScalarLines(matrix, rhs, system));
  const linearBackend = linearSolve ? (linearSolverName ?? 'custom certified linear solver') : 'scalar streamwise SLOR';
  let nodes = copy(system.initial), reason = 'iteration limit'; const history = [];
  const vectorNorm = a => Math.sqrt(a.reduce((s, v) => s + v * v, 0));
  for (let iteration = 0; iteration <= maxIterations; iteration++) {
    const e = system.evaluate(nodes, { linearize: true });
    const linear = solve(e.matrix, Float64Array.from(e.gradient, v => -v));
    if (linear.converged === false || linear.x?.length !== system.n || !linear.x.every(Number.isFinite)) {
      reason = 'linear solve failed'; history.push({ iteration, energy: e.energy, residual: e.residual }); break;
    }
    const rhsNorm = vectorNorm(e.gradient), defect = vectorNorm(sparseProduct(e.matrix, linear.x).map((v, k) => v + e.gradient[k]));
    const linearResidual = rhsNorm ? defect / rhsNorm : defect === 0 ? 0 : Infinity;
    if (!(linearResidual <= 1e-10)) { reason = 'linear solve failed'; history.push({ iteration, energy: e.energy, residual: e.residual, linearResidual }); break; }
    const maxUpdate = linear.x.reduce((m, v) => Math.max(m, Math.abs(v) / system.lengthScale), 0);
    const h = { iteration, energy: e.energy, residual: e.residual, maxUpdate, linearResidual, linearSweeps: linear.sweeps };
    history.push(h); onIteration?.(h, nodes);
    if (maxUpdate <= tolerance && e.residual <= tolerance) { reason = 'stationary'; break; }
    if (iteration === maxIterations) break;
    const descent = e.gradient.reduce((s, v, k) => s + v * linear.x[k], 0);
    if (!(descent < 0)) { reason = 'Newton direction is not descent'; break; }
    let step = 1, next;
    for (; step >= 2 ** -24; step *= .5) {
      next = system.move(nodes, linear.x, step);
      if (system.quality(next).valid && system.evaluate(next).energy <= e.energy + 1e-4 * step * descent) break;
    }
    h.step = step;
    if (step < 2 ** -24) { reason = 'line search stalled'; break; }
    nodes = next;
  }
  return { nodes, converged: reason === 'stationary', reason, history, quality: system.quality(nodes), physicsValidated: false,
    linearBackend, formulation: 'Q1 streamfunction Dirichlet energy; fixed boundaries and transverse guide lines; exact Newton Hessian. Research only.' };
}
