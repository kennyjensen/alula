// SPDX-License-Identifier: GPL-2.0-or-later
// Research inverse-potential initialization. Where w=phi+i*psi is conformal,
// x(w) and y(w) are harmonic. The supplied panel coordinates define the
// computational quadrilaterals; physical boundary nodes are held fixed.
// Constant Q1 weak-Laplace coefficients permit exact tridiagonal block SOR.
// This verifies a linear grid problem, not panel/conformal/flow accuracy.
import { quadLaplaceMatrix } from '../numerics/quad-laplace.js';

const copy = nodes => nodes.map(row => row.map(p => ({ x: p.x, y: p.y })));
export function potentialGridQuality(nodes) {
  let minCornerSine = Infinity; const invalidCells = [];
  for (let i = 0; i < nodes.length - 1; i++) for (let j = 0; j < nodes[0].length - 1; j++) {
    const p = [nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1]]; let valid = true;
    for (let k = 0; k < 4; k++) {
      const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4];
      const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
      const sine = (ux * vy - uy * vx) / (Math.hypot(ux, uy) * Math.hypot(vx, vy));
      minCornerSine = Math.min(minCornerSine, Number.isFinite(sine) ? sine : -Infinity);
      if (!(sine > 1e-12)) valid = false;
    }
    if (!valid) invalidCells.push({ i, j });
  }
  return { valid: invalidCells.length === 0, minCornerSine, invalidCells };
}

export function createPotentialPlaneGrid({ nodes, coordinates }) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  const shape = grid => Array.isArray(grid) && grid.length === nx + 1
    && grid.every(row => Array.isArray(row) && row.length === nt + 1
      && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)));
  if (!(nx >= 2 && nt >= 2) || !shape(nodes) || !shape(coordinates)) throw new Error('Invalid potential-plane grid.');
  const initial = copy(nodes), plane = copy(coordinates), origin = initial[0][0];
  let lengthScale = 0;
  for (const row of initial) for (const p of row) lengthScale = Math.max(lengthScale, Math.hypot(p.x - origin.x, p.y - origin.y));
  if (!(lengthScale > 0) || !Number.isFinite(lengthScale)) throw new Error('Degenerate physical grid.');
  const id = (i, j) => i * (nt + 1) + j, point = index => [Math.floor(index / (nt + 1)), index % (nt + 1)];
  const boundary = (i, j) => !i || i === nx || !j || j === nt;
  const coefficients = Array.from({ length: (nx + 1) * (nt + 1) }, () => new Map());
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ij = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], matrix = quadLaplaceMatrix(ij.map(([a, b]) => plane[a][b]));
    for (let a = 0; a < 4; a++) if (!boundary(...ij[a])) {
      const row = coefficients[id(...ij[a])];
      for (let b = 0; b < 4; b++) if (a !== b) {
        const index = id(...ij[b]); row.set(index, (row.get(index) ?? 0) + matrix[4 * a + b]);
      }
    }
  }
  const rows = [];
  for (let i = 1; i < nx; i++) for (let j = 1; j < nt; j++) {
    const coefficientsAtNode = coefficients[id(i, j)], diagonal = -[...coefficientsAtNode.values()].reduce((s, v) => s + v, 0);
    const scale = Math.abs(diagonal) + [...coefficientsAtNode.values()].reduce((s, v) => s + Math.abs(v), 0);
    if (!(diagonal > 0) || !Number.isFinite(scale)) throw new Error('Invalid potential-plane Laplace row.');
    rows.push({ i, j, diagonal, scale, entries: [...coefficientsAtNode].map(([index, value]) => ({ ij: point(index), value })) });
  }
  const rowAt = (i, j) => rows[(i - 1) * (nt - 1) + j - 1];
  // Factor each streamwise line once; geometry and coefficients do not
  // change during relaxation. The factors serve both physical coordinates.
  const lines = Array.from({ length: nt - 1 }, (_, index) => {
    const j = index + 1, diagonal = new Float64Array(nx - 1), upper = new Float64Array(nx - 1), factors = new Float64Array(nx - 1);
    for (let i = 1; i < nx; i++) {
      const k = i - 1, row = rowAt(i, j); diagonal[k] = row.diagonal;
      upper[k] = coefficients[id(i, j)].get(id(i + 1, j)) ?? 0;
      if (k) {
        factors[k] = (coefficients[id(i, j)].get(id(i - 1, j)) ?? 0) / diagonal[k - 1];
        diagonal[k] -= factors[k] * upper[k - 1];
      }
      if (!(diagonal[k] > 32 * Number.EPSILON * row.scale)) throw new Error('Unresolved potential-plane SLOR line.');
    }
    return { j, diagonal, upper, factors };
  });
  const validate = grid => {
    if (!shape(grid)) throw new Error('Invalid physical potential-grid state.');
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) if (boundary(i, j)
      && (grid[i][j].x !== initial[i][j].x || grid[i][j].y !== initial[i][j].y)) throw new Error('Potential-grid boundary coordinates must remain fixed.');
  };
  const residualAt = (grid, row, key) => {
    const center = grid[row.i][row.j][key];
    return row.entries.reduce((sum, e) => sum + e.value * (grid[e.ij[0]][e.ij[1]][key] - center), 0);
  };
  const residuals = grid => {
    validate(grid); let residual = 0;
    const values = rows.map(row => {
      const result = { i: row.i, j: row.j };
      for (const key of ['x', 'y']) {
        result[key] = residualAt(grid, row, key);
        residual = Math.max(residual, Math.abs(result[key]) / (row.scale * lengthScale));
      }
      return result;
    });
    if (!Number.isFinite(residual)) throw new Error('Nonfinite potential-plane residual.');
    return { residual, rows: values };
  };
  const sweep = (grid, omega = 1.3) => {
    validate(grid);
    if (!Number.isFinite(omega) || !(omega > 0 && omega < 2)) throw new Error('Potential SLOR omega must lie between zero and two.');
    const next = copy(grid); let maxUpdate = 0;
    for (const { j, diagonal, upper, factors } of lines) for (const key of ['x', 'y']) {
      const rhs = Float64Array.from({ length: nx - 1 }, (_, k) => -residualAt(next, rowAt(k + 1, j), key));
      for (let k = 1; k < rhs.length; k++) rhs[k] -= factors[k] * rhs[k - 1];
      for (let k = rhs.length - 1; k >= 0; k--) {
        rhs[k] = (rhs[k] - (k + 1 < rhs.length ? upper[k] * rhs[k + 1] : 0)) / diagonal[k];
        const delta = omega * rhs[k]; next[k + 1][j][key] += delta;
        maxUpdate = Math.max(maxUpdate, Math.abs(delta) / lengthScale);
      }
    }
    return { nodes: next, maxUpdate };
  };
  return { initial, coordinates: plane, nx, nt, lengthScale, residuals, sweep };
}

export function smoothPotentialPlaneGrid(system, { maxSweeps = 600, tolerance = 1e-11, omega = 1.3, onSweep } = {}) {
  if (!Number.isInteger(maxSweeps) || maxSweeps < 0 || !Number.isFinite(tolerance) || tolerance <= 0
    || !Number.isFinite(omega) || !(omega > 0 && omega < 2)) throw new Error('Invalid potential-plane relaxation controls.');
  let nodes = copy(system.initial), reason = 'sweep limit'; const history = [];
  for (let iteration = 0; iteration <= maxSweeps; iteration++) {
    const { residual } = system.residuals(nodes);
    history.push({ iteration, residual }); onSweep?.(history.at(-1), nodes);
    if (residual <= tolerance) { reason = 'converged'; break; }
    if (iteration === maxSweeps) break;
    const next = system.sweep(nodes, omega); nodes = next.nodes; history.at(-1).maxUpdate = next.maxUpdate;
    if (!next.maxUpdate) { reason = 'roundoff stagnation'; break; }
  }
  const quality = potentialGridQuality(nodes);
  if (reason === 'converged' && !quality.valid) reason = 'folded converged grid';
  return { nodes, history, quality, converged: reason === 'converged', reason, omega, tolerance, physicsValidated: false,
    formulation: 'Inverse potential/streamfunction map; fixed Q1 Laplace matrix, factored streamwise SLOR lines and fixed physical boundaries. Conformal and physical-grid accuracy require separate validation.' };
}
