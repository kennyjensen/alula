// SPDX-License-Identifier: GPL-2.0-or-later
// Tensor Bernstein bows fixed relative to a moving Q1 nodal grid. Exact
// coefficient conformity preserves shared edges. Cell certificates bound
// local Jacobians and transverse guides, not global overlap or flow error.
import { certifyBernsteinCell, intervalPoint, intervalAdd, intervalSub, intervalMul, intervalDiv } from './bernstein-cell-certificate.js';

const registered = new WeakSet();
export const isPolynomialGridGeometry = geometry => registered.has(geometry);
const finite = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
const copy = grid => grid.map(row => row.map(p => ({ x: p.x, y: p.y })));
const equal = (a, b) => a.x === b.x && a.y === b.y;
const vector = (a, b, f = 1) => ({ x: f * (a.x - b.x), y: f * (a.y - b.y) });
const linear = (a, b, t) => ({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) });
const casteljau = (values, t) => {
  const work = values.map(p => ({ ...p }));
  for (let k = work.length - 1; k > 0; k--) for (let i = 0; i < k; i++) {
    work[i].x = (1 - t) * work[i].x + t * work[i + 1].x;
    work[i].y = (1 - t) * work[i].y + t * work[i + 1].y;
  }
  return work[0];
};
const tensor = (control, s, t) => casteljau(control.map(row => casteljau(row, t)), s);
const bilinear = (p, s, t) => {
  const shape = [(1 - s) * (1 - t), s * (1 - t), s * t, (1 - s) * t];
  const ds = [t - 1, 1 - t, t, -t], dt = [s - 1, -s, s, 1 - s];
  const result = { point: { ...p[0] }, ds: { x: 0, y: 0 }, dt: { x: 0, y: 0 } };
  for (let k = 1; k < 4; k++) for (const key of ['x', 'y']) {
    const d = p[k][key] - p[0][key];
    result.point[key] += shape[k] * d; result.ds[key] += ds[k] * d; result.dt[key] += dt[k] * d;
  }
  return result;
};

export function createPolynomialGridGeometry({ nodes, controlPoints }) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  const shape = grid => Array.isArray(grid) && grid.length === nx + 1 && grid.every(row => Array.isArray(row)
    && row.length === nt + 1 && row.every(finite));
  if (!(nx >= 1 && nt >= 1) || !shape(nodes) || !Array.isArray(controlPoints) || controlPoints.length !== nx
    || controlPoints.some(row => !Array.isArray(row) || row.length !== nt)) throw new Error('Invalid polynomial reference grid.');
  const ns = controlPoints[0][0]?.length - 1, mt = controlPoints[0][0]?.[0]?.length - 1;
  if (![ns, mt].every(d => Number.isInteger(d) && d >= 1 && d <= 3)
    || controlPoints.some(row => row.some(cell => !Array.isArray(cell) || cell.length !== ns + 1
      || cell.some(line => !Array.isArray(line) || line.length !== mt + 1 || line.some(p => !finite(p))))))
    throw new Error('Polynomial grid needs conforming tensor degrees between one and three.');
  const initial = copy(nodes), vertices = (grid, i, j) => [grid[i][j], grid[i + 1][j], grid[i + 1][j + 1], grid[i][j + 1]];
  const cells = controlPoints.map((row, i) => row.map((control, j) => {
    const p = vertices(initial, i, j);
    if (!equal(control[0][0], p[0]) || !equal(control[ns][0], p[1]) || !equal(control[ns][mt], p[2]) || !equal(control[0][mt], p[3]))
      throw new Error('Polynomial control corners must exactly match their observation nodes.');
    if (i && control[0].some((q, b) => !equal(q, controlPoints[i - 1][j][ns][b])))
      throw new Error('Polynomial cells do not share the same complete streamwise face.');
    if (j && control.some((line, a) => !equal(line[0], controlPoints[i][j - 1][a][mt])))
      throw new Error('Polynomial cells do not share the same complete transverse face.');
    // The stored double coefficients define the fixed bow polynomial. Its
    // corners are exactly zero; Q1 node motion never changes this polynomial.
    const bow = control.map((line, a) => line.map((q, b) => {
      if ((a === 0 || a === ns) && (b === 0 || b === mt)) return { x: 0, y: 0 };
      // Use the common edge's own chord on both neighboring cells. This
      // gives identical stored face coefficients even in floating point.
      if (a === 0) return vector(q, linear(p[0], p[3], b / mt));
      if (a === ns) return vector(q, linear(p[1], p[2], b / mt));
      if (b === 0) return vector(q, linear(p[0], p[1], a / ns));
      if (b === mt) return vector(q, linear(p[3], p[2], a / ns));
      return vector(q, bilinear(p, a / ns, b / mt).point);
    }));
    const ds = bow.slice(1).map((line, a) => line.map((q, b) => vector(q, bow[a][b], ns)));
    const dt = bow.map(line => line.slice(1).map((q, b) => vector(q, line[b], mt)));
    return { bow, ds, dt };
  }));
  const checkIndex = (i, j, s, t) => {
    if (!Number.isInteger(i) || !Number.isInteger(j) || i < 0 || i >= nx || j < 0 || j >= nt
      || !Number.isFinite(s) || !Number.isFinite(t) || s < 0 || s > 1 || t < 0 || t > 1)
      throw new Error('Polynomial reference coordinates are outside the grid.');
  };
  const correction = (i, j, s, t) => {
    checkIndex(i, j, s, t); const c = cells[i][j];
    return { point: tensor(c.bow, s, t), ds: tensor(c.ds, s, t), dt: tensor(c.dt, s, t) };
  };
  const onGrid = grid => {
    if (!shape(grid)) throw new Error('Invalid polynomial grid state.');
    const saved = copy(grid);
    return { at(i, j, s, t) {
      const c = correction(i, j, s, t), q = bilinear(vertices(saved, i, j), s, t);
      for (const key of ['point', 'ds', 'dt']) for (const axis of ['x', 'y']) q[key][axis] += c[key][axis];
      return q;
    } };
  };
  const quality = (grid, directions, controls = {}) => {
    if (!shape(grid) || !shape(directions)) throw new Error('Invalid polynomial grid certificate state.');
    const reports = [], invalidCells = [], unresolvedCells = [];
    let minimumJacobian = Infinity, minimumTransversality = Infinity;
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
      const p = vertices(grid, i, j), bow = cells[i][j].bow;
      // Bernstein coefficients of Q1 interpolation are its values at the
      // degree grid. Enclose the elevation and bow addition, not merely
      // the rounded result of those operations. Subtract a cell origin so
      // translation does not degrade derivative coefficient enclosures.
      const bounds = bow.map((line, a) => line.map((q, b) => {
        const s = intervalDiv(intervalPoint(a), intervalPoint(ns)), t = intervalDiv(intervalPoint(b), intervalPoint(mt));
        const oneMinusS = intervalSub(intervalPoint(1), s), oneMinusT = intervalSub(intervalPoint(1), t);
        const weights = [intervalMul(oneMinusS, oneMinusT), intervalMul(s, oneMinusT), intervalMul(s, t), intervalMul(oneMinusS, t)];
        return Object.fromEntries(['x', 'y'].map(axis => {
          let result = intervalPoint(q[axis]);
          for (let k = 1; k < 4; k++) result = intervalAdd(result,
            intervalMul(weights[k], intervalSub(intervalPoint(p[k][axis]), intervalPoint(p[0][axis]))));
          return [axis, result];
        }));
      }));
      const certificate = certifyBernsteinCell({ controlPoints: bounds,
        directions: [[directions[i][j], directions[i][j + 1]], [directions[i + 1][j], directions[i + 1][j + 1]]] },
      { minimumTransversality: 1e-10, ...controls });
      reports.push({ i, j, ...certificate });
      minimumJacobian = Math.min(minimumJacobian, certificate.lowerBounds.jacobian);
      minimumTransversality = Math.min(minimumTransversality, certificate.lowerBounds.normalizedTransversality);
      if (certificate.status === 'rejected') invalidCells.push({ i, j });
      else if (!certificate.valid) unresolvedCells.push({ i, j });
    }
    return { valid: reports.every(r => r.valid), positiveCells: reports.every(r => r.positive),
      transverse: reports.every(r => r.transverse), minimumJacobian, minimumTransversality, invalidCells, unresolvedCells,
      cellCertificates: reports, scope: 'Whole-cell local polynomial Jacobian and guide bounds. Shared faces conform; global boundary simplicity, overlap, clearance and physical accuracy are separate.' };
  };
  const exposedInitial = Object.freeze(initial.map(row => Object.freeze(row.map(p => Object.freeze({ ...p })))));
  const result = Object.freeze({ nx, nt, initial: exposedInitial, degree: Object.freeze({ s: ns, t: mt }),
    correction, ...onGrid(initial), onGrid, quality });
  registered.add(result);
  return result;
}
