// SPDX-License-Identifier: GPL-2.0-or-later
// Inverse physical-space Q1 Laplace residual with fixed nodal eta labels.
// Geometry moves along prescribed scalar guide directions. This is R=K(z)eta,
// not stationarity of the discrete mesh Dirichlet energy. No nonlinear solve.
// An optional fixed curved correction is sampled once at construction. The
// curved path checks quadrature Jacobians only; it is not a cell certificate.
import { sparseMatrix, sparseIndex } from '../numerics/sparse.js';
import { gaussUnitRule } from '../numerics/gauss-unit.js';

const dot = (a, b) => a.x * b.x + a.y * b.y;

export function createTransverseHarmonicResidual({ nodes, massFlows, directions, geometryCorrection, quadratureOrder = 3 }) {
  const { nodes: gauss, weights } = gaussUnitRule(quadratureOrder);
  const nx = nodes?.length - 1, nt = massFlows?.length;
  const validShape = grid => Array.isArray(grid) && grid.length === nx + 1 && grid.every(row => Array.isArray(row)
    && row.length === nt + 1 && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)));
  if (!(nx >= 2 && nt >= 2) || !validShape(nodes) || !validShape(directions) || !Array.isArray(massFlows)
    || !massFlows.every(m => Number.isFinite(m) && m > 0)
    || (geometryCorrection !== undefined && typeof geometryCorrection !== 'function')) throw new Error('Invalid transverse harmonic residual grid.');
  const totalMass = massFlows.reduce((sum, m) => sum + m, 0), eta = [0];
  for (const m of massFlows) eta.push(eta.at(-1) + m / totalMass);
  eta[nt] = 1;
  if (!Number.isFinite(totalMass) || eta.some((e, j) => j && !(e > eta[j - 1])))
    throw new Error('Unresolved transverse harmonic mass coordinates.');
  const guide = directions.map(row => row.map(d => {
    const length = Math.hypot(d.x, d.y);
    if (!(length > 0) || !Number.isFinite(length)) throw new Error('A transverse harmonic guide must be nonzero.');
    return { x: d.x / length, y: d.y / length };
  }));
  const n = (nx - 1) * (nt - 1);
  const index = (i, j) => !i || i === nx || !j || j === nt ? -1 : (i - 1) * (nt - 1) + j - 1;
  const quadrature = gauss.flatMap((s, q) => gauss.map((t, r) => ({
    s, t, ds: [t - 1, 1 - t, t, -t], dt: [s - 1, -s, s, 1 - s], weight: weights[q] * weights[r],
  })));
  const finite = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
  const curvedSample = (i, j, s, t) => {
    const value = geometryCorrection(i, j, s, t);
    if (!finite(value?.point) || !finite(value?.ds) || !finite(value?.dt)) throw new Error('Nonfinite fixed harmonic geometry correction.');
    return { point: { ...value.point }, ds: { ...value.ds }, dt: { ...value.dt } };
  };
  let lengthScale = 0;
  for (const row of nodes) for (const p of row) lengthScale = Math.max(lengthScale, Math.hypot(p.x - nodes[0][0].x, p.y - nodes[0][0].y));
  if (!(lengthScale > 0) || !Number.isFinite(lengthScale)) throw new Error('Degenerate transverse harmonic reference extent.');
  const pattern = Array.from({ length: n }, () => new Set()), elements = [];
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ij = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]], ids = ij.map(([a, b]) => index(a, b));
    const dir = ij.map(([a, b]) => guide[a][b]), labels = [eta[j], eta[j], eta[j + 1], eta[j + 1]];
    for (const a of ids) if (a >= 0) for (const b of ids) if (b >= 0) pattern[a].add(b);
    let corrections;
    if (geometryCorrection) {
      for (const [s, t] of [[0, 0], [1, 0], [1, 1], [0, 1]]) {
        const { point } = curvedSample(i, j, s, t);
        if (Math.hypot(point.x, point.y) > 128 * Number.EPSILON * lengthScale)
          throw new Error('A fixed harmonic geometry correction moved an observation node.');
      }
      corrections = quadrature.map(({ s, t }) => curvedSample(i, j, s, t));
    }
    elements.push({ ij, ids, dir, labels, corrections });
  }
  const stencil = sparseMatrix(pattern);
  for (const e of elements) e.entries = e.ids.flatMap(a => e.ids.map(b => a < 0 || b < 0 ? -1 : sparseIndex(stencil, a, b)));

  const evaluate = (grid, { linearize = true } = {}) => {
    if (!validShape(grid) || typeof linearize !== 'boolean') throw new Error('Invalid transverse harmonic residual state.');
    const residual = new Float64Array(n);
    const matrix = linearize ? { ...stencil, values: new Float64Array(stencil.values.length) } : null;
    for (const { ij, ids, dir, labels, entries, corrections } of elements) {
      const vertices = ij.map(([i, j]) => grid[i][j]), origin = vertices[0];
      const p = vertices.map(v => ({ x: v.x - origin.x, y: v.y - origin.y }));
      for (let k = 0; !corrections && k < 4; k++) {
        const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4];
        if (!((b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x) > 0))
          throw new Error('Invalid quadrilateral in transverse harmonic residual.');
      }
      for (let q = 0; q < quadrature.length; q++) {
        const { ds, dt, weight } = quadrature[q], c = corrections?.[q];
        let xs = c?.ds.x ?? 0, ys = c?.ds.y ?? 0, xt = c?.dt.x ?? 0, yt = c?.dt.y ?? 0;
        for (let k = 1; k < 4; k++) {
          xs += ds[k] * p[k].x; ys += ds[k] * p[k].y;
          xt += dt[k] * p[k].x; yt += dt[k] * p[k].y;
        }
        const jacobian = xs * yt - ys * xt;
        if (!(jacobian > 0) || !Number.isFinite(jacobian)) throw new Error('Invalid transverse harmonic Jacobian.');
        const gradients = ds.map((s, k) => ({ x: (yt * s - ys * dt[k]) / jacobian, y: (-xt * s + xs * dt[k]) / jacobian }));
        const flow = { x: 0, y: 0 };
        // Difference form is invariant under an additive label offset and
        // avoids cancellation when a tube has a small mass fraction.
        for (let k = 1; k < 4; k++) {
          flow.x += (labels[k] - labels[0]) * gradients[k].x;
          flow.y += (labels[k] - labels[0]) * gradients[k].y;
        }
        const w = weight * jacobian;
        for (let a = 0; a < 4; a++) if (ids[a] >= 0) {
          const ga = gradients[a], gaFlow = dot(ga, flow);
          residual[ids[a]] += w * gaFlow;
          if (linearize) for (let k = 0; k < 4; k++) if (ids[k] >= 0) {
            const gk = gradients[k], d = dir[k];
            // Material derivative for V=N_k D_k. No finite-difference
            // columns and no assumption that the general Jacobian is SPD.
            matrix.values[entries[4 * a + k]] += w * (gaFlow * dot(d, gk)
              - dot(ga, d) * dot(gk, flow) - dot(ga, gk) * dot(d, flow));
          }
        }
      }
    }
    if (!residual.every(Number.isFinite) || (linearize && !matrix.values.every(Number.isFinite)))
      throw new Error('Nonfinite transverse harmonic residual or Jacobian.');
    if (!linearize) return { residual };
    const rowScale = new Float64Array(n);
    for (let row = 0; row < n; row++) for (let p = matrix.rowPtr[row]; p < matrix.rowPtr[row + 1]; p++)
      rowScale[row] += Math.abs(matrix.values[p]);
    return { residual, matrix, rowScale };
  };
  return { nx, nt, n, evaluate, quadratureOrder,
    geometryModel: geometryCorrection ? 'Q1 node displacement plus a fixed curved correction' : 'bilinear Q1',
    scope: 'Residual and analytic material derivative only. Curved geometry has quadrature checks, not a whole-cell validity certificate.' };
}
