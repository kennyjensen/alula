// SPDX-License-Identifier: GPL-2.0-or-later
// Mathematical reconstruction of the SMOVE coordinate equation and movement
// formula (Giles pp.280–282). Body/cut banks have zero coordinate correction;
// farfield banks adjoining the excluded logical strip use the natural weak
// condition, not an artificial fixed coordinate. Boundary roles belong to
// the caller; no mesh clipping, line search or density remap occurs.
import { quadLaplaceMatrix } from '../numerics/quad-laplace.js';
import { sparseMatrix, sparseAdd } from '../numerics/sparse.js';
import { solveAlternatingScalarLines } from '../numerics/alternating-scalar-lines.js';

export function assembleTangentialCoordinate(nodes, { referenceBank = 0, fixedBanks = [true, true], quadratureDomain = 'convex' } = {}) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  if (!Array.isArray(nodes) || nx < 2 || nt < 2 || !nodes.every(row => Array.isArray(row) && row.length === nt + 1
    && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))) || ![0, nt].includes(referenceBank)
    || !Array.isArray(fixedBanks) || fixedBanks.length !== 2 || !fixedBanks.every(v => typeof v === 'boolean'))
    throw new Error('Invalid tangential redistribution grid or reference bank.');
  const increments = new Float64Array(nx), baseline = new Float64Array(nx + 1);
  for (let i = 0; i < nx; i++) {
    const a = nodes[i][referenceBank], b = nodes[i + 1][referenceBank];
    increments[i] = Math.sqrt(Math.sqrt(Math.hypot(b.x - a.x, b.y - a.y)));
    if (!(increments[i] > 0)) throw new Error('Degenerate reference-bank segment.');
    baseline[i + 1] = baseline[i] + increments[i];
  }
  const firstStreamline = fixedBanks[0] ? 1 : 0, lastStreamline = fixedBanks[1] ? nt - 1 : nt;
  const count = lastStreamline - firstStreamline + 1, n = (nx - 1) * count;
  const id = (i, j) => i > 0 && i < nx && j >= firstStreamline && j <= lastStreamline ? (i - 1) * count + j - firstStreamline : null;
  const rows = Array.from({ length: n }, () => new Map()), rhs = new Float64Array(n);
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ids = [id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)];
    const k = quadLaplaceMatrix([nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1]], { quadratureOrder: 2, quadratureDomain });
    const local = [0, increments[i], increments[i], 0];
    for (let a = 0; a < 4; a++) if (ids[a] !== null) for (let b = 0; b < 4; b++) {
      // Subtract the row value using partition of unity, avoiding cancellation
      // of large cumulative coordinates. Dirichlet correction values are zero.
      rhs[ids[a]] -= k[4 * a + b] * (local[b] - local[a]);
      if (ids[b] !== null) rows[ids[a]].set(ids[b], (rows[ids[a]].get(ids[b]) ?? 0) + k[4 * a + b]);
    }
  }
  const matrix = sparseMatrix(rows.map(row => [...row.keys()]));
  rows.forEach((row, r) => { for (const [c, v] of row) sparseAdd(matrix, r, c, v); });
  return { matrix, rhs, nx, nt, lineDimensions: { nx, nt: count + 1 }, firstStreamline, lastStreamline,
    increments, baseline, referenceBank, fixedBanks: fixedBanks.slice() };
}

export function redistributeStreamtubeTangentially(nodes, { referenceBank = 0, fixedBanks = [true, true], pairs = 5, quadratureDomain = 'convex', correctionScale = 1 } = {}) {
  if (!Number.isFinite(correctionScale) || !(correctionScale > 0 && correctionScale <= 1))
    throw new Error('Invalid tangential coordinate correction scale.');
  const coordinate = assembleTangentialCoordinate(nodes, { referenceBank, fixedBanks, quadratureDomain });
  const { nx, increments, firstStreamline, lastStreamline, lineDimensions } = coordinate;
  const solution = solveAlternatingScalarLines(coordinate.matrix, coordinate.rhs, { ...lineDimensions, pairs });
  // Apply the common scale to the coordinate field, including the span in
  // the movement denominator. The five-pair equation itself is unchanged.
  const correction = (i, j) => i > 0 && i < nx && j >= firstStreamline && j <= lastStreamline
    ? correctionScale * solution.x[(i - 1) * (lineDimensions.nt - 1) + j - firstStreamline] : 0;
  const moved = nodes.map(row => row.map(p => ({ ...p })));
  let maxDisplacement = 0, minCoordinateSpan = Infinity;
  for (let i = 1; i < nx; i++) for (let j = firstStreamline; j <= lastStreamline; j++) {
    const span = increments[i] + increments[i - 1] + correction(i + 1, j) - correction(i - 1, j);
    minCoordinateSpan = Math.min(minCoordinateSpan, span);
    if (!(span > 0) || !Number.isFinite(span)) throw new Error(`Nonpositive tangential coordinate span at station ${i}, streamline ${j}.`);
    const fraction = -correction(i, j) / span, a = nodes[i - 1][j], b = nodes[i + 1][j];
    const dx = fraction * (b.x - a.x), dy = fraction * (b.y - a.y);
    moved[i][j].x += dx; moved[i][j].y += dy;
    if (!Number.isFinite(moved[i][j].x) || !Number.isFinite(moved[i][j].y)) throw new Error('Nonfinite tangential grid movement.');
    maxDisplacement = Math.max(maxDisplacement, Math.hypot(dx, dy));
  }
  return { nodes: moved, solution, referenceBank, fixedBanks: fixedBanks.slice(), maxDisplacement, minCoordinateSpan };
}
