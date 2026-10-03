// SPDX-License-Identifier: GPL-2.0-or-later
// Source-comparison control only. Giles SMOVE pp.280–282 uses a single
// reference-bank coordinate distribution over every transverse row. Reuse
// the production stiffness/line solver, changing only that prescribed RHS.
import { assembleTangentialCoordinate } from '../../src/geometry/streamtube-tangential-redistribution.js';
import { quadLaplaceMatrix } from '../../src/numerics/quad-laplace.js';
import { solveAlternatingScalarLines } from '../../src/numerics/alternating-scalar-lines.js';

export function redistributeWithSharedCoordinate(nodes, referenceNodes, options = {}) {
  const coordinate = assembleTangentialCoordinate(nodes, options);
  const { nx, nt, firstStreamline, lastStreamline, lineDimensions, matrix } = coordinate;
  if (!Array.isArray(referenceNodes) || referenceNodes.length !== nx + 1
    || !referenceNodes.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))
    throw new Error('Shared SMOVE coordinate requires a complete finite reference bank.');
  const increments = Float64Array.from({ length: nx }, (_, i) => Math.sqrt(Math.sqrt(
    Math.hypot(referenceNodes[i + 1].x - referenceNodes[i].x, referenceNodes[i + 1].y - referenceNodes[i].y))));
  if (!increments.every(v => Number.isFinite(v) && v > 0)) throw new Error('Degenerate shared reference interval.');
  const id = (i, j) => i > 0 && i < nx && j >= firstStreamline && j <= lastStreamline
    ? (i - 1) * (lineDimensions.nt - 1) + j - firstStreamline : null;
  const rhs = new Float64Array(matrix.n);
  for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
    const ids = [id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)];
    const k = quadLaplaceMatrix([nodes[i][j], nodes[i + 1][j], nodes[i + 1][j + 1], nodes[i][j + 1]],
      { quadratureOrder: 2, quadratureDomain: options.quadratureDomain ?? 'convex' });
    const local = [0, increments[i], increments[i], 0];
    for (let a = 0; a < 4; a++) if (ids[a] !== null) for (let b = 0; b < 4; b++)
      rhs[ids[a]] -= k[4 * a + b] * (local[b] - local[a]);
  }
  const solution = solveAlternatingScalarLines(matrix, rhs, { ...lineDimensions, pairs: options.pairs ?? 5 });
  const correction = (i, j) => id(i, j) === null ? 0 : solution.x[id(i, j)];
  const moved = nodes.map(row => row.map(p => ({ ...p })));
  let maxDisplacement = 0;
  for (let i = 1; i < nx; i++) for (let j = firstStreamline; j <= lastStreamline; j++) {
    const span = increments[i] + increments[i - 1] + correction(i + 1, j) - correction(i - 1, j);
    if (!(span > 0) || !Number.isFinite(span)) throw new Error('Nonpositive shared-coordinate span.');
    const fraction = -correction(i, j) / span;
    const dx = fraction * (nodes[i + 1][j].x - nodes[i - 1][j].x);
    const dy = fraction * (nodes[i + 1][j].y - nodes[i - 1][j].y);
    moved[i][j].x += dx; moved[i][j].y += dy;
    if (![moved[i][j].x, moved[i][j].y].every(Number.isFinite)) throw new Error('Nonfinite shared-coordinate movement.');
    maxDisplacement = Math.max(maxDisplacement, Math.hypot(dx, dy));
  }
  return { nodes: moved, increments, solution, maxDisplacement };
}
