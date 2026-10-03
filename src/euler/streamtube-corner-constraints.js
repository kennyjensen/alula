// SPDX-License-Identifier: GPL-2.0-or-later
// Signed corner areas of the actual primal quadrilaterals. These are step
// constraints, not residual rows or a replacement for nonlinear admissibility.
export function streamtubeCornerConstraints(euler, state, { displacementMap, minimumSine = 0, derivatives = true } = {}) {
  if (!Number.isFinite(minimumSine) || minimumSine < 0 || minimumSine >= 1) throw new Error('Invalid corner-sine margin.');
  const { nodes } = euler.decode(state), ne = euler.layout.n;
  const maps = derivatives ? euler.geometryDerivatives(state, { includeDisplacement: displacementMap !== undefined }) : null;
  const corners = [];
  nodes.forEach((region, g) => {
    for (let i = 0; i + 1 < region.length; i++) for (let j = 0; j + 1 < region[i].length; j++) {
      const ids = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]];
      const points = ids.map(([u, v]) => region[u][v]);
      for (let k = 0; k < 4; k++) {
        const a = points[k], b = points[(k + 1) % 4], c = points[(k + 2) % 4];
        const ab = { x: b.x - a.x, y: b.y - a.y }, bc = { x: c.x - b.x, y: c.y - b.y };
        const area = ab.x * bc.y - ab.y * bc.x, lengthAB = Math.hypot(ab.x, ab.y), lengthBC = Math.hypot(bc.x, bc.y);
        const value = area - minimumSine * lengthAB * lengthBC, gradient = new Map();
        // d[(b-a) x (c-b)]/d(a,b,c), before the normal-chart and BL chain rules.
        const pointGradients = [{ x: -bc.y, y: bc.x }, { x: c.y - a.y, y: a.x - c.x }, { x: -ab.y, y: ab.x }];
        if (minimumSine > 0) {
          const u = minimumSine * lengthBC / lengthAB, v = minimumSine * lengthAB / lengthBC;
          pointGradients[0].x += u * ab.x; pointGradients[0].y += u * ab.y;
          pointGradients[1].x -= u * ab.x - v * bc.x; pointGradients[1].y -= u * ab.y - v * bc.y;
          pointGradients[2].x -= v * bc.x; pointGradients[2].y -= v * bc.y;
        }
        const add = (col, v) => gradient.set(col, (gradient.get(col) ?? 0) + v);
        for (let m = 0; m < 3; m++) {
          const [u, v] = ids[(k + m) % 4], p = pointGradients[m];
          for (const [col, d] of maps?.[g][u][v] ?? []) {
            const derivative = p.x * d.x + p.y * d.y;
            if (col < ne) add(col, derivative);
            else {
              if (!displacementMap?.[col - ne]) throw new Error('Missing corner displacement chain rule.');
              for (const [target, factor] of displacementMap[col - ne]) add(target, factor * derivative);
            }
          }
        }
        for (const [col, v] of gradient) if (v === 0) gradient.delete(col);
        corners.push({ g, i, j, k, value, area, sine: area / (lengthAB * lengthBC), gradient });
      }
    }
  });
  return corners;
}
