// SPDX-License-Identifier: GPL-2.0-or-later
import { quadrilateralCornerShape } from './quadrilateral-shape.js';

// Geometry-quality initialization only. Constant tube masses and the flow
// equations are not part of this minimization. Begin from a positive mesh;
// each accepted vertex update decreases shape energy and retains every
// incident corner's orientation. Fixed boundaries are never moved.
export function smoothQuadrilaterals({ vertices, cells, fixed = new Set() }, {
  aspects, maxSweeps = 50, tolerance = 1e-8,
} = {}) {
  if (!Array.isArray(vertices) || !vertices.length || !vertices.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))
    || !Array.isArray(cells) || !cells.length || !cells.every(c => Array.isArray(c) && c.length === 4 && new Set(c).size === 4
      && c.every(i => Number.isInteger(i) && i >= 0 && i < vertices.length))
    || !(fixed instanceof Set) || [...fixed].some(i => !Number.isInteger(i) || i < 0 || i >= vertices.length)
    || !Number.isInteger(maxSweeps) || maxSweeps < 0 || !Number.isFinite(tolerance) || tolerance <= 0
    || (aspects !== undefined && (!Array.isArray(aspects) || aspects.length !== cells.length
      || !aspects.every(row => Array.isArray(row) && row.length === 4 && row.every(r => Number.isFinite(r) && r > 0)))))
    throw new Error('Invalid quadrilateral smoothing controls or mesh.');
  const original = vertices.map(p => ({ ...p })), points = vertices.map(p => ({ ...p }));
  const corners = [], incident = points.map(() => []), lengths = points.map(() => Infinity);
  const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
  cells.forEach((cell, i) => cell.forEach((a, k) => {
    const b = cell[(k + 1) % 4], c = cell[(k + 2) % 4], ab = distance(points[a], points[b]), bc = distance(points[b], points[c]);
    const corner = { ids: [a, b, c], aspect: aspects?.[i][k] ?? ab / bc };
    quadrilateralCornerShape(corner.ids.map(id => points[id]), corner.aspect, 0);
    const index = corners.length; corners.push(corner);
    corner.ids.forEach((id, slot) => incident[id].push({ index, slot }));
    lengths[a] = Math.min(lengths[a], ab); lengths[b] = Math.min(lengths[b], ab);
  }));
  const local = (id, trial = points[id]) => {
    let value = 0; const gradient = [0, 0], hessian = [0, 0, 0, 0], constraints = [];
    for (const { index, slot } of incident[id]) {
      const c = corners[index], r = quadrilateralCornerShape(c.ids.map(j => j === id ? trial : points[j]), c.aspect, slot);
      value += r.value; gradient.forEach((_, k) => { gradient[k] += r.gradient[k]; });
      hessian.forEach((_, k) => { hessian[k] += r.hessian[k]; }); constraints.push(r);
    }
    return { value, gradient, hessian, constraints };
  };
  const diagnostics = () => {
    let energy = 0, minCornerSine = Infinity, maxShape = 0, gradient = 0;
    for (const c of corners) {
      const p = c.ids.map(id => points[id]), r = quadrilateralCornerShape(p, c.aspect, 0);
      energy += r.value; maxShape = Math.max(maxShape, r.value);
      minCornerSine = Math.min(minCornerSine, r.determinant / (distance(p[0], p[1]) * distance(p[1], p[2])));
    }
    for (let id = 0; id < points.length; id++) if (!fixed.has(id) && incident[id].length) {
      const r = local(id); gradient = Math.max(gradient, Math.hypot(...r.gradient) * lengths[id] / r.value);
    }
    return { energy, maxShape, minCornerSine, gradient };
  };
  let q = diagnostics(), reason = 'sweep limit'; const history = [{ sweep: 0, ...q }];
  for (let sweep = 0; sweep < maxSweeps && q.gradient > tolerance; sweep++) {
    let updates = 0;
    for (let id = 0; id < points.length; id++) {
      if (fixed.has(id) || !incident[id].length) continue;
      const r = local(id), [gx, gy] = r.gradient, [hxx, hxy, , hyy] = r.hessian;
      if (Math.hypot(gx, gy) * lengths[id] / r.value <= tolerance) continue;
      // Scale the 2x2 system before inversion to avoid unnecessary products
      // of large physical-coordinate curvatures on very short edges.
      const scale = Math.max(Math.abs(hxx), Math.abs(hxy), Math.abs(hyy));
      const a = hxx / scale, b = hxy / scale, c = hyy / scale, det = a * c - b * b;
      if (!(det > 0) || !(scale > 0)) continue;
      const dx = (-c * gx + b * gy) / (scale * det), dy = (b * gx - a * gy) / (scale * det);
      const slope = gx * dx + gy * dy; if (!(slope < 0) || ![dx, dy].every(Number.isFinite)) continue;
      let step = 1;
      for (const corner of r.constraints) {
        const change = corner.determinantGradient[0] * dx + corner.determinantGradient[1] * dy;
        if (change < 0) step = Math.min(step, .9 * corner.determinant / -change);
      }
      for (let backtrack = 0; backtrack < 30; backtrack++, step *= .5) {
        const trial = { x: points[id].x + step * dx, y: points[id].y + step * dy };
        let candidate;
        try { candidate = local(id, trial); } catch { continue; }
        if (candidate.value <= r.value + 1e-4 * step * slope) { points[id] = trial; updates++; break; }
      }
    }
    q = diagnostics(); history.push({ sweep: sweep + 1, ...q, updates });
    if (!updates) { reason = 'Vertex updates made no progress.'; break; }
  }
  const converged = q.gradient <= tolerance;
  return { converged, reason: converged ? 'shape gradient' : reason, vertices: points, history, quality: { valid: q.minCornerSine > 0, ...q },
    aspects: cells.map((_, i) => corners.slice(4 * i, 4 * i + 4).map(c => c.aspect)),
    maxDisplacement: Math.max(...points.map((p, i) => distance(p, original[i]))),
    status: 'geometry quality only; flow equations have not been solved' };
}
