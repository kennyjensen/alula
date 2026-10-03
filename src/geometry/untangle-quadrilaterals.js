// SPDX-License-Identifier: GPL-2.0-or-later
// Geometry-only feasibility repair. For one moving vertex, every incident
// signed corner determinant is affine in its two coordinates. Intersect
// those half-planes in a bounded local frame, then project toward the original
// position. If infeasible, minimize the maximum normalized constraint defect.
// This is not a fluid solve, and no repaired grid is certified as such.
const determinant = (a, b, c) => (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);

function clip(polygon, { x, y, bound }, slack) {
  const output = [];
  for (let i = 0; i < polygon.length; i++) {
    const a = polygon[i], b = polygon[(i + 1) % polygon.length];
    const fa = x * a.x + y * a.y - bound + slack, fb = x * b.x + y * b.y - bound + slack;
    if (fa >= 0) output.push(a);
    if ((fa >= 0) !== (fb >= 0)) {
      const t = fa / (fa - fb); output.push({ x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) });
    }
  }
  return output;
}

function closestToOrigin(polygon) {
  let best = null, squared = Infinity;
  for (let i = 0; i < polygon.length; i++) {
    const a = polygon[i], b = polygon[(i + 1) % polygon.length], dx = b.x - a.x, dy = b.y - a.y;
    const lengthSquared = dx * dx + dy * dy;
    const t = lengthSquared > 0 ? Math.max(0, Math.min(1, -(a.x * dx + a.y * dy) / lengthSquared)) : 0;
    const p = { x: a.x + t * dx, y: a.y + t * dy }, value = p.x * p.x + p.y * p.y;
    if (value < squared) { best = p; squared = value; }
  }
  return best;
}

export function untangleQuadrilaterals({ vertices, cells, fixed = new Set() }, {
  minimumCorner = 1e-3, displacementLimit = 2, maxSweeps = 100, minimumNormalAreas,
} = {}) {
  if (!Array.isArray(vertices) || !vertices.length || !vertices.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))
    || !Array.isArray(cells) || !cells.length || !cells.every(c => Array.isArray(c) && c.length === 4 && new Set(c).size === 4
      && c.every(i => Number.isInteger(i) && i >= 0 && i < vertices.length))
    || !(fixed instanceof Set) || [...fixed].some(i => !Number.isInteger(i) || i < 0 || i >= vertices.length)
    || ![minimumCorner, displacementLimit].every(Number.isFinite) || minimumCorner <= 0 || minimumCorner >= 1 || displacementLimit <= 0
    || !Number.isInteger(maxSweeps) || maxSweeps < 0
    || (minimumNormalAreas !== undefined && (!Array.isArray(minimumNormalAreas) || minimumNormalAreas.length !== cells.length
      || minimumNormalAreas.some(a => !Number.isFinite(a) || a < 0)))) throw new Error('Invalid quadrilateral repair controls or mesh.');
  const original = vertices.map(p => ({ ...p })), points = vertices.map(p => ({ ...p }));
  // This tolerance applies only to the auxiliary repair margin. Corner
  // positivity is checked independently. Avoid activating different vertices
  // after rigid translations merely from roundoff in short edge differences.
  const marginTolerance = minimumCorner * 1e-7;
  const corners = [], incident = vertices.map(() => []), neighbors = vertices.map(() => new Set());
  cells.forEach((cell, cellIndex) => cell.forEach((a, k) => {
    const b = cell[(k + 1) % 4], c = cell[(k + 2) % 4];
    const scale = distance(original[a], original[b]) * distance(original[b], original[c]);
    if (!(scale > 0) || !Number.isFinite(scale)) throw new Error('Unresolved reference edge in quadrilateral repair.');
    const value = determinant(original[a], original[b], original[c]) / scale;
    const corner = { ids: [a, b, c], scale, target: value > 1e-12 ? Math.min(minimumCorner, .5 * value) : minimumCorner, cell: cellIndex };
    const index = corners.length; corners.push(corner);
    for (const id of [a, b, c]) incident[id].push(index);
    neighbors[a].add(b); neighbors[b].add(a);
  }));
  // A section's normal width is A/|d|, with A its signed quadrilateral
  // area and d the vector between its upstream/downstream edge centers.
  // With one moving vertex, A and d are affine: A - amin*|d| >= 0 is
  // therefore a convex feasible set. Tangent cuts of this concave function
  // give outer approximations; every accepted local update is checked
  // against the nonlinear constraint, not just its tangent.
  const areaIncident = vertices.map(() => []), areaConstraints = [];
  const areaValue = (constraint, moving = -1, replacement) => {
    const p = constraint.ids.map(i => i === moving && replacement ? replacement : points[i]);
    const dx = .5 * ((p[1].x - p[0].x) + (p[2].x - p[3].x));
    const dy = .5 * ((p[1].y - p[0].y) + (p[2].y - p[3].y));
    const length = Math.hypot(dx, dy);
    const area = .5 * (determinant(p[0], p[1], p[2]) + determinant(p[0], p[2], p[3]));
    const value = (area - constraint.minimum * length) / constraint.scale;
    if (moving < 0) return { value, ratio: length > 0 ? area / (constraint.minimum * length) : 0 };
    const slot = constraint.ids.indexOf(moving), previous = p[(slot + 3) % 4], next = p[(slot + 1) % 4];
    const sign = slot === 1 || slot === 2 ? .5 : -.5;
    // At d=0, zero is a valid subgradient of its Euclidean norm.
    return { value, gradient: {
      x: (.5 * (next.y - previous.y) - constraint.minimum * sign * (length > 0 ? dx / length : 0)) / constraint.scale,
      y: (.5 * (previous.x - next.x) - constraint.minimum * sign * (length > 0 ? dy / length : 0)) / constraint.scale,
    } };
  };
  minimumNormalAreas?.forEach((minimum, cell) => {
    if (minimum === 0) return;
    const ids = cells[cell], p = ids.map(i => original[i]);
    const length = .5 * Math.hypot(p[1].x - p[0].x + p[2].x - p[3].x, p[1].y - p[0].y + p[2].y - p[3].y);
    if (!(length > 0)) throw new Error('Unresolved reference section direction in quadrilateral repair.');
    const index = areaConstraints.length;
    areaConstraints.push({ ids, minimum, scale: minimum * length, cell });
    ids.forEach(id => areaIncident[id].push(index));
  });
  const frames = original.map((p, i) => {
    if (fixed.has(i) || !neighbors[i].size) return null;
    const adjacent = [...neighbors[i]], q = original[adjacent[0]], length = distance(p, q);
    const radius = displacementLimit * Math.min(...adjacent.map(j => distance(p, original[j])));
    if (!(radius > 0) || !Number.isFinite(radius)) throw new Error('Unresolved quadrilateral repair movement scale.');
    return { x: (q.x - p.x) / length, y: (q.y - p.y) / length, radius };
  });
  const quality = () => {
    let maximumDefect = 0, minCornerSine = Infinity, minNormalAreaRatio = Infinity;
    const invalid = new Set(), active = new Set();
    corners.forEach((corner, k) => {
      const [a, b, c] = corner.ids.map(id => points[id]), det = determinant(a, b, c);
      const value = det / corner.scale, defect = corner.target - value;
      const sine = det / (distance(a, b) * distance(b, c));
      minCornerSine = Math.min(minCornerSine, sine); maximumDefect = Math.max(maximumDefect, defect);
      if (!(sine > 1e-12)) invalid.add(corner.cell);
      if (defect > marginTolerance || !(sine > 1e-12)) corner.ids.forEach(id => { if (!fixed.has(id)) active.add(id); });
    });
    const invalidNormalAreaCells = [];
    areaConstraints.forEach(constraint => {
      const r = areaValue(constraint);
      minNormalAreaRatio = Math.min(minNormalAreaRatio, r.ratio); maximumDefect = Math.max(maximumDefect, -r.value);
      if (r.value < -marginTolerance || !(r.ratio > 0)) {
        invalidNormalAreaCells.push(constraint.cell);
        constraint.ids.forEach(id => { if (!fixed.has(id)) active.add(id); });
      }
    });
    return { valid: invalid.size === 0 && invalidNormalAreaCells.length === 0, maximumDefect, minCornerSine,
      invalidCells: [...invalid], invalidNormalAreaCells, minNormalAreaRatio, active: [...active] };
  };
  let q = quality(), reason = 'sweep limit';
  const history = [{ sweep: 0, invalidCells: q.invalidCells.length, maximumDefect: q.maximumDefect, minCornerSine: q.minCornerSine }];
  for (let sweep = 0; sweep < maxSweeps && (!q.valid || q.maximumDefect > marginTolerance); sweep++) {
    let maxMove = 0;
    if (!q.active.length) { reason = 'Fixed boundary corners prevent repair.'; break; }
    for (const id of q.active) {
      const frame = frames[id], origin = original[id], point = points[id];
      const coordinates = p => ({ x: ((p.x - origin.x) * frame.x + (p.y - origin.y) * frame.y) / frame.radius,
        y: (-(p.x - origin.x) * frame.y + (p.y - origin.y) * frame.x) / frame.radius });
      const current = coordinates(point);
      const constraints = incident[id].map(k => {
        const corner = corners[k], [a, b, c] = corner.ids.map(j => points[j]), slot = corner.ids.indexOf(id);
        const gradients = [{ x: b.y - c.y, y: c.x - b.x }, { x: c.y - a.y, y: a.x - c.x }, { x: a.y - b.y, y: b.x - a.x }];
        const gradient = gradients[slot];
        const x = frame.radius * (gradient.x * frame.x + gradient.y * frame.y) / corner.scale;
        const y = frame.radius * (-gradient.x * frame.y + gradient.y * frame.x) / corner.scale;
        return { x, y, bound: corner.target - determinant(a, b, c) / corner.scale + x * current.x + y * current.y };
      });
      const feasible = slack => {
        let polygon = [{ x: -1, y: -1 }, { x: 1, y: -1 }, { x: 1, y: 1 }, { x: -1, y: 1 }];
        for (const constraint of constraints) { polygon = clip(polygon, constraint, slack); if (!polygon.length) break; }
        return polygon;
      };
      const locate = () => {
        let polygon = feasible(0), usedSlack = 0;
        if (!polygon.length) {
          let lower = 0, upper = Math.max(...constraints.map(c => c.bound - c.x * current.x - c.y * current.y), 0);
          // Add only a roundoff allowance for the current point's feasibility.
          upper += 64 * Number.EPSILON * Math.max(1, upper);
          polygon = feasible(upper);
          if (!polygon.length) return null;
          for (let iteration = 0; iteration < 48; iteration++) {
            const middle = .5 * (lower + upper), trial = feasible(middle);
            if (trial.length) { upper = middle; polygon = trial; } else lower = middle;
          }
          usedSlack = upper;
        }
        // Test origin membership against the half-planes themselves. A
        // minimax polygon can collapse to a segment/point; its zero edge
        // cross-products must not incorrectly classify the origin as inside.
        const local = constraints.every(c => c.bound <= usedSlack) ? { x: 0, y: 0 } : closestToOrigin(polygon);
        return local && { local, usedSlack };
      };
      const toPhysical = local => ({ x: origin.x + frame.radius * (local.x * frame.x - local.y * frame.y),
        y: origin.y + frame.radius * (local.x * frame.y + local.y * frame.x) });
      const areaCut = (index, local) => {
        const r = areaValue(areaConstraints[index], id, toPhysical(local));
        const x = frame.radius * (r.gradient.x * frame.x + r.gradient.y * frame.y);
        const y = frame.radius * (-r.gradient.x * frame.y + r.gradient.y * frame.x);
        return { value: r.value, cut: { x, y, bound: -r.value + x * local.x + y * local.y } };
      };
      areaIncident[id].forEach(index => constraints.push(areaCut(index, current).cut));
      let projection;
      for (let iteration = 0; iteration < 32; iteration++) {
        projection = locate(); if (!projection) break;
        const violated = areaIncident[id].map(index => areaCut(index, projection.local))
          .filter(r => r.value < -projection.usedSlack - marginTolerance);
        if (!violated.length) break;
        constraints.push(...violated.map(r => r.cut)); projection = null;
      }
      if (!projection) continue;
      const { local } = projection;
      const next = { x: origin.x + frame.radius * (local.x * frame.x - local.y * frame.y),
        y: origin.y + frame.radius * (local.x * frame.y + local.y * frame.x) };
      maxMove = Math.max(maxMove, distance(next, point) / frame.radius); points[id] = next;
    }
    q = quality(); history.push({ sweep: sweep + 1, invalidCells: q.invalidCells.length, maximumDefect: q.maximumDefect, minCornerSine: q.minCornerSine, maxMove });
    if (maxMove === 0) { reason = 'Bounded vertex updates made no progress.'; break; }
  }
  const converged = q.valid && q.maximumDefect <= marginTolerance;
  const maxDisplacement = Math.max(...points.map((p, i) => distance(p, original[i])));
  return { converged, reason: converged ? areaConstraints.length ? 'positive corners and minimum normal areas' : 'positive corner constraints' : reason, vertices: points, history,
    quality: { valid: q.valid, invalidCells: q.invalidCells, minCornerSine: q.minCornerSine,
      ...(areaConstraints.length ? { invalidNormalAreaCells: q.invalidNormalAreaCells, minNormalAreaRatio: q.minNormalAreaRatio } : {}) }, maxDisplacement,
    status: 'geometry feasibility only; flow equations have not been solved' };
}
