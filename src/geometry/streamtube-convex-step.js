// SPDX-License-Identifier: GPL-2.0-or-later
// A fraction-to-boundary safeguard for a proposed coordinate update. Every
// quad corner orientation is quadratic along an affine coordinate path.
// Its first positive zero bounds the connected convex domain containing
// the current grid; checking only the proposed endpoint misses interior
// folds (including a double zero). This does not certify a nonlinear wall
// reconstruction, gas admissibility, or nonlocal mesh intersections.
export const STREAMTUBE_CONVEX_MINIMUM_CORNER_SINE = 1e-12; // Same domain as streamtubeMeshSnapshot.
const MINIMUM_CORNER_SINE = STREAMTUBE_CONVEX_MINIMUM_CORNER_SINE;
const cross = (ux, uy, vx, vy) => ux * vy - uy * vx;

export class StreamtubeGridStepError extends Error {
  constructor(message, code, details = {}) {
    super(message); this.name = 'StreamtubeGridStepError'; this.code = code; this.details = details;
    const { group, i, tube, corner } = details;
    this.diagnostics = { code, ...details, ...(Number.isInteger(group) && Number.isInteger(i) && Number.isInteger(tube)
      ? { cell: { group, i, tube, ...(Number.isInteger(corner) ? { corner } : {}) } } : {}) };
  }
}

function validateShape(from, to) {
  if (!Array.isArray(from) || !from.length || !Array.isArray(to) || to.length !== from.length)
    throw new StreamtubeGridStepError('Grid update requires matching nonempty passages.', 'streamtube-grid-step-shape');
  for (let g = 0; g < from.length; g++) {
    const a = from[g], b = to[g], width = a?.[0]?.length;
    if (!Array.isArray(a) || a.length < 2 || !Array.isArray(b) || b.length !== a.length || !(width >= 2))
      throw new StreamtubeGridStepError('Grid update requires matching structured passages.', 'streamtube-grid-step-shape', { group: g });
    for (let i = 0; i < a.length; i++) {
      if (!Array.isArray(a[i]) || !Array.isArray(b[i]) || a[i].length !== width || b[i].length !== width)
        throw new StreamtubeGridStepError('Grid update requires matching streamtube rows.', 'streamtube-grid-step-shape', { group: g, i });
      for (let j = 0; j < width; j++) for (const [nodes, label] of [[a, 'starting'], [b, 'proposed']]) {
        if (!Number.isFinite(nodes[i][j]?.x) || !Number.isFinite(nodes[i][j]?.y))
          throw new StreamtubeGridStepError(`Nonfinite ${label} grid update coordinate.`, 'streamtube-grid-step-coordinate', { group: g, i, tube: j, state: label });
      }
    }
  }
}

// Exact endpoint coordinates matter when a full update is already safe.
const mix = (a, b, t) => t === 0 ? { ...a } : t === 1 ? { ...b }
  : Object.is(a.x, b.x) && Object.is(a.y, b.y) ? { ...a }
  : { ...a, x: a.x + t * (b.x - a.x), y: a.y + t * (b.y - a.y) };

export function interpolateStreamtubeGridNodes(fromNodes, toNodes, step) {
  validateShape(fromNodes, toNodes);
  if (!Number.isFinite(step) || step < 0 || step > 1)
    throw new StreamtubeGridStepError('Grid interpolation step must lie between zero and one.', 'streamtube-grid-step-control');
  return fromNodes.map((group, g) => group.map((row, i) => row.map((p, j) => mix(p, toNodes[g][i][j], step))));
}

function cornerSine(a, b, c) {
  const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
  const u = Math.hypot(ux, uy), v = Math.hypot(vx, vy);
  return cross(ux / u, uy / u, vx / v, vy / v);
}

/** Check actual supplied coordinates; no interpolation or geometry adoption. */
export function streamtubeGridConvexity(nodes) {
  validateShape(nodes, nodes);
  let minCornerSine = Infinity, cellsChecked = 0;
  const invalidCells = [];
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length - 1; i++)
    for (let j = 0; j < nodes[g][i].length - 1; j++) {
      const p = [nodes[g][i][j], nodes[g][i + 1][j], nodes[g][i + 1][j + 1], nodes[g][i][j + 1]];
      let invalid = null;
      cellsChecked++;
      for (let k = 0; k < 4; k++) {
        const sine = cornerSine(p[k], p[(k + 1) % 4], p[(k + 2) % 4]);
        minCornerSine = Math.min(minCornerSine, Number.isFinite(sine) ? sine : -Infinity);
        if (!(sine > MINIMUM_CORNER_SINE)) invalid ??= { group: g, i, tube: j, corner: k, sine };
      }
      if (invalid) invalidCells.push(invalid);
    }
  return { valid: invalidCells.length === 0, minCornerSine, invalidCells, cellsChecked };
}

export function assertConvexStreamtubeGrid(nodes) {
  const quality = streamtubeGridConvexity(nodes);
  if (!quality.valid) {
    const { group, i, tube, corner } = quality.invalidCells[0];
    throw new StreamtubeGridStepError(`Nonconvex grid cell i=${i}, group=${group}, tube=${tube}, corner=${corner}.`,
      'streamtube-grid-nonconvex', quality.invalidCells[0]);
  }
  return quality;
}

// Normalize the two edge families separately. This avoids determinant
// overflow/underflow and preserves roots under changes of physical units,
// including highly unequal streamwise and transverse edge lengths.
function cornerPolynomial(p, q, k) {
  const a = p[k], b = p[(k + 1) % 4], c = p[(k + 2) % 4];
  const A = q[k], B = q[(k + 1) % 4], C = q[(k + 2) % 4];
  const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
  const Ux = B.x - A.x, Uy = B.y - A.y, Vx = C.x - B.x, Vy = C.y - B.y;
  const su = Math.max(Math.abs(ux), Math.abs(uy), Math.abs(Ux), Math.abs(Uy));
  const sv = Math.max(Math.abs(vx), Math.abs(vy), Math.abs(Vx), Math.abs(Vy));
  const u0 = ux / su, u1 = uy / su, v0 = vx / sv, v1 = vy / sv;
  const du0 = Ux / su - u0, du1 = Uy / su - u1, dv0 = Vx / sv - v0, dv1 = Vy / sv - v1;
  return { coefficients: [cross(du0, du1, dv0, dv1), cross(du0, du1, v0, v1) + cross(u0, u1, dv0, dv1), cross(u0, u1, v0, v1)],
    edges: { u0, u1, v0, v1, du0, du1, dv0, dv1 } };
}

function firstPositiveRoot(coefficients, edges, maximumStep) {
  let [a, b, c] = coefficients;
  const scale = Math.max(Math.abs(a), Math.abs(b), Math.abs(c));
  a /= scale; b /= scale; c /= scale;
  if (a === 0) return b < 0 ? -c / b : Infinity;
  let discriminant = b * b - 4 * a * c;
  const vertex = -b / (2 * a);
  if (a > 0 && vertex > 0 && vertex <= maximumStep
    && Math.abs(discriminant) <= Math.sqrt(Number.EPSILON) * (b * b + Math.abs(4 * a * c))) {
    // Nearly coincident roots lose their separating discriminant in the
    // expanded polynomial. Recenter using the actual affine edge vectors
    // instead. Merely clamping a small negative discriminant would invent
    // a collapse for a rotation arbitrarily close to 180 degrees, whose
    // interpolated quads remain orthogonal and strictly convex.
    const { u0, u1, v0, v1, du0, du1, dv0, dv1 } = edges;
    const U0 = u0 + vertex * du0, U1 = u1 + vertex * du1;
    const V0 = v0 + vertex * dv0, V1 = v1 + vertex * dv1;
    a = coefficients[0]; b = cross(du0, du1, V0, V1) + cross(U0, U1, dv0, dv1); c = cross(U0, U1, V0, V1);
    const localScale = Math.max(Math.abs(a), Math.abs(b), Math.abs(c));
    a /= localScale; b /= localScale; c /= localScale;
    discriminant = b * b - 4 * a * c;
    if (discriminant < 0 && discriminant >= -8 * Number.EPSILON * (b * b + Math.abs(4 * a * c))) discriminant = 0;
    if (discriminant < 0) return Infinity;
    const root = Math.sqrt(discriminant), q = -.5 * (b + (b < 0 ? -root : root));
    const x = vertex + q / a, y = vertex + (q === 0 ? 0 : c / q);
    return Math.min(x > 0 ? x : Infinity, y > 0 ? y : Infinity);
  }
  if (discriminant < 0) return Infinity;
  const root = Math.sqrt(discriminant);
  // q/a and c/q avoid subtracting nearly equal numbers for either root.
  const q = -.5 * (b + (b < 0 ? -root : root));
  const x = q / a, y = q === 0 ? Infinity : c / q;
  return Math.min(x > 0 ? x : Infinity, y > 0 ? y : Infinity);
}

/** Limit an affine coordinate step without changing any coordinates.
 *
 * The returned step is absolute on from + step * (to - from). A bound
 * found on the requested path is reduced by safetyFraction; unrestricted
 * steps are returned exactly. The actual interpolated final grid also
 * passes the existing strict corner-sine threshold. Callers that rebuild
 * nodes nonlinearly must check that rebuilt candidate separately.
 */
export function limitStreamtubeGridStep(fromNodes, toNodes, { maximumStep = 1, safetyFraction = .95 } = {}) {
  validateShape(fromNodes, toNodes);
  if (!Number.isFinite(maximumStep) || maximumStep < 0 || maximumStep > 1
    || !Number.isFinite(safetyFraction) || !(safetyFraction > 0 && safetyFraction < 1))
    throw new StreamtubeGridStepError('Invalid convex grid step controls.', 'streamtube-grid-step-control');
  let step = maximumStep, boundaryStep = Infinity, limiter = null, cellsChecked = 0;
  const cells = [];
  for (let g = 0; g < fromNodes.length; g++) for (let i = 0; i < fromNodes[g].length - 1; i++)
    for (let j = 0; j < fromNodes[g][i].length - 1; j++) {
      const cell = grid => [grid[g][i][j], grid[g][i + 1][j], grid[g][i + 1][j + 1], grid[g][i][j + 1]];
      const p = cell(fromNodes), q = cell(toNodes), location = { group: g, i, tube: j };
      cells.push({ p, q, location }); cellsChecked++;
      for (let k = 0; k < 4; k++) {
        const sine = cornerSine(p[k], p[(k + 1) % 4], p[(k + 2) % 4]);
        if (!(sine > MINIMUM_CORNER_SINE))
          throw new StreamtubeGridStepError('Convex grid update requires an admissible starting cell.', 'streamtube-grid-invalid-start', { ...location, corner: k, sine });
        const { coefficients, edges } = cornerPolynomial(p, q, k);
        if (!coefficients.every(Number.isFinite) || !(coefficients[2] > 0))
          throw new StreamtubeGridStepError('Unresolved grid update orientation polynomial.', 'streamtube-grid-step-resolution', { ...location, corner: k });
        const root = firstPositiveRoot(coefficients, edges, maximumStep);
        if (root <= maximumStep && root < boundaryStep) {
          boundaryStep = root; limiter = { ...location, corner: k, boundaryStep: root, reason: 'corner-orientation' };
        }
      }
    }
  if (limiter) step = safetyFraction * boundaryStep;
  const qualityAt = t => {
    let minCornerSine = Infinity, invalid = null;
    for (const { p, q, location } of cells) {
      const points = p.map((a, k) => mix(a, q[k], t));
      for (let k = 0; k < 4; k++) {
        const sine = cornerSine(points[k], points[(k + 1) % 4], points[(k + 2) % 4]);
        minCornerSine = Math.min(minCornerSine, Number.isFinite(sine) ? sine : -Infinity);
        if (!(sine > MINIMUM_CORNER_SINE)) invalid ??= { ...location, corner: k, sine };
      }
    }
    return { valid: invalid === null, minCornerSine, invalid };
  };
  let quality = qualityAt(step), validationBacktracks = 0;
  // Root roundoff and near-collinear endpoints still face the actual
  // geometry predicate. This safeguard never changes the domain tolerance.
  while (!quality.valid && validationBacktracks < 60 && step > 0) {
    limiter = { ...quality.invalid, boundaryStep: step, reason: 'corner-sine' };
    step *= .5; validationBacktracks++; quality = qualityAt(step);
  }
  if (!quality.valid || maximumStep > 0 && step === 0)
    throw new StreamtubeGridStepError('No resolved positive convex grid update step.', 'streamtube-grid-step-resolution', limiter ?? {});
  return { step, limited: step < maximumStep, limiter, quality, cellsChecked, validationBacktracks };
}
