// SPDX-License-Identifier: GPL-2.0-or-later
// A potential-ordered guide may turn in physical x. Select a unique resolved
// x-plane intersection without sorting, clipping, or discarding guide points.
// Uniqueness concerns the resolved polyline brackets, not an unproved claim
// about oscillations hidden between trace samples.
function failure(code, message, diagnostics) {
  return Object.assign(new Error(message), { code, diagnostics });
}

export function findUniqueGuideXBracket(path, x) {
  if (!Array.isArray(path) || path.length < 2 || !Number.isFinite(x)
    || !path.every((p, i) => p && [p.x, p.y, p.potential].every(Number.isFinite)
      && (!i || p.potential > path[i - 1].potential)))
    throw new Error('A guide needs finite points in strictly increasing potential order.');
  const vertices = [], intervals = [];
  for (let i = 0; i < path.length; i++) {
    if (path[i].x === x) vertices.push(i);
    if (!i) continue;
    const a = path[i - 1].x, b = path[i].x;
    if (a === b && a === x)
      throw failure('GUIDE_X_AMBIGUOUS', 'The requested x-plane contains an unresolved vertical guide interval.',
        { x, interval: [i - 1, i] });
    // Open intervals exclude endpoints, so adjacent brackets sharing an
    // exact vertex are counted once by the vertex list above.
    if (Math.min(a, b) < x && x < Math.max(a, b)) intervals.push([i - 1, i]);
  }
  const count = vertices.length + intervals.length;
  if (count !== 1)
    throw failure(count ? 'GUIDE_X_AMBIGUOUS' : 'GUIDE_X_UNRESOLVED', count
      ? 'The requested x-plane has multiple resolved guide intersections.'
      : 'The requested x-plane has no resolved guide intersection.', { x, vertices, intervals });
  if (vertices.length) return { vertex: vertices[0], lower: vertices[0], upper: vertices[0], orientation: 0 };
  const [lower, upper] = intervals[0];
  return { vertex: null, lower, upper, orientation: Math.sign(path[upper].x - path[lower].x) };
}

export function sampleUniqueGuideX(path, x, { sample, velocity, tolerance, maxIterations = 40 }) {
  if (typeof sample !== 'function' || typeof velocity !== 'function' || !(tolerance > 0)
    || !Number.isFinite(tolerance) || !Number.isInteger(maxIterations) || maxIterations < 1)
    throw new Error('Invalid unique guide-intersection controls.');
  const bracket = findUniqueGuideXBracket(path, x);
  const directed = (lower, upper) => {
    const orientation = Math.sign(path[upper].x - path[lower].x);
    const a = velocity(path[lower]), b = velocity(path[upper]);
    if (!(orientation * a.u > 0 && orientation * b.u > 0))
      throw failure('GUIDE_X_UNRESOLVED', 'The requested guide intersection touches an unresolved x turning interval.',
        { x, interval: [lower, upper], orientation, endpointVelocityX: [a.u, b.u] });
  };
  if (bracket.vertex !== null) {
    // A sampled extremum need not be the continuous trace extremum. For
    // example, an interval leaving a vertex with u>0 and ending at smaller
    // x contains another crossing of that vertex's x-plane. Reject such
    // unresolved turns rather than counting the common vertex as unique.
    if (bracket.vertex > 0) directed(bracket.vertex - 1, bracket.vertex);
    if (bracket.vertex < path.length - 1) directed(bracket.vertex, bracket.vertex + 1);
    return { ...path[bracket.vertex] };
  }
  directed(bracket.lower, bracket.upper);
  const a = path[bracket.lower], b = path[bracket.upper];
  let lower = a.potential, upper = b.potential;
  let potential = lower + (upper - lower) * (x - a.x) / (b.x - a.x);
  for (let k = 0; k < maxIterations; k++) {
    const point = sample(path, potential);
    if (![point?.x, point?.y].every(Number.isFinite)) throw new Error('Nonfinite guide-intersection sample.');
    const error = point.x - x;
    if (Math.abs(error) <= tolerance) return point;
    if (bracket.orientation * error < 0) lower = potential;
    else upper = potential;
    // The existing potential parameter satisfies dr/dphi = q / |q|².
    // Bracketing remains valid for either x orientation; a singular or
    // out-of-bracket Newton proposal falls back to the bracket midpoint.
    const q = velocity(point), next = potential - error * (q.u * q.u + q.v * q.v) / q.u;
    potential = next > lower && next < upper ? next : .5 * (lower + upper);
  }
  throw new Error('Unique streamline/x-plane intersection did not converge.');
}
