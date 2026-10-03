// SPDX-License-Identifier: GPL-2.0-or-later
// Streamline distance s is the integration coordinate:
// dr/ds = V/|V|, dphi/ds = |V|. Potential therefore remains a measured
// coordinate, rather than prescribing physical station spacing near stagnation.
const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);

export function traceArcCurve({ seed, initialArc = 0, endArc, initialPotential = 0, velocity,
  plane, tolerance = 1e-8, potentialTolerance = tolerance, maxStep = .1, minStep = 1e-12,
  maxSteps = 20000, admissible = () => true, admissibleSegment = () => true }) {
  if (!finitePoint(seed) || ![initialArc, endArc, initialPotential, tolerance, potentialTolerance, maxStep, minStep].every(Number.isFinite)
    || initialArc === endArc || Math.min(tolerance, potentialTolerance, maxStep, minStep) <= 0 || minStep > maxStep
    || !Number.isInteger(maxSteps) || maxSteps < 1 || typeof velocity !== 'function' || !admissible(seed)
    || (plane && (!finitePoint(plane.normal) || !Number.isFinite(plane.offset) || Math.hypot(plane.normal.x, plane.normal.y) === 0)))
    throw new Error('Invalid arc-curve controls or initial point.');
  const sign = Math.sign(endArc - initialArc), planeLength = plane ? Math.hypot(plane.normal.x, plane.normal.y) : 1;
  const planeDistance = p => plane ? (plane.normal.x * p.x + plane.normal.y * p.y - plane.offset) / planeLength : null;
  const slope = p => {
    if (!finitePoint(p) || !admissible(p)) throw new Error('Arc curve left the admissible fluid region.');
    const q = velocity(p), speed = Math.hypot(q.u, q.v);
    if (!(speed > 0) || !Number.isFinite(speed)) throw new Error('Streamline direction is singular at stagnation.');
    return { x: q.u / speed, y: q.v / speed, potential: speed };
  };
  const rk4 = (p, h) => {
    const advance = (v, fraction) => ({ x: p.x + h * fraction * v.x, y: p.y + h * fraction * v.y });
    const a = slope(p), b = slope(advance(a, .5)), c = slope(advance(b, .5)), d = slope(advance(c, 1));
    return Object.fromEntries(['x', 'y', 'potential'].map(k => [k, p[k] + h * (a[k] + 2 * b[k] + 2 * c[k] + d[k]) / 6]));
  };
  const trial = (point, step) => {
    const whole = rk4(point, step), half = rk4(rk4(point, step / 2), step / 2);
    const correction = Object.fromEntries(['x', 'y', 'potential'].map(k => [k, (half[k] - whole[k]) / 15]));
    const next = Object.fromEntries(['x', 'y', 'potential'].map(k => [k, half[k] + correction[k]]));
    if (!finitePoint(next) || !Number.isFinite(next.potential) || !admissible(next) || !admissibleSegment(point, next))
      throw new Error('Arc curve left the admissible fluid region.');
    return { next, error: Math.max(Math.hypot(correction.x, correction.y) / tolerance, Math.abs(correction.potential) / potentialTolerance) };
  };
  let point = { x: seed.x, y: seed.y, potential: initialPotential }, arc = initialArc;
  let stepSize = Math.min(maxStep, Math.abs(endArc - initialArc)), accepted = 0, rejected = 0, lastRejection;
  const points = [{ ...point, arc }], startSide = Math.sign(planeDistance(point));
  if (plane && Math.abs(planeDistance(point)) <= tolerance) return { converged: true, reason: 'plane', points, accepted, rejected };
  for (let attempt = 0; attempt < maxSteps; attempt++) {
    const remaining = Math.abs(endArc - arc);
    if (remaining <= 4 * Number.EPSILON * Math.max(Math.abs(initialArc), Math.abs(endArc), Number.MIN_VALUE))
      return { converged: !plane, reason: plane ? 'arc limit before plane' : 'arc', points, accepted, rejected };
    const step = sign * Math.min(stepSize, remaining);
    if (Math.abs(step) < minStep) return { converged: false, reason: lastRejection ?? 'minimum arc step', points, accepted, rejected };
    let result;
    try { result = trial(point, step); }
    catch (error) { lastRejection = error.message; stepSize = Math.abs(step) / 2; rejected++; continue; }
    if (result.error > 1) { stepSize = Math.abs(step) * Math.max(.1, Math.min(.5, .8 * result.error ** -.2)); rejected++; continue; }
    let taken = step;
    if (plane && startSide * planeDistance(result.next) <= 0) {
      // Reintegrate every event trial. Projecting a chord onto the plane
      // would lose both physical arc and potential accuracy on curved flow.
      let lower = 0, upper = 1;
      for (let iteration = 0; iteration < 60 && Math.abs(planeDistance(result.next)) > tolerance; iteration++) {
        const fraction = (lower + upper) / 2; taken = step * fraction; result = trial(point, taken);
        if (startSide * planeDistance(result.next) > 0) lower = fraction; else upper = fraction;
      }
      if (result.error > 1 || Math.abs(planeDistance(result.next)) > tolerance)
        return { converged: false, reason: 'arc-curve plane event did not converge', points, accepted, rejected };
      points.push({ ...result.next, arc: arc + taken }); accepted++;
      return { converged: true, reason: 'plane', points, accepted, rejected };
    }
    point = result.next; arc += taken; points.push({ ...point, arc }); accepted++;
    stepSize = Math.min(maxStep, Math.abs(step) * Math.max(.5, Math.min(2, .9 * Math.max(result.error, 1e-30) ** -.2)));
  }
  return { converged: false, reason: 'arc-curve step limit', points, accepted, rejected };
}
