// SPDX-License-Identifier: GPL-2.0-or-later
// Trace a streamline with velocity potential as its parameter:
// d(position)/d(phi) = velocity / |velocity|^2. Increasing phi follows the
// flow even where physical x turns backward near a leading edge.
const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);

export function tracePotentialCurve({ seed, initialPotential = 0, endPotential, velocity,
  plane, tolerance = 1e-8, maxStep = .05, maxSpatialStep = .1, minStep = 1e-12, maxSteps = 20000,
  admissible = () => true, admissibleSegment = () => true }) {
  if (!finitePoint(seed) || ![initialPotential, endPotential, tolerance, maxStep, maxSpatialStep, minStep].every(Number.isFinite)
    || endPotential === initialPotential || Math.min(tolerance, maxStep, maxSpatialStep, minStep) <= 0 || minStep > maxStep
    || !Number.isInteger(maxSteps) || maxSteps < 1 || typeof velocity !== 'function' || !admissible(seed)
    || (plane && (!finitePoint(plane.normal) || !Number.isFinite(plane.offset) || Math.hypot(plane.normal.x, plane.normal.y) === 0)))
    throw new Error('Invalid potential-curve controls or initial point.');
  const sign = Math.sign(endPotential - initialPotential), planeLength = plane ? Math.hypot(plane.normal.x, plane.normal.y) : 1;
  const planeDistance = p => plane ? (plane.normal.x * p.x + plane.normal.y * p.y - plane.offset) / planeLength : null;
  const slope = p => {
    if (!finitePoint(p) || !admissible(p)) throw new Error('Potential curve left the admissible fluid region.');
    const q = velocity(p), speedSquared = q.u * q.u + q.v * q.v;
    if (!(speedSquared > 0) || !Number.isFinite(speedSquared)) throw new Error('Potential parameter is singular at stagnation.');
    return { x: q.u / speedSquared, y: q.v / speedSquared };
  };
  const rk4 = (point, step) => {
    const advance = (v, factor) => ({ x: point.x + step * factor * v.x, y: point.y + step * factor * v.y });
    const a = slope(point), b = slope(advance(a, .5)), c = slope(advance(b, .5)), d = slope(advance(c, 1));
    return { x: point.x + step * (a.x + 2 * b.x + 2 * c.x + d.x) / 6,
      y: point.y + step * (a.y + 2 * b.y + 2 * c.y + d.y) / 6 };
  };
  const trial = (point, step) => {
    const whole = rk4(point, step), half = rk4(rk4(point, step / 2), step / 2);
    const correction = { x: (half.x - whole.x) / 15, y: (half.y - whole.y) / 15 };
    const next = { x: half.x + correction.x, y: half.y + correction.y };
    if (!finitePoint(next) || !admissible(next) || !admissibleSegment(point, next)) throw new Error('Potential curve left the admissible fluid region.');
    return { next, error: Math.hypot(correction.x, correction.y) };
  };
  let point = { ...seed }, potential = initialPotential, stepSize = Math.min(maxStep, Math.abs(endPotential - initialPotential));
  const points = [{ ...point, potential }], startSide = Math.sign(planeDistance(point));
  if (plane && Math.abs(planeDistance(point)) <= tolerance) return { converged: true, reason: 'plane', points, accepted: 0, rejected: 0 };
  let accepted = 0, rejected = 0, lastRejection;
  for (let attempt = 0; attempt < maxSteps; attempt++) {
    const remaining = Math.abs(endPotential - potential);
    if (remaining <= 4 * Number.EPSILON * Math.max(1, Math.abs(endPotential)))
      return { converged: !plane, reason: plane ? 'potential limit before plane' : 'potential', points, accepted, rejected };
    // Restrict physical distance independently of the potential step, which
    // becomes a poor spatial scale near a stagnation point.
    let direction;
    try { direction = slope(point); } catch (error) { return { converged: false, reason: error.message, points, accepted, rejected }; }
    const step = sign * Math.min(stepSize, remaining, maxSpatialStep / Math.hypot(direction.x, direction.y));
    if (Math.abs(step) < minStep) return { converged: false, reason: lastRejection ?? 'minimum potential step', points, accepted, rejected };
    let result;
    try { result = trial(point, step); }
    catch (error) { lastRejection = error.message; stepSize = Math.abs(step) / 2; rejected++; continue; }
    if (result.error > tolerance || distance(point, result.next) > maxSpatialStep) {
      stepSize = Math.abs(step) * Math.max(.1, Math.min(.5, .8 * (tolerance / result.error) ** .2)); rejected++; continue;
    }
    let taken = step;
    if (plane && startSide * planeDistance(result.next) <= 0) {
      // Locate the event by reintegrating, not by clipping or linearly
      // projecting a curved streamline onto the requested plane.
      let lower = 0, upper = 1;
      for (let iteration = 0; iteration < 60 && Math.abs(planeDistance(result.next)) > tolerance; iteration++) {
        const fraction = (lower + upper) / 2; result = trial(point, step * fraction); taken = step * fraction;
        if (startSide * planeDistance(result.next) > 0) lower = fraction; else upper = fraction;
      }
      if (result.error > tolerance || Math.abs(planeDistance(result.next)) > tolerance)
        return { converged: false, reason: 'potential-curve plane event did not converge', points, accepted, rejected };
      points.push({ ...result.next, potential: potential + taken }); accepted++;
      return { converged: true, reason: 'plane', points, accepted, rejected };
    }
    point = result.next; potential += taken; points.push({ ...point, potential }); accepted++;
    stepSize = Math.min(maxStep, Math.abs(step) * Math.max(.5, Math.min(2, .9 * (tolerance / Math.max(result.error, 1e-30)) ** .2)));
  }
  return { converged: false, reason: 'potential-curve step limit', points, accepted, rejected };
}
