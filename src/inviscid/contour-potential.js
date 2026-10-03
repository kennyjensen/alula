// SPDX-License-Identifier: GPL-2.0-or-later
import { potentialDifference } from './streamfunction.js';
import { velocityAt } from './linear-vortex.js';

const same = (a, b) => a.x === b.x && a.y === b.y;
const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);
const MAX_DEPTH = 12, MAX_SEGMENTS = 65536;

// A potential branch is defined by its path, not by an independently fitted
// surface constant. Follow the supplied contour outward from stagnation on
// its two sides. Each accepted chord uses the analytic linear-vortex-panel
// primitive. Subdivision resolves visibility, not quadrature accuracy.
//
// The curve must remain unchanged during use. Field data are copied because
// changing strengths after construction would invalidate the prefix cache.
export function createContourPotential({ curve, stagnation, stagnationPotential, field, subdivisions = 2 } = {}) {
  if (!curve || !Number.isFinite(curve.length) || !(curve.length > 0) || typeof curve.evaluate !== 'function'
    || !curve.knots || curve.knots.length < 2 || !Number.isFinite(stagnationPotential)
    || !(Number.isFinite(stagnation) && stagnation > 0 && stagnation < curve.length)
    || !Number.isInteger(subdivisions) || subdivisions < 1 || subdivisions > 32)
    throw new Error('Invalid contour-potential curve, stagnation, or subdivision count.');
  const knots = Array.from(curve.knots);
  if (knots[0] !== 0 || knots.at(-1) !== curve.length || !knots.every((s, i) => Number.isFinite(s) && (!i || s > knots[i - 1])))
    throw new Error('Contour-potential knots must strictly span the curve.');
  if (!field || !Array.isArray(field.panels) || !field.gamma || ![field.u, field.v].every(Number.isFinite)
    || !Array.from(field.gamma).every(Number.isFinite) || !field.panels.every(p => finitePoint(p.a) && finitePoint(p.b)
      && [p.length, p.tx, p.ty].every(Number.isFinite) && p.length > 0 && Number.isInteger(p.node)
      && p.node >= 0 && p.node + 1 < field.gamma.length)) throw new Error('Invalid contour-potential panel field.');
  if (field.basePanels !== undefined && (!Array.isArray(field.basePanels) || !field.basePanels.every(p => finitePoint(p.a) && finitePoint(p.b)
    && [p.length, p.tx, p.ty, p.sourceStrength, p.vortexStrength].every(Number.isFinite) && p.length > 0
    && finitePoint(p.cutDirection) && (p.cutOrigin === undefined || finitePoint(p.cutOrigin)))))
    throw new Error('Invalid contour-potential finite-base field.');
  field = { u: field.u, v: field.v, gamma: Float64Array.from(field.gamma),
    panels: field.panels.map(p => ({ ...p, a: { ...p.a }, b: { ...p.b } })),
    ...(field.basePanels !== undefined ? { basePanels: field.basePanels.map(p => ({ ...p, a: { ...p.a }, b: { ...p.b },
      cutDirection: { ...p.cutDirection }, ...(p.cutOrigin !== undefined ? { cutOrigin: { ...p.cutOrigin } } : {}) })) } : {}) };
  const sheets = field.basePanels?.length ? [...field.panels, ...field.basePanels] : field.panels;

  // Closed panel chains define solid interiors. Open sheets still receive
  // potentialDifference's crossing check. These floating-point geometric
  // checks detect incompatible paths; they are not a whole-curve certificate.
  const groups = new Map();
  for (const p of sheets) {
    if (!groups.has(p.element)) groups.set(p.element, []);
    groups.get(p.element).push(p);
  }
  // A finite-base field contains open wetted sheets followed by the actual
  // retained base segments. Their union, not the wetted chain alone, bounds
  // the solid. Reject missing/misordered base geometry rather than silently
  // treating such a body as an open sheet with no interior.
  for (const element of new Set((field.basePanels ?? []).map(p => p.element))) {
    const ps = groups.get(element);
    if (!ps || ps.length < 3 || !ps.every((p, i) => same(p.b, ps[(i + 1) % ps.length].a)))
      throw new Error('Finite-base contour-potential sheets must form the complete ordered closed body boundary.');
  }
  const contours = [...groups.values()].filter(ps => ps.length >= 3
    && ps.every((p, i) => same(p.b, ps[(i + 1) % ps.length].a))).map(ps => ps.map(p => p.a));
  const inside = point => contours.some(points => {
    let result = false;
    for (let i = 0, j = points.length - 1; i < points.length; j = i++) {
      const a = points[j], b = points[i];
      const cross = (b.x - a.x) * (point.y - a.y) - (b.y - a.y) * (point.x - a.x);
      if (same(point, a) || same(point, b) || (cross === 0 && point.x >= Math.min(a.x, b.x)
        && point.x <= Math.max(a.x, b.x) && point.y >= Math.min(a.y, b.y) && point.y <= Math.max(a.y, b.y))) return false;
      if ((a.y > point.y) !== (b.y > point.y)
        && point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x) result = !result;
    }
    return result;
  });
  const diagnostics = { method: 'analytic panel potential along knot-split contour paths', subdivisions,
    sampleCount: 0, analyticSegments: 0, adaptiveSubdivisions: 0, maximumSubdivisionDepth: 0,
    geometryCheck: 'sampled fluid membership and chord sheet crossings; not a whole-curve certificate' };
  const incompatible = (message, interval) => {
    const error = new Error(`Contour-potential path is incompatible with the panel sheets: ${message}`);
    error.code = 'CONTOUR_POTENTIAL_PATH_INCOMPATIBLE';
    error.interval = interval;
    error.diagnostics = { ...diagnostics };
    return error;
  };
  const pointAt = s => {
    const p = curve.evaluate(s)?.point;
    if (!finitePoint(p)) throw new Error('Nonfinite contour-potential curve point.');
    if (inside(p)) throw incompatible('a sampled contour point lies inside a panel body.', [s, s]);
    for (const panel of sheets) {
      if (same(p, panel.a) || same(p, panel.b)) continue;
      const dx = panel.b.x - panel.a.x, dy = panel.b.y - panel.a.y;
      const px = p.x - panel.a.x, py = p.y - panel.a.y, dot = px * dx + py * dy;
      if (px * dy - py * dx === 0 && dot > 0 && dot < dx * dx + dy * dy)
        throw incompatible('a sampled contour point lies on a sheet interior with an unspecified one-sided limit.', [s, s]);
    }
    return p;
  };
  // Budget applies to a single path evaluation as well as to initialization;
  // normal use can evaluate arbitrarily many separate phase queries.
  let segmentBudget = 0;
  const increment = (a, pa, b, pb, depth = 0) => {
    if (a === b) return 0;
    if (++segmentBudget > MAX_SEGMENTS) throw incompatible('path subdivision budget exhausted.', [a, b]);
    diagnostics.maximumSubdivisionDepth = Math.max(diagnostics.maximumSubdivisionDepth, depth);
    const middle = a + .5 * (b - a), pm = pointAt(middle);
    let split = inside({ x: .5 * (pa.x + pb.x), y: .5 * (pa.y + pb.y) });
    if (!split) {
      try {
        const value = potentialDifference(pa, pb, field);
        diagnostics.analyticSegments++;
        return value;
      } catch (error) {
        if (!/crosses a (?:vortex|finite-base) sheet/.test(error.message)) throw error;
        split = true;
      }
    }
    if (split && (depth >= MAX_DEPTH || middle === a || middle === b))
      throw incompatible('a chord still crosses a sheet at the subdivision limit.', [a, b]);
    diagnostics.adaptiveSubdivisions++;
    return increment(a, pa, middle, pm, depth + 1) + increment(middle, pm, b, pb, depth + 1);
  };

  const breaks = knots.includes(stagnation) ? knots : [...knots, stagnation].sort((a, b) => a - b);
  if ((breaks.length - 1) * subdivisions + 1 > MAX_SEGMENTS)
    throw new Error('Contour-potential prefix cache exceeds its point budget.');
  const parameters = [breaks[0]];
  for (let i = 1; i < breaks.length; i++) for (let j = 1; j <= subdivisions; j++)
    parameters.push(j === subdivisions ? breaks[i] : breaks[i - 1] + j / subdivisions * (breaks[i] - breaks[i - 1]));
  if (!parameters.every((s, i) => !i || s > parameters[i - 1])) throw new Error('Contour-potential subdivisions are below parameter resolution.');
  const points = parameters.map(pointAt), relative = new Float64Array(parameters.length), index = parameters.indexOf(stagnation);
  diagnostics.sampleCount = parameters.length;
  const accumulate = direction => {
    let sum = 0, correction = 0;
    for (let i = index + direction; i >= 0 && i < parameters.length; i += direction) {
      const previous = i - direction;
      const term = increment(parameters[previous], points[previous], parameters[i], points[i]) - correction;
      const next = sum + term;
      correction = (next - sum) - term;
      relative[i] = sum = next;
    }
  };
  accumulate(-1); accumulate(1);
  const phase = s => {
    if (!Number.isFinite(s) || s < 0 || s > curve.length) throw new Error('Contour-potential parameter is outside the curve.');
    if (s === stagnation) return stagnationPotential;
    let lo = 0, hi = parameters.length;
    while (hi - lo > 1) { const mid = (hi + lo) >> 1; if (parameters[mid] > s) hi = mid; else lo = mid; }
    if (parameters[lo] === s) return stagnationPotential + relative[lo];
    segmentBudget = 0;
    return stagnationPotential + relative[lo] + increment(parameters[lo], points[lo], s, pointAt(s));
  };
  const branch = direction => {
    let minimumOutwardIncrement = Infinity, negativeIncrements = 0;
    for (let i = index + direction; i >= 0 && i < relative.length; i += direction) {
      const delta = relative[i] - relative[i - direction];
      minimumOutwardIncrement = Math.min(minimumOutwardIncrement, delta);
      if (delta < 0) negativeIncrements++;
    }
    return { minimumOutwardIncrement, negativeIncrements, sampledMonotone: negativeIncrements === 0 };
  };
  diagnostics.branches = { upper: branch(-1), lower: branch(1) };
  const derivative = s => {
    const value = curve.evaluate(s);
    if (!finitePoint(value.derivative)) return null;
    // Potential has a finite branch limit at a panel vertex, although the
    // unit-sheet velocity queried by velocityAt is undefined there. Keep
    // that distinction; adjacent probes still assess branch monotonicity.
    if (sheets.some(p => same(value.point, p.a) || same(value.point, p.b))) return null;
    const q = velocityAt(value.point, field), d = q.u * value.derivative.x + q.v * value.derivative.y;
    return Number.isFinite(d) ? d : null;
  };
  const h = Math.min(stagnation - parameters[index - 1], parameters[index + 1] - stagnation) / 256;
  const upperDerivative = derivative(stagnation - h), lowerDerivative = derivative(stagnation + h);
  diagnostics.nearStagnation = { parameterStep: h, derivative: derivative(stagnation),
    upperOutwardDerivative: upperDerivative === null ? null : -upperDerivative,
    lowerOutwardDerivative: lowerDerivative,
    upperIncrement: phase(stagnation - h) - stagnationPotential,
    lowerIncrement: phase(stagnation + h) - stagnationPotential,
    sampledMonotone: upperDerivative !== null && lowerDerivative !== null && upperDerivative <= 0 && lowerDerivative >= 0 };
  diagnostics.sampledMonotone = diagnostics.branches.upper.sampledMonotone && diagnostics.branches.lower.sampledMonotone
    && diagnostics.nearStagnation.sampledMonotone;
  diagnostics.contourPotentialIncrement = relative.at(-1) - relative[0];
  return { phase, diagnostics };
}
