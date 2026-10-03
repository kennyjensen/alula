// SPDX-License-Identifier: GPL-2.0-or-later
import { velocityAt } from './linear-vortex.js';
import { createContourPotential } from './contour-potential.js';

const same = (a, b) => a.x === b.x && a.y === b.y;
const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);

// The supplied contour is piecewise cubic between its knots. Its restricted
// cubic Bezier control polygon bounds the physical arc between two parameters;
// endpoint distance alone would underestimate uncertainty on a curved wall.
function physicalBracketWidth(curve, left, right) {
  const breaks = [left, ...curve.knots.filter(s => s > left && s < right), right];
  let width = 0;
  for (let i = 1; i < breaks.length; i++) {
    const a = curve.evaluate(breaks[i - 1]), b = curve.evaluate(breaks[i]), h = (breaks[i] - breaks[i - 1]) / 3;
    width += h * (Math.hypot(a.derivative.x, a.derivative.y) + Math.hypot(b.derivative.x, b.derivative.y))
      + Math.hypot(b.point.x - a.point.x - h * (a.derivative.x + b.derivative.x),
        b.point.y - a.point.y - h * (a.derivative.y + b.derivative.y));
  }
  return width;
}

// Negative/positive endpoint derivatives locate an interior minimum of a
// continuous surface-potential branch. At a panel vertex that minimum need
// not have a defined classical velocity. Report a narrow minimum bracket
// separately; never replace an undefined velocity with zero or move a query
// off the prescribed wall. Ordinary off-vertex roots retain the velocity
// residual criterion used by the initializer.
export function refineSurfaceStagnation({ curve, field, left, right, chord, relativeTolerance = 2e-9, maxIterations = 56 }) {
  if (!curve || !Array.isArray(curve.knots) || typeof curve.evaluate !== 'function'
    || ![left, right, chord, relativeTolerance].every(Number.isFinite) || !(left >= 0 && left < right && right <= curve.length)
    || !(chord > 0 && relativeTolerance > 0) || !Number.isInteger(maxIterations) || maxIterations < 1
    || !field || !Array.isArray(field.panels) || !Number.isFinite(Math.hypot(field.u, field.v)) || !(Math.hypot(field.u, field.v) > 0))
    throw new Error('Invalid analytic surface-stagnation bracket.');
  const sheets = field.basePanels?.length ? [...field.panels, ...field.basePanels] : field.panels;
  const derivative = s => {
    const value = curve.evaluate(s);
    if (!finitePoint(value.point) || !finitePoint(value.derivative)) throw new Error('Nonfinite surface-stagnation curve.');
    if (sheets.some(p => same(value.point, p.a) || same(value.point, p.b))) return null;
    const q = velocityAt(value.point, field), d = q.u * value.derivative.x + q.v * value.derivative.y;
    if (!Number.isFinite(d)) throw new Error('Nonfinite surface tangential velocity.');
    return d;
  };
  const original = { left, right }, physicalTolerance = relativeTolerance * chord;
  let lo = left, hi = right, dlo = derivative(lo), dhi = derivative(hi), vertex = null, comparison = null, iterations = 0;
  if (!(dlo !== null && dhi !== null && dlo < 0 && dhi > 0))
    throw new Error('Analytic surface stagnation is not bracketed by finite signed panel velocities.');
  const bracket = () => ({ left: lo, right: hi, leftDerivative: dlo, rightDerivative: dhi,
    physicalWidth: physicalBracketWidth(curve, lo, hi), physicalTolerance });
  const compareVertex = s => {
    // Use the original finite bracket, where potential changes are resolved,
    // rather than subtracting nearly equal potentials after shrinking it.
    // The contour-potential helper checks its knot-split paths for panel-sheet
    // incompatibilities. The zero here is only a potential gauge at the vertex.
    const length = original.right - original.left;
    const local = { length, knots: [0, ...curve.knots.filter(k => k > original.left && k < original.right)
      .map(k => k - original.left), length], evaluate: t => curve.evaluate(original.left + t) };
    const branch = createContourPotential({ curve: local, field, stagnation: s - original.left, stagnationPotential: 0 });
    const leftPotential = branch.phase(0), rightPotential = branch.phase(length);
    return { left: original.left, vertex: s, right: original.right, leftPotential, vertexPotential: 0, rightPotential,
      vertexBelowEndpoints: leftPotential > 0 && rightPotential > 0,
      pathCheck: branch.diagnostics.geometryCheck, analyticSegments: branch.diagnostics.analyticSegments };
  };
  const resolvedVertex = () => {
    if (vertex === null || !(lo < vertex && vertex < hi) || !(dlo < 0 && dhi > 0) || !comparison?.vertexBelowEndpoints) return null;
    const interval = bracket();
    if (!(interval.physicalWidth <= physicalTolerance)) return null;
    return { parameter: vertex, kind: 'potential-minimum-bracket', derivative: null,
      derivativeStatus: 'undefined at exact panel vertex', bracket: interval, comparison, iterations };
  };
  for (; iterations < maxIterations; iterations++) {
    const resolved = resolvedVertex();
    if (resolved) return resolved;
    const mid = lo + .5 * (hi - lo);
    if (mid === lo || mid === hi) break;
    const dmid = derivative(mid);
    if (dmid !== null) {
      if (dmid < 0) { lo = mid; dlo = dmid; } else { hi = mid; dhi = dmid; }
      continue;
    }
    if (vertex !== mid || comparison === null) comparison = compareVertex(mid);
    vertex = mid;
    const resolvedAtMid = resolvedVertex();
    if (resolvedAtMid) return resolvedAtMid;
    const a = lo + .5 * (mid - lo), b = mid + .5 * (hi - mid);
    if (!(lo < a && a < mid && mid < b && b < hi)) break;
    const da = derivative(a), db = derivative(b);
    if (da === null || db === null) throw new Error('Surface-stagnation quarter probes coincide with unresolved panel vertices.');
    if (da >= 0) { hi = a; dhi = da; }
    else if (db <= 0) { lo = b; dlo = db; }
    else { lo = a; dlo = da; hi = b; dhi = db; }
  }
  const resolved = resolvedVertex();
  if (resolved) return resolved;
  if (vertex !== null && lo < vertex && vertex < hi)
    throw new Error('Analytic surface potential minimum bracket did not resolve its physical width and endpoint comparison.');
  const parameter = lo + .5 * (hi - lo), d = derivative(parameter), tangent = curve.evaluate(parameter).derivative;
  const velocityTolerance = relativeTolerance * Math.hypot(field.u, field.v) * Math.hypot(tangent.x, tangent.y);
  if (d === null || Math.abs(d) > velocityTolerance)
    throw new Error('Analytic surface stagnation did not resolve a continuous tangential-velocity zero.');
  return { parameter, kind: 'velocity-zero', derivative: d, derivativeStatus: 'finite', velocityTolerance,
    bracket: bracket(), comparison, iterations };
}
