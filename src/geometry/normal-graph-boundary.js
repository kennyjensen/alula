// SPDX-License-Identifier: GPL-2.0-or-later
// A fixed C1 cubic-Hermite graph and its normal-foot equation. Applying this
// to a farfield is a mathematical extension, not Giles's literal x-only copy.
import { intervalPoint as pointInterval, intervalAdd, intervalSub, intervalMul, intervalDiv,
  nextUp, nextDown } from './bernstein-cell-certificate.js';

const finitePoint = p => Number.isFinite(p?.x) && Number.isFinite(p?.y);
const magnitude = interval => Math.max(Math.abs(interval[0]), Math.abs(interval[1]));
const hull = intervals => [Math.min(...intervals.map(v => v[0])), Math.max(...intervals.map(v => v[1]))];

export function createNormalGraphBoundary({ points, slopes } = {}) {
  if (!Array.isArray(points) || points.length < 2 || !Array.from(points).every(finitePoint)
    || points.some((p, i) => i && !(p.x > points[i - 1].x))
    || !Array.isArray(slopes) || slopes.length !== points.length || !Array.from(slopes).every(Number.isFinite))
    throw new Error('Normal graph requires finite, strictly increasing x knots and matching finite slopes.');
  points = Object.freeze(points.map(p => Object.freeze({ x: p.x, y: p.y })));
  slopes = Object.freeze([...slopes]);
  const xmin = points[0].x, xmax = points.at(-1).x, span = xmax - xmin;
  if (!(span > 0) || !Number.isFinite(span)) throw new Error('Normal graph extent is not representable.');
  const pieces = points.slice(1).map((p, i) => {
    const origin = points[i], h = p.x - origin.x, delta = p.y - origin.y, secant = delta / h;
    const left = slopes[i] - secant, right = slopes[i + 1] - secant;
    const cubic = h * (left + right), quadratic = -h * (2 * left + right), rightQuadratic = h * (left + 2 * right);
    if (![h, secant, cubic, quadratic, rightQuadratic, h * slopes[i], h * slopes[i + 1]].every(Number.isFinite))
      throw new Error('Normal graph Hermite coefficients are not representable.');
    // Bernstein coefficients relative to the first y avoid cancellation
    // under vertical translation. Intervals enclose exact-input Hermite data,
    // including the arithmetic needed to form h and the endpoint controls.
    const hi = intervalSub(pointInterval(p.x), pointInterval(origin.x));
    const dy = intervalSub(pointInterval(p.y), pointInterval(origin.y));
    const controls = [[0, 0], intervalDiv(intervalMul(hi, pointInterval(slopes[i])), [3, 3]),
      intervalSub(dy, intervalDiv(intervalMul(hi, pointInterval(slopes[i + 1])), [3, 3])), dy];
    const first = controls.slice(1).map((c, k) => intervalDiv(intervalMul([3, 3], intervalSub(c, controls[k])), hi));
    const second = first.slice(1).map((c, k) => intervalDiv(intervalMul([2, 2], intervalSub(c, first[k])), hi));
    const yRange = hull(controls.map(c => intervalAdd(pointInterval(origin.y), c))), slopeRange = hull(first), secondRange = hull(second);
    if (![...yRange, ...slopeRange, ...secondRange].every(Number.isFinite))
      throw new Error('Normal graph derivative bounds are not representable.');
    const slopeMinimum = slopeRange[0] > 0 ? slopeRange[0] : slopeRange[1] < 0 ? -slopeRange[1] : 0;
    const minimumSlopeSquared = slopeMinimum === 0 ? 0 : Math.max(0, nextDown(slopeMinimum * slopeMinimum));
    if (!Number.isFinite(minimumSlopeSquared)) throw new Error('Normal graph slope scale is not representable.');
    return { h, cubic, quadratic, rightQuadratic, yRange, curvature: magnitude(secondRange), minimumSlopeSquared };
  });
  const global = {
    yRange: [Math.min(...pieces.map(p => p.yRange[0])), Math.max(...pieces.map(p => p.yRange[1]))],
    curvature: Math.max(...pieces.map(p => p.curvature)),
    minimumSlopeSquared: Math.min(...pieces.map(p => p.minimumSlopeSquared)),
  };
  const indexAt = x => {
    if (!Number.isFinite(x) || x < xmin || x > xmax) throw new Error('Normal graph argument lies outside its knots.');
    let lo = 0, hi = points.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (x < points[mid].x) hi = mid; else lo = mid; }
    return lo;
  };
  const evaluate = x => {
    const i = indexAt(x), p = pieces[i], t = (x - points[i].x) / p.h;
    let y, slope;
    // Anchor each half at its nearby endpoint, preserving exact knot values
    // and slopes without subtracting nearly equal absolute coordinates.
    if (t <= .5) {
      y = points[i].y + t * (p.h * slopes[i] + t * (p.quadratic + t * p.cubic));
      slope = slopes[i] + t * (2 * p.quadratic + 3 * t * p.cubic) / p.h;
    } else {
      const v = (x - points[i + 1].x) / p.h;
      y = points[i + 1].y + v * (p.h * slopes[i + 1] + v * (p.rightQuadratic + v * p.cubic));
      slope = slopes[i + 1] + v * (2 * p.rightQuadratic + 3 * v * p.cubic) / p.h;
    }
    const second = (2 * p.quadratic / p.h + 6 * t * p.cubic / p.h) / p.h;
    if (![y, slope, second].every(Number.isFinite)) throw new Error('Normal graph evaluation is not representable.');
    return { point: { x, y }, derivative: { x: 1, y: slope }, secondDerivative: { x: 0, y: second } };
  };
  const derivativeLowerBound = (bounds, py) => {
    const distance = Math.max(Math.abs(nextDown(bounds.yRange[0] - py)), Math.abs(nextUp(bounds.yRange[1] - py)));
    const adverse = bounds.curvature === 0 ? 0 : nextUp(distance * bounds.curvature);
    return nextDown(nextDown(1 + bounds.minimumSlopeSquared) - adverse);
  };
  const projectNormal = p => {
    if (!finitePoint(p)) throw new Error('Normal projection requires a finite point.');
    // F'=1+f'^2+(f-py)f''. Positive lower bounds on every C1 piece imply
    // strict global monotonicity, even where the one-sided f'' values jump.
    // Try the cheap whole-domain bound before the sharper per-piece bounds.
    let lowerBound = derivativeLowerBound(global, p.y);
    if (!(lowerBound > 0)) {
      lowerBound = Infinity;
      for (const piece of pieces) lowerBound = Math.min(lowerBound, derivativeLowerBound(piece, p.y));
    }
    if (!(lowerBound > 0) || !Number.isFinite(lowerBound))
      throw new Error('Cannot establish a unique normal projection from the sufficient derivative bound.');
    const equation = x => {
      const q = evaluate(x), dx = x - p.x, dy = q.point.y - p.y, slope = q.derivative.y;
      const residual = dx + dy * slope, derivative = 1 + slope * slope + dy * q.secondDerivative.y;
      if (![residual, derivative].every(Number.isFinite) || !(derivative > 0))
        throw new Error('Normal projection equation is numerically unresolved.');
      return { ...q, residual, equationDerivative: derivative,
        geometricResidual: Math.abs(residual) / Math.hypot(1, slope) };
    };
    let lo = xmin, hi = xmax;
    const left = equation(lo), right = equation(hi);
    if (left.residual > 0 || right.residual < 0)
      throw new Error('Normal projection root is not bracketed inside the graph; extrapolation and endpoint clipping are not allowed.');
    if (left.residual === 0) return left.point;
    if (right.residual === 0) return right.point;
    let x = p.x > lo && p.x < hi ? p.x : lo + (hi - lo) / 2;
    for (let iteration = 0; iteration < 128; iteration++) {
      const q = equation(x);
      // Account for coordinate representation and evaluation roundoff. This
      // is a checked floating-point root, not an interval enclosure of it.
      const scale = Math.max(span, Math.abs(p.x), Math.abs(p.y), Math.abs(q.point.x), Math.abs(q.point.y));
      const tolerance = 64 * Number.EPSILON * scale;
      if (!(tolerance > 0) || !Number.isFinite(tolerance)) throw new Error('Normal projection residual scale is not representable.');
      if (q.geometricResidual <= tolerance) return q.point;
      if (q.residual < 0) lo = x; else hi = x;
      const width = hi - lo, candidate = x - q.residual / q.equationDerivative;
      // Restrict Newton to the central half of its current sign bracket;
      // otherwise bisect. This bounds interval reduction without changing F.
      const next = candidate > lo + .25 * width && candidate < hi - .25 * width ? candidate : lo + width / 2;
      if (!(next > lo && next < hi)) throw new Error('Normal projection bracket reached unresolved floating-point spacing.');
      x = next;
    }
    throw new Error('Normal projection exceeded its bounded iteration limit.');
  };
  const descriptor = Object.freeze({ kind: 'piecewise-cubic-hermite-graph', points, slopes });
  return Object.freeze({ points, slopes, descriptor, evaluate, projectNormal });
}
