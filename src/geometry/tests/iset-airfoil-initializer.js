// SPDX-License-Identifier: GPL-2.0-or-later
// Independent implementation of the rounded-LE, closed-TE isolated-airfoil
// NORMIT/SPLNIT/OUTLIN/NORLIN/RESPLI mathematics in Giles (1985), pp175–181.
// This builds geometry from supplied distributions, without a panel/flow solve.
import { createContourCurve } from '../contour-curve.js';
import { gilesStreamwiseCoordinates } from '../elliptic-streamtube-grid.js';
import { signedArea } from '../airfoil.js';

function normalizeStations(values, name, duplicate = false) {
  if (!Array.isArray(values) || values.length < 2 || !values.every(Number.isFinite)
    || values.some((v, i) => i && (duplicate ? v < values[i - 1] : v <= values[i - 1])))
    throw new Error(`ISET ${name} stations must be finite and ordered.`);
  const span = values.at(-1) - values[0];
  if (!(span > 0) || !Number.isFinite(span)) throw new Error(`Unresolved ISET ${name} extent.`);
  const result = values.map((v, i) => !i ? 0 : i === values.length - 1 ? 1 : (v - values[0]) / span);
  if (result.some((v, i) => i && (duplicate
    ? v < result[i - 1] || (v === result[i - 1]) !== (values[i] === values[i - 1])
    : v <= result[i - 1])))
    throw new Error(`ISET ${name} stations collapse when normalized.`);
  return result;
}

export function locateIsetLeadingEdge(curve, inletSlope) {
  if (!curve || !Array.isArray(curve.knots) || typeof curve.evaluate !== 'function' || !Number.isFinite(inletSlope))
    throw new Error('ISET leading-edge location needs a curve and a finite inlet slope.');
  const tangentProjection = s => {
    const d = curve.evaluate(s).derivative;
    return d.x + inletSlope * d.y;
  };
  // SPLNIT interpolates a sign change in sampled tangent projection; it
  // does not perform a nonlinear root solve or locate panel stagnation.
  for (let i = 1; i < curve.knots.length; i++) {
    const a = curve.knots[i - 1], b = curve.knots[i], da = tangentProjection(a), db = tangentProjection(b);
    if (da < 0 && db >= 0) {
      const parameter = a + (b - a) * da / (da - db);
      return { parameter, point: curve.evaluate(parameter).point, bracket: [a, b],
        tangentProjection: tangentProjection(parameter), method: 'linear interpolation of sampled tangent projections' };
    }
  }
  throw new Error('ISET did not find the incoming normal point on the rounded leading edge.');
}

export function createIsetAirfoilInitialGrid({ points, inletStations, outletStations, surfaceStations, transverseStations,
  inletLength = 4, outletLength = 4, width = 4, inletSlope = 0, outletSlope = 0, leadingEdge = 'rounded' } = {}) {
  if (leadingEdge !== 'rounded') throw new Error('This ISET initializer currently requires a rounded leading edge.');
  if (!Array.isArray(points) || points.length < 4 || !points.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))
    || ![inletLength, outletLength, width, inletSlope, outletSlope].every(Number.isFinite)
    || Math.min(inletLength, outletLength, width) <= 0)
    throw new Error('Invalid ISET airfoil geometry or domain controls.');
  if (points[0].x !== points.at(-1).x || points[0].y !== points.at(-1).y)
    throw new Error('This ISET initializer needs an exactly closed trailing edge; the finite-TE branch is not implemented.');
  // NORMIT uses extrema of the supplied polygon, not the spline or a rotated
  // chord. y's origin is the y coordinate of the first minimum-x sample.
  const nose = points.reduce((a, p) => p.x < a.x ? p : a), xmax = Math.max(...points.map(p => p.x));
  const chord = xmax - nose.x;
  if (!(chord > 0) || !Number.isFinite(chord)) throw new Error('Unresolved ISET chord.');
  const normalizedPoints = points.map(p => ({ x: (p.x - nose.x) / chord, y: (p.y - nose.y) / chord }));
  if (!(signedArea(normalizedPoints) > 0))
    throw new Error('ISET expects positive contour orientation: suction TE to LE to pressure TE.');
  const curve = createContourCurve(normalizedPoints), leading = locateIsetLeadingEdge(curve, inletSlope);
  const inlet = normalizeStations(inletStations, 'inlet'), outlet = normalizeStations(outletStations, 'outlet');
  if (!Array.isArray(surfaceStations?.upper) || !Array.isArray(surfaceStations?.lower))
    throw new Error('Supply both ISET surface station arrays.');
  const bodyCount = Math.min(surfaceStations.upper.length, surfaceStations.lower.length);
  if (bodyCount < 3) throw new Error('ISET needs at least three stations on each surface.');
  // READIN retains the shorter count and NORMIT independently normalizes
  // those prefixes. Record any discarded points instead of hiding them.
  const surface = Object.fromEntries(['upper', 'lower'].map(side => [side,
    normalizeStations(surfaceStations[side].slice(0, bodyCount), `${side} surface`)]));
  const eta = normalizeStations(transverseStations, 'transverse', true);
  const breaks = eta.flatMap((v, j) => j && v === eta[j - 1] ? [j - 1] : []);
  if (breaks.length !== 1) throw new Error('An isolated ISET airfoil needs exactly one repeated transverse station for the farfield break.');
  const split = breaks[0];
  if (split < 2 || eta.length - split - 2 < 2 || !(eta[split] > 0 && eta[split] < 1))
    throw new Error('Each ISET exterior passage needs at least two positive mass intervals.');
  const leadingIndex = inlet.length - 1, trailingIndex = leadingIndex + bodyCount - 1;
  const nx = trailingIndex + outlet.length - 1, te = curve.evaluate(0).point;
  const sLE = leading.parameter, outline = { upper: [], lower: [] };
  for (const side of ['upper', 'lower']) {
    const upstream = inlet.slice(0, -1).map(v => {
      const dx = inletLength * (v - 1);
      return { x: leading.point.x + dx, y: leading.point.y + inletSlope * dx };
    });
    const wall = surface[side].map(v => curve.evaluate(side === 'upper'
      ? sLE * (1 - v) : sLE + (curve.length - sLE) * v).point);
    const wake = outlet.slice(1).map(v => ({ x: te.x + outletLength * v, y: te.y + outletSlope * outletLength * v }));
    outline[side] = [...upstream, ...wall, ...wake];
  }
  // NORLIN's outer x values are uniform in station INDEX. Its two y rows
  // inherit the suction cut's endpoint heights and are separated by width.
  // They need not be horizontal, nor share exact x endpoints with the cuts.
  const farUpper = Array.from({ length: nx + 1 }, (_, i) => {
    const t = i / nx;
    return { x: -inletLength + (1 + inletLength + outletLength) * t,
      y: (1 - t) * outline.upper[0].y + t * outline.upper[nx].y + width * eta[split] };
  });
  const farLower = farUpper.map(p => ({ x: p.x, y: p.y - width }));
  const blend = (a, b, t) => !t ? { ...a } : t === 1 ? { ...b }
    : { x: (1 - t) * a.x + t * b.x, y: (1 - t) * a.y + t * b.y };
  const upperEta = eta.slice(0, split + 1).map(v => v / eta[split]);
  const lowerEta = eta.slice(split + 1).map(v => (v - eta[split]) / (1 - eta[split]));
  const upper = outline.upper.map((p, i) => upperEta.map(t => blend(p, farUpper[i], t)));
  const lower = outline.lower.map((p, i) => lowerEta.map(t => blend(farLower[i], p, t)));
  // RESPLI precedes ELLIP: take the actual grid's suction TE→LE samples,
  // then pressure LE→TE samples, recompute polygon length and respline.
  // A repeated TE remains two independent natural endpoints (p278).
  const resampledPoints = [...outline.upper.slice(leadingIndex, trailingIndex + 1).toReversed(),
    ...outline.lower.slice(leadingIndex + 1, trailingIndex + 1)].map(p => ({ ...p }));
  const resampledCurve = createContourCurve(resampledPoints), resampledLeading = resampledCurve.knots[bodyCount - 1];
  const resampledFractions = { upper: [], lower: [] };
  for (let k = 0; k < bodyCount; k++) {
    resampledFractions.upper.push((resampledLeading - resampledCurve.knots[bodyCount - 1 - k]) / resampledLeading);
    resampledFractions.lower.push((resampledCurve.knots[bodyCount - 1 + k] - resampledLeading) / (resampledCurve.length - resampledLeading));
  }
  const upperMass = eta.slice(1, split + 1).map((v, j) => v - eta[j]);
  const lowerMass = eta.slice(split + 2).map((v, j) => v - eta[split + 1 + j]);
  return { nodes: [lower, upper], massFlows: [lowerMass, upperMass],
    outline, farfield: { lower: farLower, upper: farUpper }, leadingIndex, trailingIndex,
    xi: gilesStreamwiseCoordinates(outline.upper),
    curve, leading, resampled: { points: resampledPoints, curve: resampledCurve, leadingParameter: resampledLeading, surfaceFractions: resampledFractions },
    distributions: { inlet, outlet, surface, transverse: eta, upperEta, lowerEta },
    diagnostics: { stages: ['NORMIT', 'SPLNIT', 'OUTLIN', 'NORLIN', 'RESPLI'],
      originalFortranExecuted: false, flowSolved: false, physicalAcceptance: false,
      geometryNormalization: { origin: { x: nose.x, y: nose.y }, chord },
      surfacePrefixTruncation: Object.fromEntries(['upper', 'lower'].map(side => [side, surfaceStations[side].length - bodyCount])),
      originalTransverseBreak: split, regionOrder: ['lower farfield to pressure outline', 'suction outline to upper farfield'],
      curveParameter: 'cumulative input polygon length; natural cubic, independent TE endpoints',
      farfieldBoundary: 'giles-indexed-y',
      limitations: 'Rounded LE and closed TE, isolated airfoil only. No cascade, finite-TE flap, BL displacement or flow initialization.' } };
}
