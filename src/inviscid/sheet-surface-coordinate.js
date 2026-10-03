// SPDX-License-Identifier: GPL-2.0-or-later
import { prepareContour } from '../geometry/airfoil.js';
import { createContourTopology } from '../geometry/contour-topology.js';

const same = (a, b) => a.x === b.x && a.y === b.y;
const finitePoint = p => p && Number.isFinite(p.x) && Number.isFinite(p.y);

// XFOIL GGCALC, third_party/Xfoil/src/xpanel.f:1113–1117, assigns
// QINVU = GAMU; QISET:1636–1650 forms the surface tangential speeds.
// Our CCW panel convention has q_exterior - q_interior = gamma, so using
// gamma as surface speed is the zero-interior-velocity approximation.
//
// This is a material surface COORDINATE: integrate linear nodal gamma on
// the original cumulative panel-length parameter, then pull that scalar
// back to the C2 wall using the same parameter. It is not the panel field's
// physical potential at a point on a different spline boundary, and it
// does not extend a vortex-sheet field through the solid. No geometry,
// nodal strength or circulation is modified here. The optional conservative
// reconstruction changes only unresolved within-panel coordinate variation.
export function createSheetSurfaceCoordinate({ points, knots, field, element = 0, stagnationPotential = 0, trailingEdge,
  reconstruction = 'nodal-linear' } = {}) {
  if (!['nodal-linear', 'conservative-panel'].includes(reconstruction)) throw new Error('Unknown sheet-potential reconstruction.');
  const topology = trailingEdge?.kind === 'finite-base' ? createContourTopology(points, { trailingEdge }) : null;
  if (topology) points = topology.surface.points;
  if (!Array.isArray(points) || !points.every(finitePoint)
    || !Array.isArray(knots) && !(knots instanceof Float64Array)
    || knots.length !== points.length || !Number.isInteger(element) || element < 0
    || !Number.isFinite(stagnationPotential)) throw new Error('Invalid sheet-surface coordinate inputs.');
  // Reuse the existing contour validity rules, but refuse normalization:
  // the ordered, exactly closed CCW contour must already match the field.
  const prepared = topology ? topology.surface.points : prepareContour(points);
  if (!points.every((p, i) => same(p, prepared[i])))
    throw new Error('Sheet-surface points must already be exactly closed and counterclockwise.');
  points = points.map(p => ({ x: p.x, y: p.y }));
  knots = Array.from(knots);
  if (knots[0] !== 0) throw new Error('Sheet-surface knots must start at zero.');
  let length = 0;
  const lengths = [];
  for (let i = 1; i < points.length; i++) {
    const distance = Math.hypot(points[i].x - points[i - 1].x, points[i].y - points[i - 1].y);
    lengths.push(distance); length += distance;
    if (!Number.isFinite(knots[i]) || !(knots[i] > knots[i - 1]) || knots[i] !== length)
      throw new Error('Sheet-surface knots must equal the original cumulative panel lengths.');
  }
  if (!field || !Array.isArray(field.panels)
    || !Array.isArray(field.gamma) && !(field.gamma instanceof Float64Array)
    || !Array.from(field.gamma).every(Number.isFinite)) throw new Error('Invalid linear-vortex sheet field.');
  const panels = field.panels.filter(p => p?.element === element);
  if (panels.length !== points.length - 1 || !Number.isInteger(panels[0]?.node) || panels[0].node < 0)
    throw new Error('Sheet-surface element does not match the ordered panel group.');
  const start = panels[0].node, gamma = Array.from(field.gamma).slice(start, start + points.length);
  if (gamma.length !== points.length) throw new Error('Sheet-surface strengths do not cover every independent TE-to-TE node.');
  panels.forEach((p, i) => {
    const tx = (points[i + 1].x - points[i].x) / lengths[i];
    const ty = (points[i + 1].y - points[i].y) / lengths[i];
    if (!finitePoint(p.a) || !finitePoint(p.b) || !same(p.a, points[i]) || !same(p.b, points[i + 1])
      || p.node !== start + i || p.length !== lengths[i] || p.tx !== tx || p.ty !== ty)
      throw new Error('Sheet-surface panel geometry, orientation or nodal ordering does not match its contour.');
  });

  // A linear segment can cross zero only once. A zero interval has no
  // unique stagnation location and no invertible potential coordinate.
  const roots = [];
  for (let i = 0; i < gamma.length - 1; i++) {
    if (gamma[i] === 0 && gamma[i + 1] === 0)
      throw new Error('A flat zero-strength segment makes the sheet-surface coordinate undefined.');
    if (gamma[i] < 0 && gamma[i + 1] > 0) {
      // Normalization avoids overflow for otherwise finite strengths.
      const scale = Math.max(-gamma[i], gamma[i + 1]);
      const fraction = (-gamma[i] / scale) / ((-gamma[i] / scale) + gamma[i + 1] / scale);
      roots.push({ parameter: knots[i] + fraction * (knots[i + 1] - knots[i]), segment: i, fraction, kind: 'interior-linear-zero' });
    }
    if (i > 0 && gamma[i] === 0 && gamma[i - 1] < 0 && gamma[i + 1] > 0)
      roots.push({ parameter: knots[i], segment: i, fraction: 0, kind: 'nodal-zero' });
  }
  const linearPanels = new Set();
  let root = roots[0];
  if (reconstruction === 'conservative-panel') {
    // A nodal interpolation undershoot can reverse the derivative inside
    // a panel whose measured potential increment is strictly forward. Keep
    // that integral and the incoming-stagnation quadratic; use its secant
    // only on such unresolved intervals. Never erase a reversed panel flux.
    const candidates = roots.filter(candidate => gamma[0] <= 0 && gamma.at(-1) >= 0
      && lengths.every((_, i) => {
        if (i === candidate.segment || candidate.fraction === 0 && i === candidate.segment - 1) return true;
        const mean = .5 * gamma[i] + .5 * gamma[i + 1];
        return knots[i + 1] <= candidate.parameter ? mean < 0 : mean > 0;
      }));
    if (candidates.length !== 1) throw new Error('Panel potential increments do not define one resolved incoming branch.');
    root = candidates[0];
    lengths.forEach((_, i) => {
      if (i === root.segment || root.fraction === 0 && i === root.segment - 1) return;
      if (knots[i + 1] <= root.parameter ? gamma[i] > 0 || gamma[i + 1] > 0 : gamma[i] < 0 || gamma[i + 1] < 0)
        linearPanels.add(i);
    });
  } else if (roots.length !== 1) throw Object.assign(new Error('Sheet-surface coordinate requires exactly one incoming negative-to-positive strength zero.'),
    { code: 'SHEET_SURFACE_NONMONOTONE' });
  const stagnationParameter = root.parameter;
  if (!(stagnationParameter > 0 && stagnationParameter < length))
    throw new Error('Sheet-surface incoming zero is outside resolvable interior parameter space.');
  if (root.fraction !== 0 && !(stagnationParameter > knots[root.segment] && stagnationParameter < knots[root.segment + 1]))
    throw new Error('Sheet-surface incoming zero is below the parameter resolution.');
  if (reconstruction === 'nodal-linear' && gamma.some((q, i) => knots[i] < stagnationParameter ? q > 0 : knots[i] > stagnationParameter && q < 0))
    throw Object.assign(new Error('Sheet-surface strengths are not monotone in potential on both stagnation-to-TE branches.'),
      { code: 'SHEET_SURFACE_NONMONOTONE' });

  const slopes = lengths.map((_, i) => (gamma[i + 1] - gamma[i]) / (knots[i + 1] - knots[i]));
  if (!slopes.every(Number.isFinite)) throw new Error('Sheet-surface strength slope is not finite.');
  const stagnationSlopes = root.fraction === 0 ? { upper: slopes[root.segment - 1], lower: slopes[root.segment] }
    : { upper: slopes[root.segment], lower: slopes[root.segment] };
  const stagnationSlope = .5 * stagnationSlopes.upper + .5 * stagnationSlopes.lower;
  const segmentIntegral = (i, a, b) => {
    const h = knots[i + 1] - knots[i];
    if (linearPanels.has(i)) return (b - a) * (.5 * gamma[i] + .5 * gamma[i + 1]);
    const ga = gamma[i] + (a - knots[i]) / h * (gamma[i + 1] - gamma[i]);
    const gb = gamma[i] + (b - knots[i]) / h * (gamma[i + 1] - gamma[i]);
    return (b - a) * (.5 * ga + .5 * gb);
  };
  // Integrate outward from the actual zero. This retains the O(ds^2)
  // stagnation behavior without subtracting two O(chord) primitives.
  const relative = new Float64Array(knots.length), r = root.segment;
  let rightIndex;
  if (root.fraction === 0) {
    relative[r] = 0; rightIndex = r;
  } else {
    relative[r] = .5 * slopes[r] * (knots[r] - stagnationParameter) ** 2;
    relative[r + 1] = .5 * slopes[r] * (knots[r + 1] - stagnationParameter) ** 2;
    rightIndex = r + 1;
  }
  const accumulate = (begin, end, direction) => {
    let sum = relative[begin], correction = 0;
    for (let i = begin + direction; direction > 0 ? i <= end : i >= end; i += direction) {
      const previous = i - direction, interval = Math.min(i, previous);
      const term = segmentIntegral(interval, knots[previous], knots[i]) - correction;
      const next = sum + term; correction = (next - sum) - term;
      relative[i] = sum = next;
    }
  };
  accumulate(r, 0, -1); accumulate(rightIndex, knots.length - 1, 1);
  if (!Array.from(relative).every(Number.isFinite)) throw new Error('Sheet-surface primitive is not finite.');
  const locate = s => {
    if (!Number.isFinite(s) || s < 0 || s > length) throw new Error('Sheet-surface parameter is outside the contour.');
    let lo = 0, hi = knots.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (knots[mid] > s) hi = mid; else lo = mid; }
    return lo;
  };
  const derivative = s => {
    const i = locate(s);
    if (s === stagnationParameter) return 0;
    if (linearPanels.has(i)) return .5 * gamma[i] + .5 * gamma[i + 1];
    if (s === length) return gamma.at(-1);
    if (s === knots[i]) return gamma[i];
    if (i === r || root.fraction === 0 && i === r - 1) return slopes[i] * (s - stagnationParameter);
    const f = (s - knots[i]) / (knots[i + 1] - knots[i]);
    return (1 - f) * gamma[i] + f * gamma[i + 1];
  };
  const phase = s => {
    const i = locate(s);
    if (s === stagnationParameter) return stagnationPotential;
    let value;
    if (s === length) value = relative.at(-1);
    else if (s === knots[i]) value = relative[i];
    else if (i === r || root.fraction === 0 && i === r - 1) value = .5 * slopes[i] * (s - stagnationParameter) ** 2;
    else value = relative[i] + segmentIntegral(i, knots[i], s);
    const result = stagnationPotential + value;
    if (!Number.isFinite(result)) throw new Error('Sheet-surface potential gauge overflows.');
    return result;
  };
  let circulation = 0, correction = 0;
  for (let i = 0; i < lengths.length; i++) {
    const term = lengths[i] * (.5 * gamma[i] + .5 * gamma[i + 1]) - correction;
    const next = circulation + term; correction = (next - circulation) - term; circulation = next;
  }
  const contourPotentialIncrement = relative.at(-1) - relative[0];
  if (![circulation, contourPotentialIncrement].every(Number.isFinite)) throw new Error('Sheet-surface circulation is not finite.');
  return { length, stagnationParameter, phase, derivative,
    diagnostics: {
      method: 'integrated nodal-sheet-strength material surface coordinate',
      ...(linearPanels.size ? { reconstruction: 'conservative-panel', reconstructedPanels: [...linearPanels].map(i => ({ index: i,
        from: knots[i], to: knots[i + 1], originalEndpointStrengths: [gamma[i], gamma[i + 1]],
        potentialIncrement: lengths[i] * (.5 * gamma[i] + .5 * gamma[i + 1]) })),
        reconstructionScope: 'Initializer coordinate only: original panel potential increments, strengths, circulation and geometry retained; unresolved within-panel reversal removed.' } : {}),
      interpretation: 'Initializer coordinate on unchanged chord-length knots; not physical panel potential on a different C2 wall.',
      source: 'third_party/Xfoil/src/xpanel.f:1113-1117 and 1636-1650', element,
      sourceConvention: 'Surface speed equals nodal gamma (zero interior velocity approximation); not exact off-sheet field velocity.',
      panelCount: panels.length, nodeCount: points.length, stagnationPotential,
      stagnation: { ...root, derivative: 0, strengthSlope: slopes[r] },
      stagnationSlope, stagnationSlopes,
      stagnationSlopeConvention: root.fraction === 0 ? 'Arithmetic mean of the two positive one-sided sheet slopes at a node.' : 'Unique linear sheet slope through the interior zero.',
      sampledMonotone: true, monotonicity: linearPanels.size ? 'Strict panel potential increments; exact stagnation quadratic; conservative secants on unresolved intervals.'
        : 'Exact piecewise-linear sign check; no zero-strength interval.',
      upperPotentialIncrement: relative[0], lowerPotentialIncrement: relative.at(-1),
      circulation, contourPotentialIncrement,
      circulationRoundoffDifference: contourPotentialIncrement - circulation
    }
  };
}
