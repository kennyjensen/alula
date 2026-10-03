// SPDX-License-Identifier: GPL-2.0-or-later
// Explicit solid-boundary topology. A finite base is a real solid polyline,
// separate from the smooth upper-TE -> nose -> lower-TE surface chain.
// No geometric corner inference, source-point removal or base straightening.
import { prepareContour, signedArea } from './airfoil.js';
import { createContourCurve } from './contour-curve.js';

const copy = p => ({ x: p.x, y: p.y });
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const dot = (a, b) => a.x * b.x + a.y * b.y;
const norm = a => Math.hypot(a.x, a.y);
const right = a => ({ x: a.y, y: -a.x });
function unit(a) {
  const length = norm(a);
  if (!(length > 0) || !Number.isFinite(length)) throw new Error('Degenerate trailing-edge direction.');
  return { x: a.x / length, y: a.y / length };
}
function chain(points, indices, panelCount) {
  const selected = indices.map(i => copy(points[i]));
  const panels = selected.slice(1).map((end, k) => {
    const start = selected[k], difference = sub(end, start), tangent = unit(difference);
    return { sourcePanelIndex: indices[k] % panelCount, start: copy(start), end: copy(end),
      length: norm(difference), tangent, outwardNormal: right(tangent) };
  });
  return { points: selected, indices: indices.slice(), panels,
    length: panels.reduce((sum, p) => sum + p.length, 0) };
}

export function createContourTopology(input, { trailingEdge } = {}) {
  if (trailingEdge !== undefined && (!trailingEdge || !['sharp', 'finite-base'].includes(trailingEdge.kind)))
    throw new Error('Supply explicit sharp or finite-base trailing-edge topology.');
  const finite = trailingEdge?.kind === 'finite-base';
  if (finite) {
    if (!Array.isArray(input) || input.length < 9 || !input.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))
      throw new Error('Finite-base topology requires a finite closed solid contour.');
    if (input[0].x !== input.at(-1).x || input[0].y !== input.at(-1).y)
      throw new Error('Finite-base topology requires the retained repeated closing point.');
    // Explicit indices refer to the supplied CCW source contour. Reversing it
    // silently would reinterpret its physical corner and base designations.
    if (!(signedArea(input) > 0)) throw new Error('Finite-base source indices require a counterclockwise contour.');
  }
  const points = prepareContour(input, finite ? { lifting: false } : undefined);
  const n = points.length - 1;
  let upperIndex = 0, lowerIndex = n, surfaceIndices = Array.from({ length: n + 1 }, (_, i) => i), baseIndices = [n];
  if (finite) {
    ({ upperIndex, lowerIndex } = trailingEdge);
    if (![upperIndex, lowerIndex].every(i => Number.isInteger(i) && i >= 0 && i < n) || upperIndex === lowerIndex)
      throw new Error('Finite-base upper/lower corners require distinct indices of original nonrepeated vertices.');
    const walk = (from, to) => Array.from({ length: (to - from + n) % n + 1 }, (_, k) => {
      const i = (from + k) % n;
      return i === 0 && k > 0 ? n : i;
    });
    surfaceIndices = walk(upperIndex, lowerIndex); baseIndices = walk(lowerIndex, upperIndex);
    if (surfaceIndices.length < 4) throw new Error('Finite-base surface chain needs at least four points.');
  }
  const surface = chain(points, surfaceIndices, n), base = chain(points, baseIndices, n);
  const first = surface.panels[0], last = surface.panels.at(-1);
  const upper = { index: upperIndex, point: copy(surface.points[0]), contourTangent: copy(first.tangent),
    downstreamTangent: { x: -first.tangent.x, y: -first.tangent.y }, outwardNormal: copy(first.outwardNormal) };
  const lower = { index: lowerIndex, point: copy(surface.points.at(-1)), contourTangent: copy(last.tangent),
    downstreamTangent: copy(last.tangent), outwardNormal: copy(last.outwardNormal) };
  const wakeTangent = unit({ x: upper.downstreamTangent.x + lower.downstreamTangent.x,
    y: upper.downstreamTangent.y + lower.downstreamTangent.y });
  const wakeNormal = { x: -wakeTangent.y, y: wakeTangent.x }, gapVector = sub(upper.point, lower.point);
  const gapLength = norm(gapVector), normalGap = dot(gapVector, wakeNormal);
  if (finite && (!(gapLength > 0) || !(normalGap > 0)))
    throw new Error('Finite-base corner ordering must give a positive gap across the downstream surface bisector.');
  return { kind: finite ? 'finite-base' : 'sharp', points, surface, base,
    sourcePanelCount: n, orientation: 'counterclockwise',
    units: { coordinates: 'input length units', lengths: 'input length units', directions: 'unit vectors' },
    indexing: finite ? 'original closed CCW input; closing vertex retained'
      : 'prepared sharp contour; existing orientation normalization retained',
    trailingEdge: { kind: finite ? 'finite-base' : 'sharp', upper, lower,
      center: { x: .5 * (upper.point.x + lower.point.x), y: .5 * (upper.point.y + lower.point.y) },
      gapVector, gapLength, normalGap, tangentialOffset: dot(gapVector, wakeTangent), wakeTangent, wakeNormal,
      directionDefinition: 'one-sided source-panel surface tangents; downstream bisector is geometric, not a solved wake direction' } };
}

// Same curve API and branch parameter convention as the existing sharp
// contour. Its two finite-TE branch endpoints are distinct physical corners.
// The retained base belongs to createContourTopology(...).base, never this fit.
export function createSurfaceContourCurve(points, options = {}) {
  const topology = createContourTopology(points, options);
  return createContourCurve(topology.surface.points, { allowOpenEndpoints: topology.kind === 'finite-base' });
}
