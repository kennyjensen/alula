// SPDX-License-Identifier: GPL-2.0-or-later
// Standard airfoil files supply an open wetted TE-to-TE chain when the TE
// has finite thickness. Close only its implied straight base for solvers;
// retain the original ordering/coordinates for standard-file round trips.
import { prepareContour, signedArea } from './airfoil.js';
import { createContourTopology } from './contour-topology.js';

const copy = points => points.map(p => ({ x: p.x, y: p.y }));
const same = (a, b) => a.x === b.x && a.y === b.y;

export function prepareAirfoilElement(element) {
  if (!element || !Array.isArray(element.points) || element.points.length < 9
    || !element.points.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)))
    throw new Error('Airfoil elements require at least nine finite supplied coordinates.');
  if (element.sourcePoints !== undefined && (!Array.isArray(element.sourcePoints)
    || !element.sourcePoints.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y))))
    throw new Error('Invalid original airfoil source coordinates.');
  const points = copy(element.points), result = { ...element,
    ...(element.sourcePoints === undefined ? {} : { sourcePoints: copy(element.sourcePoints) }) };
  if (points.some((p, i) => i > 0 && same(p, points[i - 1])))
    throw new Error('Repeated successive airfoil coordinates mark surface corners; this solver does not yet support those corner/split-spline conditions.');
  if (element.trailingEdge?.kind === 'finite-base') {
    // Explicit indices identify a measured base. Do not infer a different
    // pair of corners or discard any point of that retained solid polyline.
    const topology = createContourTopology(points, { trailingEdge: element.trailingEdge });
    return { ...result, points: topology.points, trailingEdge: { ...element.trailingEdge } };
  }
  if (element.trailingEdge !== undefined && element.trailingEdge?.kind !== 'sharp')
    throw new Error('Unknown explicit airfoil trailing-edge topology.');
  if (same(points[0], points.at(-1))) {
    // Keep the existing sharp preparation/output shape. A closed polyline
    // alone cannot tell us where a measured finite base begins or ends.
    return { ...result, points: prepareContour(points),
      ...(element.trailingEdge === undefined ? {} : { trailingEdge: { ...element.trailingEdge } }) };
  }
  if (element.trailingEdge?.kind === 'sharp')
    throw new Error('Explicit sharp trailing edges require coincident supplied endpoints; distinct endpoints are not snapped together.');
  const wetted = signedArea(points) < 0 ? points.slice().reverse() : points;
  const closed = [...wetted, { ...wetted[0] }];
  const trailingEdge = { kind: 'finite-base', upperIndex: 0, lowerIndex: wetted.length - 1 };
  const topology = createContourTopology(closed, { trailingEdge });
  return { ...result, points: topology.points, trailingEdge,
    sourcePoints: copy(element.sourcePoints ?? points) };
}
