// SPDX-License-Identifier: GPL-2.0-or-later
// Read-only physical dead-air gap at an existing BL wake station. Decoded
// final states are authoritative. Live callbacks carry only static BL IDs;
// reconstruct their gap with exactly the production BL arc convention.
import { createSurfaceContourCurve } from '../geometry/contour-topology.js';
import { streamtubeBaseGeometry } from '../euler/streamtube-geometry.js';
import { createXfoilDeadAirGap } from '../viscous/xfoil-dead-air-gap.js';

export function quadCoupledWakeGap({ body, bodyIndex, flow, station, solverLength }) {
  if (station?.wakeGap !== undefined) {
    if (!Number.isFinite(station.wakeGap) || station.wakeGap < 0 || !(solverLength > 0) || !Number.isFinite(solverLength))
      throw new Error('Invalid decoded wake gap or physical length normalization.');
    return station.wakeGap * solverLength;
  }
  if (body?.trailingEdge?.kind !== 'finite-base') return 0;
  if (!Number.isInteger(bodyIndex) || bodyIndex < 0 || !Number.isInteger(station?.i) || station.i < body.trailingIndex)
    throw new Error('Incomplete finite-base wake station mapping.');
  const curve = createSurfaceContourCurve(body.points, body);
  const base = streamtubeBaseGeometry([body], [curve])[0];
  let previous = base.center, distance = 0;
  // Production BL geometry assigns s_wake=0 to the TE station, starts at
  // the solid TE midpoint, and sums current bank-center legs from TE+1.
  for (let i = body.trailingIndex + 1; i <= station.i; i++) {
    const lower = flow?.nodes?.[bodyIndex]?.[i]?.at(-1), upper = flow?.nodes?.[bodyIndex + 1]?.[i]?.[0];
    if (![lower?.x, lower?.y, upper?.x, upper?.y].every(Number.isFinite))
      throw new Error('Missing current finite-base wake geometry.');
    const point = { x: .5 * (lower.x + upper.x), y: .5 * (lower.y + upper.y) };
    const length = Math.hypot(point.x - previous.x, point.y - previous.y);
    if (!(length > 0)) throw new Error('Collapsed current finite-base wake interval.');
    distance += length; previous = point;
  }
  return createXfoilDeadAirGap({ normalGap: base.width, upperDerivative: base.upperDerivative,
    lowerDerivative: base.lowerDerivative }).at(distance).gap;
}
