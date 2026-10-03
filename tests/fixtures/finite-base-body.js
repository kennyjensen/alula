// SPDX-License-Identifier: GPL-2.0-or-later
// Manufactured symmetric thin base, not an approximation of measured NLR.
import { intrinsicBodyFixture } from './intrinsic-body.js';
import { createSurfaceContourCurve } from '../../src/geometry/contour-topology.js';

export function finiteBaseBodyFixture({ bodySegments = 8, tubes = 3 } = {}) {
  const input = intrinsicBodyFixture({ bodySegments, tubes, contourPanels: 40 });
  const body = input.bodies[0], lowerIndex = body.points.length - 1;
  const surface = body.points.map((p, i) => ({ x: p.x,
    y: p.y + (i < lowerIndex / 2 ? 1 : -1) * .0005 * p.x }));
  body.points = [...surface, { x: 1, y: 0 }, { ...surface[0] }];
  body.trailingEdge = { kind: 'finite-base', upperIndex: 0, lowerIndex };
  body.stagnationParameter = createSurfaceContourCurve(body.points, body).length / 2;
  return input;
}
