// SPDX-License-Identifier: GPL-2.0-or-later
// Manufactured stacked bodies, not measured airfoils or a solved flow.
// Unequal base widths and curved near-wake centers exercise two active
// XFOIL cubic tails on opposite banks of the same middle passage.
import { intrinsicBodyFixture } from './intrinsic-body.js';
import { createSurfaceContourCurve } from '../../src/geometry/contour-topology.js';
import { createStreamtubeBodySystem } from '../../src/euler/streamtube-body.js';
import { createStreamtubeBoundaryLayers } from '../../src/euler/streamtube-boundary-layers.js';
import { createCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';
import { initialStreamtubeDisplacement, streamtubeBaseGeometry } from '../../src/euler/streamtube-geometry.js';

export function twoActiveFiniteBaseWakes() {
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 2, contourPanels: 40, mach: .03 });
  const source = input.bodies[0];
  input.bodies = [-.25, .25].map((offset, b) => {
    const body = structuredClone(source), lowerIndex = body.points.length - 1, halfGap = [.002, .003][b];
    const wetted = body.points.map((p, i) => ({ x: p.x,
      y: p.y + offset + (i < lowerIndex / 2 ? 1 : -1) * halfGap * p.x }));
    body.points = [...wetted, { x: 1, y: offset }, { ...wetted[0] }];
    body.trailingEdge = { kind: 'finite-base', upperIndex: 0, lowerIndex };
    body.element = b;
    body.stagnationParameter = createSurfaceContourCurve(body.points, body).length / 2;
    return body;
  });
  input.cutPaths = [-.25, .25].map(y => input.outerLower.map(p => ({ x: p.x, y })));
  input.weights = [[1, 2.3], [1.4, 3.1], [2.7, 1]];
  const te = source.trailingIndex;
  for (const row of [input.outerLower, input.outerUpper, ...input.cutPaths])
    [.0015, .004, .007].forEach((dx, k) => { row[te + k + 1].x = 1 + dx; });
  input.wakeGeometry = 'independent-banks'; input.wakeOutlet = 'banks';
  const displacement = initialStreamtubeDisplacement({ bodies: input.bodies, nx: input.outerLower.length - 1 },
    streamtubeBaseGeometry(input.bodies, input.bodies.map(b => createSurfaceContourCurve(b.points, b))));
  const euler = createStreamtubeBodySystem({ ...input, displacement }), state = euler.initial.slice();
  // Each body's two banks move together here, with different curvature.
  // All middle-passage unknowns subsequently remain independent FD columns.
  for (const p of euler.layout.positions) if (p.i > te)
    state[p.column] += .00012 * Math.sin(1.3 * p.i + (p.body ?? p.group ?? 0) * .7);
  const bl = createStreamtubeBoundaryLayers(euler, state), geo = bl.geometry(state);
  const values = new Float64Array(4 * bl.stations.length);
  for (const s of bl.stations) {
    const theta = (s.kind === 'wake' ? .0004 : .0002) * (1 + .15 * s.body);
    const delta = 2.5 * theta + (geo.coordinates[s.id].wakeGap ?? 0);
    values.set([['similarity', 'laminar'].includes(s.regime) ? .2 : .03,
      theta / bl.scale, delta / bl.scale, 1 + .008 * Math.sin(s.id + .3)], 4 * s.id);
  }
  const system = createCoupledStreamtubeBody(input, { initialEuler: { x: state, nodes: euler.decode(state).nodes },
    initialBL: values, reynolds: 1e6, ncrit: 9, edgeMatching: 'section-velocity' });
  return { input, system, x: system.initial.slice() };
}
