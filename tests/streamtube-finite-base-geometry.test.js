import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createSurfaceContourCurve } from '../src/geometry/contour-topology.js';
import { streamtubeMeshConnectivity } from '../src/geometry/streamtube-mesh-connectivity.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

// Manufactured thin finite-base shape. This deliberately modifies a NACA
// fixture; it does not alter or substitute for either measured NLR element.
function fixture() {
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 2, contourPanels: 40 });
  const body = input.bodies[0], lowerIndex = body.points.length - 1;
  const surface = body.points.map((p, i) => ({ x: p.x,
    y: p.y + (i < lowerIndex / 2 ? 1 : -1) * .0005 * p.x }));
  body.points = [...surface, { x: 1, y: 0 }, { ...surface[0] }];
  body.trailingEdge = { kind: 'finite-base', upperIndex: 0, lowerIndex };
  body.stagnationParameter = createSurfaceContourCurve(body.points, body).length / 2;
  input.displacement = {
    surfaces: [{ upper: Array(5).fill(.0001), lower: Array(5).fill(.0001) }],
    wakes: [Array(input.outerLower.length - 1 - body.trailingIndex).fill(.0012)],
  };
  return input;
}
const error = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);

test('finite-base Euler geometry preserves both solid TE corners and separate physical connectivity', () => {
  const input = fixture(), before = structuredClone(input), system = createStreamtubeBodySystem(input);
  const { nodes, undisplacedNodes } = system.decode(system.initial), b = input.bodies[0], te = b.trailingIndex;
  assert.deepEqual(input, before);
  assert.deepEqual(undisplacedNodes[1][te][0], b.points[0]);
  assert.deepEqual(undisplacedNodes[0][te].at(-1), b.points[b.trailingEdge.lowerIndex]);
  assert.equal(system.curves[0].knots.length, b.trailingEdge.lowerIndex + 1);
  const mesh = streamtubeMeshSnapshot({ system, nodes });
  assert.equal(mesh.quality.valid, true);
  const { indices } = streamtubeMeshConnectivity(system.layout, nodes);
  assert.notEqual(indices[0][te].at(-1), indices[1][te][0]);
  for (let i = te + 1; i < nodes[0].length; i++)
    assert.ok(Math.abs(error(nodes[0][i].at(-1), nodes[1][i][0]) - .0012) < 1e-15);
  const adopted = system.adoptGeometry(system.initial, nodes);
  assert.deepEqual(system.decode(adopted).nodes, nodes);
});

test('undisplaced finite TE vertices remain separate while the leading edge is shared', () => {
  const input = fixture(); delete input.displacement;
  const system = createStreamtubeBodySystem(input), nodes = system.decode(system.initial).nodes;
  const { indices } = streamtubeMeshConnectivity(system.layout, nodes), b = input.bodies[0];
  assert.equal(indices[0][b.leadingIndex].at(-1), indices[1][b.leadingIndex][0]);
  assert.notEqual(indices[0][b.trailingIndex].at(-1), indices[1][b.trailingIndex][0]);
  assert.equal(streamtubeMeshSnapshot({ system, nodes }).quality.valid, true);
});

test('complete prescribed-gap Euler Jacobian includes finite TE surface motion and both wake banks', () => {
  const input = fixture(), system = createStreamtubeBodySystem(input);
  const state = system.initial.map((v, i) => v + 1e-6 * Math.sin(i));
  const jacobian = system.jacobian(state), h = 1e-7;
  for (const phase of [.3, 1.1, 2.4]) {
    const direction = state.map((_, i) => Math.sin((i + 1) * phase));
    const plus = system.residual(state.map((v, i) => v + h * direction[i]));
    const minus = system.residual(state.map((v, i) => v - h * direction[i]));
    for (let row = 0; row < system.layout.n; row++) {
      let exact = 0;
      for (let col = 0; col < system.layout.n; col++) exact += jacobian[row * system.layout.n + col] * direction[col];
      const fd = (plus[row] - minus[row]) / (2 * h);
      assert.ok(Math.abs(fd - exact) < 3e-6 * Math.max(1, Math.abs(exact)), `row ${row}: ${fd} != ${exact}`);
    }
  }
});
