import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { streamtubeCellGeometry } from '../src/euler/streamtube-cell.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { refineStreamtubeBody } from '../src/euler/streamtube-refinement.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { requirePositiveSimplePolygon } from '../src/geometry/simple-polygon.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directPolygonGeometry, directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

test('production volume checks agree with independent geometry for bent and crossed conservation polygons', () => {
  const polygons = [
    [[0, 0], [1, 0], [2, 0], [2, 1], [1, .4], [0, 1]],
    [[0, 0], [1, 0], [2, 0], [2, 1], [1, -.1], [0, 1]],
    [[0, 0], [1, 0], [.5, 0], [2, 1], [1, 1], [0, 1]],
  ].map(row => row.map(([x, y]) => ({ x, y })));
  for (const p of polygons) for (const vertices of [p, p.toReversed()]) {
    const oracle = directPolygonGeometry(vertices);
    if (oracle.valid) assert.ok(Math.abs(requirePositiveSimplePolygon(vertices) - oracle.area) < 1e-14);
    else assert.throws(() => requirePositiveSimplePolygon(vertices));
  }
  // Both displayed quads are positive simple polygons, but a midpoint
  // subvolume crosses. This state must fail before any flux is evaluated.
  const lower = [{ x: -1, y: 0 }, { x: 0, y: 0 }, { x: 1, y: 0 }];
  const upper = [{ x: -1, y: 3 }, { x: .6, y: 3 }, { x: .5, y: 4 }];
  assert.throws(() => streamtubeCellGeometry(lower, upper, { geometryDomain: 'positive-simple' }), /half-volume/);
});

test('the conservative retained root remains admissible while the crossed isentropic root is rejected', () => {
  const rows = JSON.parse(readFileSync(new URL('../docs/default-euler-ises-sampled.json', import.meta.url))).cases;
  for (const row of rows) {
    const { input, initialEuler } = row.restart, system = createStreamtubeBodySystem(input);
    const state = system.adoptGeometry(Float64Array.from(initialEuler.x), initialEuler.nodes);
    if (row.streamwiseMode === 'momentum') {
      assert.equal(directStreamtubeVolumeGeometry(initialEuler.nodes).valid, true);
      assert.ok(system.evaluate(state).diagnostics.residual < 1e-10);
    } else {
      assert.equal(directStreamtubeVolumeGeometry(initialEuler.nodes).valid, false);
      assert.throws(() => system.evaluate(state), /half-volume|conservation-volume/);
    }
  }
});

test('common scalar retry accepts a valid refined step and zero retries retain the last accepted state', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'walls-only', normalStencil: 'body-stations' };
  const source = createStreamtubeBodySystem(input), parent = solveStreamtubeBody(source, { maxIterations: 8 });
  assert.equal(parent.converged, true);
  const fine = refineStreamtubeBody(input, source, { initial: parent.x });
  const start = solveStreamtubeIses(fine.input, { initialEuler: fine.initialEuler, maxIterations: 0 });
  const blocked = solveStreamtubeIses(fine.input, { initialEuler: fine.initialEuler, maxIterations: 1, stepAcceptance: 'admissible', maxBacktracks: 0 });
  assert.equal(blocked.converged, false); assert.equal(blocked.history.length, 1);
  assert.deepEqual(blocked.x, start.x); assert.deepEqual(blocked.nodes, start.nodes);
  const meshes = [], retried = solveStreamtubeIses(fine.input, { initialEuler: fine.initialEuler, maxIterations: 1,
    stepAcceptance: 'admissible', onMesh: m => meshes.push(m.nodes) });
  assert.equal(retried.history.length, 2, retried.reason);
  assert.ok(retried.history[1].backtracks > 0);
  assert.equal(retried.linearDiagnostics.solves, 1, 'a retry must keep the same Newton direction');
  assert.equal(meshes.length, 2, 'publish accepted grids only');
  for (const nodes of meshes) assert.equal(directStreamtubeVolumeGeometry(nodes).valid, true);
  assert.equal(retried.converged, false, 'one valid update is not a solution');
  assert.throws(() => solveStreamtubeIses(input, { stepAcceptance: 'unknown' }), /controls/);
});
