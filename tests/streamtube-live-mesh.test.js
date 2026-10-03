import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('every accepted Euler iteration publishes its actual updated grid before solve completion', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
  const snapshots = [], events = [];
  const result = solveStreamtubeBody(system, { maxIterations: 3, stepMethod: 'dogleg',
    onIteration: h => events.push(`iteration ${h.iteration}`),
    onMesh: state => {
      events.push(`mesh ${state.iteration.iteration}`);
      snapshots.push(streamtubeMeshSnapshot(state));
    } });
  assert.equal(snapshots.length, 3);
  assert.deepEqual(events, ['iteration 0', 'iteration 1', 'mesh 1', 'iteration 2', 'mesh 2', 'iteration 3', 'mesh 3']);
  for (let k = 0; k < snapshots.length; k++) {
    assert.deepEqual(snapshots[k].iteration, result.history[k + 1]);
    assert.equal(snapshots[k].initialization.flowSolved, false);
    if (k) { assert.notDeepEqual(snapshots[k].vertices, snapshots[k - 1].vertices); assert.deepEqual(snapshots[k].cells, snapshots[k - 1].cells); }
  }
  assert.deepEqual(snapshots.at(-1).vertices, streamtubeMeshSnapshot({ system, nodes: result.nodes }).vertices);
});
