import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { directPolygonGeometry, directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

const points = a => a.map(([x, y]) => ({ x, y }));

test('polygon oracle distinguishes concavity, crossed edges, reversal, touching and collinear wall pieces', () => {
  const concave = points([[0, 0], [1, 0], [.49, .49], [0, 1]]);
  const q = directPolygonGeometry(concave);
  assert.equal(q.valid, true); assert.equal(q.area, .49); assert.ok(q.minCornerSine < 0);
  assert.equal(directPolygonGeometry(points([[0, 0], [.5, 0], [1, 0], [1, 1], [.5, 1], [0, 1]])).valid, true);
  for (const p of [concave.toReversed(), points([[0, 0], [2, 1], [0, 1], [1, 0]]),
    points([[0, 0], [1, 0], [.5, 0], [0, 1]]), points([[0, 0], [1, 0], [1, 0], [0, 1]])])
    assert.equal(directPolygonGeometry(p).valid, false);
});

test('positive simple displayed cells can hide a crossed midpoint control volume', () => {
  // Exact rational-coordinate counterexample; independent of an airfoil solve.
  const nodes = [[points([[-1, 0], [-1, 3]]), points([[0, 0], [.6, 3]]), points([[1, 0], [.5, 4]])]];
  const q = directStreamtubeVolumeGeometry(nodes);
  assert.equal(q.primal.failures.length, 0);
  assert.equal(q.embedding.crossings.length, 0);
  assert.equal(q.valid, false); assert.ok(q.halves.failures.length > 0); assert.ok(q.dual.failures.length > 0);
});

test('global edge audit rejects duplicate overlapping regions and partial edges sharing an endpoint', () => {
  const first = [points([[0, 0], [0, 1]]), points([[1, 0], [1, 1]])];
  const duplicate = directStreamtubeVolumeGeometry([first, structuredClone(first)]);
  assert.equal(duplicate.valid, false); assert.ok(duplicate.embedding.nonmanifoldEdges.length > 0);
  const partial = directStreamtubeVolumeGeometry([first, [points([[0, -.5], [0, 0]]), points([[.5, -.5], [.5, 0]])]]);
  assert.equal(partial.primal.failures.length, 0);
  assert.ok(partial.embedding.crossings.some(c => c.type === 'shared-endpoint-overlap'));
});

test('retained default roots distinguish a valid polygon mesh from a crossed flap control volume', () => {
  const saved = JSON.parse(readFileSync(new URL('../docs/default-euler-ises-sampled.json', import.meta.url)));
  for (const row of saved.cases) {
    const { input, initialEuler } = row.restart, q = directStreamtubeVolumeGeometry(initialEuler.nodes);
    assert.equal(q.primal.count, 2919); assert.equal(q.halves.count, 5838); assert.equal(q.dual.count, 2898);
    assert.equal(q.primal.failures.length, 0); assert.equal(q.concavePrimal.length, 1);
    assert.equal(q.embedding.crossings.length, 0); assert.equal(q.embedding.nonmanifoldEdges.length, 0);
    assert.equal(input.bodies[0].element, 1); assert.equal(input.bodies[0].leadingIndex, 58);
    if (row.streamwiseMode === 'momentum') assert.equal(q.valid, true);
    else {
      assert.equal(q.valid, false); assert.equal(q.halves.failures.length, 1); assert.equal(q.dual.failures.length, 1);
      const bad = q.dual.failures[0];
      assert.equal(bad.group, 0); assert.equal(bad.station, 57); assert.equal(bad.tube, 6);
      assert.equal(bad.intersections[0].type, 'proper');
    }
  }
});
