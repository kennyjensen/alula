import test from 'node:test';
import assert from 'node:assert/strict';
import { redistributeWithSharedCoordinate } from '../scripts/validation/shared-smove-coordinate.js';
import { redistributeStreamtubeTangentially } from '../src/geometry/streamtube-tangential-redistribution.js';

test('shared coordinate reproduces production when the reference is unchanged and respects coordinate scaling', () => {
  const nx = 9, nt = 5;
  const nodes = Array.from({ length: nx + 1 }, (_, i) => Array.from({ length: nt + 1 }, (_, j) => {
    const u = i / nx, v = j / nt;
    return { x: u + .07 * Math.sin(Math.PI * u) * (1 + .2 * v), y: v + .02 * Math.sin(2 * Math.PI * u) * v * (1 - v) };
  }));
  const saved = structuredClone(nodes);
  for (const [referenceBank, fixedBanks] of [[0, [true, false]], [nt, [false, true]]]) {
    const options = { referenceBank, fixedBanks }, reference = nodes.map(row => row[referenceBank]);
    const production = redistributeStreamtubeTangentially(nodes, options);
    const control = redistributeWithSharedCoordinate(nodes, reference, options);
    assert.deepEqual(control.nodes, production.nodes);
    assert.deepEqual(control.solution, production.solution);
    // Scaling the coordinate labels alone cannot change the physical move.
    const rescaled = redistributeWithSharedCoordinate(nodes, reference.map(p => ({ x: 16 * p.x, y: 16 * p.y })), options);
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) {
      assert.ok(Math.hypot(rescaled.nodes[i][j].x - control.nodes[i][j].x,
        rescaled.nodes[i][j].y - control.nodes[i][j].y) < 2e-14);
      if (i === 0 || i === nx || fixedBanks[0] && j === 0 || fixedBanks[1] && j === nt)
        assert.deepEqual(control.nodes[i][j], nodes[i][j]);
    }
  }
  assert.deepEqual(nodes, saved);
  assert.throws(() => redistributeWithSharedCoordinate(nodes, [{ x: 0, y: 0 }]), /complete/);
});
