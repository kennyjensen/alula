import test from 'node:test';
import assert from 'node:assert/strict';
import { gaussUnitRule } from '../src/numerics/gauss-unit.js';

test('unit Gauss rules integrate every polynomial through their stated exactness degree', () => {
  for (const order of [3, 5]) {
    const { nodes, weights } = gaussUnitRule(order);
    assert.equal(nodes.length, order); assert.equal(weights.length, order);
    assert.ok(nodes.every((x, k) => x > 0 && x < 1 && (!k || x > nodes[k - 1])));
    assert.ok(weights.every(w => w > 0));
    for (let degree = 0; degree < 2 * order; degree++) {
      const value = nodes.reduce((sum, x, k) => sum + weights[k] * x ** degree, 0);
      assert.ok(Math.abs(value - 1 / (degree + 1)) < 4e-16, `Order ${order}, degree ${degree}: ${value}`);
    }
  }
  assert.throws(() => gaussUnitRule(4), /must be/);
});
