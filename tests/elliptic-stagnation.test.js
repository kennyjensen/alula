import test from 'node:test';
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { stagnationPoint, stagnationCenterOracle, runStagnationCase, stagnationRefinement } from '../scripts/validation/elliptic-stagnation.js';

test('square-root reference inverts analytic stagnation potential and streamfunction, including the exact corner', () => {
  for (const phi of [-1, -.2, 0, .2, 1]) for (const psi of [0, 1e-14, .125, 1]) {
    const { x, y } = stagnationPoint(phi, psi);
    assert.ok(Math.abs(x * x - y * y - phi) < 8e-16);
    assert.ok(Math.abs(2 * x * y - psi) < 4e-16);
    assert.ok(x >= 0 && y >= 0);
  }
  assert.deepEqual(stagnationPoint(0, 0), { x: 0, y: 0 });
});

test('Giles residual at stagnation agrees with an independently derived singular-stencil coefficient', () => {
  for (const nt of [4, 8, 16]) for (const ratio of [1, 2]) {
    const nx = 2 * nt / ratio, h = 1 / nt;
    const nodes = Array.from({ length: nx + 1 }, (_, i) =>
      Array.from({ length: nt + 1 }, (_, j) => stagnationPoint(2 * i / nx - 1, j * h)));
    const system = createEllipticStreamtubeGrid({ nodes, massFlows: Array(nt).fill(h), discretization: 'giles-1985' });
    const row = system.residuals(nodes).rows.find(r => r.i === nx / 2 && r.j === 1);
    const oracle = stagnationCenterOracle(h, ratio);
    const expected = oracle.residualX / ((oracle.alpha + oracle.gamma) * system.lengthScale);
    assert.ok(Math.abs(row.x / expected - 1) < 2e-13);
    assert.ok(Math.abs(row.y / expected - 1) < 2e-13);
    assert.ok(expected < 0);
  }
});

test('documented SLOR refines the exact stagnation flow while local-tube error need not vanish', () => {
  for (const settings of [{ offset: .25, stretch: 0 }, { offset: 0, stretch: 0 }, { offset: 0, stretch: 3 }]) {
    const cases = [8, 16].map(nt => runStagnationCase({ nt, ...settings }));
    for (const result of cases) {
      assert.equal(result.converged, true, result.reason);
      assert.equal(result.invalidCells, 0);
      assert.equal(result.changedBoundaries, 0);
      assert.ok(result.residual < 1e-9);
    }
    const order = stagnationRefinement(...cases);
    assert.ok(order.absolute > .8 && order.awayFromStagnation > 1.7, JSON.stringify(order));
    if (settings.offset) assert.ok(order.absolute > 1.8);
    else {
      assert.ok(cases.every(c => c.localPeak.i === c.nx / 2 && c.localPeak.j === 1));
      // The denominator shrinks with refinement. This remains a diagnostic,
      // not a new airfoil acceptance limit or permission to ignore error.
      assert.ok(cases.every(c => c.errors.localIntervals > .05));
      assert.ok(Math.abs(order.localIntervals) < .15);
    }
  }
});
