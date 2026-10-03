import test from 'node:test';
import assert from 'node:assert/strict';
import { createInitialStreamtubeTopology } from '../src/geometry/streamtube-topology.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const input = { elements: [{ points: naca4('2412', 80) },
  { points: transform(naca4('0012', 80), { chord: .3, x: .94, y: -.08, angle: -15 }) }], alpha: 4 };
const fractions = weights => {
  const total = weights.reduce((a, b) => a + b, 0);
  return weights.map(w => w / total);
};
const cumulative = weights => fractions(weights).reduce((a, w) => [...a, a.at(-1) + w], [0]);

test('default normal refinement subdivides the same mass coordinates in every fluid region', () => {
  const grids = [7, 14, 28].map(tubes => createInitialStreamtubeTopology(input, { tubes }));
  // Independent refinement identity: each coarse tube's mass is the sum
  // of its children, and all previous streamline levels are retained.
  grids.slice(1).forEach((grid, level) => {
    grid.weights.forEach((weights, g) => {
      const coarse = fractions(grids[level].weights[g]), fine = fractions(weights);
      const coarseLevels = cumulative(grids[level].weights[g]), fineLevels = cumulative(weights);
      coarse.forEach((mass, i) => {
        assert.ok(Math.abs(mass - fine[2 * i] - fine[2 * i + 1]) < 1e-14);
        assert.ok(Math.abs(coarseLevels[i] - fineLevels[2 * i]) < 1e-14);
      });
      assert.ok(fine.every(w => w > 0));
    });
    assert.deepEqual(grid.bodies, grids[0].bodies);
    assert.deepEqual(grid.outerLower, grids[0].outerLower);
  });
});

test('11 tubes refine the original seven-tube profile without exponential wall crowding; explicit growth remains supported', () => {
  const baseline = createInitialStreamtubeTopology(input, { tubes: 7 });
  assert.deepEqual(baseline, createInitialStreamtubeTopology(input, { tubes: 7, tubeGrowth: 3 }));
  const refined = createInitialStreamtubeTopology(input, { tubes: 11 });
  const old = createInitialStreamtubeTopology(input, { tubes: 11, tubeGrowth: 3 });
  for (const g of [0, 2]) {
    const first = Math.min(...fractions(baseline.weights[g]));
    const next = Math.min(...fractions(refined.weights[g]));
    assert.ok(first / next > 1.5 && first / next < 2.5);
    assert.ok(first / Math.min(...fractions(old.weights[g])) > 80);
  }
  const uniform = createInitialStreamtubeTopology(input, { tubes: 11, tubeGrowth: 1 });
  for (const g of [0, 2]) assert.deepEqual(uniform.weights[g], Array(11).fill(1));
});

test('inlet and wake station counts are independent, graded, and preserve the body and far boundaries', () => {
  const a = createInitialStreamtubeTopology(input), b = createInitialStreamtubeTopology(input, { inletIntervals: 16, outletIntervals: 64 });
  for (const [grid, expected] of [[a, [8, 8]], [b, [16, 64]]]) {
    const leading = Math.min(...grid.bodies.map(b => b.leadingIndex));
    const trailing = Math.max(...grid.bodies.map(b => b.trailingIndex));
    assert.equal(leading, expected[0]); assert.equal(grid.outerLower.length - 1 - trailing, expected[1]);
    const x = grid.outerLower.map(p => p.x);
    const inlet = x.slice(0, leading + 1).reverse().map(v => x[leading] - v);
    const outlet = x.slice(trailing).map(v => v - x[trailing]);
    for (const row of [inlet, outlet]) {
      const widths = row.slice(1).map((v, i) => v - row[i]);
      widths.slice(1).forEach((v, i) => assert.ok(v >= widths[i] - 2e-12));
    }
  }
  assert.deepEqual(a.outerLower[0], b.outerLower[0]); assert.deepEqual(a.outerLower.at(-1), b.outerLower.at(-1));
  a.bodies.forEach((body, i) => {
    assert.deepEqual(body.points, b.bodies[i].points);
    assert.deepEqual(body.surfaceFractions, b.bodies[i].surfaceFractions);
    assert.equal(body.trailingIndex - body.leadingIndex, b.bodies[i].trailingIndex - b.bodies[i].leadingIndex);
  });
  for (const controls of [{ inletIntervals: 0 }, { outletIntervals: NaN }, { outletIntervals: 257 }])
    assert.throws(() => createInitialStreamtubeTopology(input, controls), /controls/);
});

test('64 and 128 surface intervals refine each element while fixed inlet/wake counts remain unchanged', () => {
  for (const surfaceIntervals of [64, 128]) {
    const grid = createInitialStreamtubeTopology(input, { surfaceIntervals, inletIntervals: 32, outletIntervals: 32 });
    assert.equal(grid.gridSpacing.surfaceIntervals, surfaceIntervals);
    assert.equal(Math.min(...grid.bodies.map(b => b.leadingIndex)), 32);
    assert.equal(grid.outerLower.length - 1 - Math.max(...grid.bodies.map(b => b.trailingIndex)), 32);
    grid.bodies.forEach((body, i) => {
      assert.ok(body.trailingIndex - body.leadingIndex >= surfaceIntervals);
      assert.deepEqual(grid.gridSpacing.surfaceIntervalsByElement[i], { element: body.element, intervals: body.trailingIndex - body.leadingIndex });
      assert.equal(body.surfaceFractions.length, body.trailingIndex - body.leadingIndex + 1);
      assert.ok(body.surfaceFractions.every((f, i, row) => i === 0 || f > row[i - 1]));
    });
  }
  assert.throws(() => createInitialStreamtubeTopology(input, { surfaceIntervals: 129 }), /controls/);
});

test('automatic inlet and wake counts track half the surface intervals', () => {
  for (const surfaceIntervals of [8, 16, 32, 64, 128]) {
    const grid = createInitialStreamtubeTopology(input, { surfaceIntervals });
    assert.equal(Math.min(...grid.bodies.map(b => b.leadingIndex)), surfaceIntervals / 2);
    assert.equal(grid.outerLower.length - 1 - Math.max(...grid.bodies.map(b => b.trailingIndex)), surfaceIntervals / 2);
  }
});
