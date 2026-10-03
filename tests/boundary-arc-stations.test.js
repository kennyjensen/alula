import test from 'node:test';
import assert from 'node:assert/strict';
import { createBoundaryArcStationMap } from '../src/geometry/boundary-arc-stations.js';
const points = intervals => { const p = [{ x: 0, y: 0 }]; intervals.forEach(h => p.push({ x: p.at(-1).x + h, y: 0 })); return p; };
const ratio = positions => {
  const h = positions.slice(1).map((s, i) => s - positions[i]);
  return h.slice(1).reduce((r, value, i) => Math.max(r, value / h[i], h[i] / value), 1);
};

test('convex boundary arc transfer retains ordered stations, endpoint values and the growth bound', () => {
  const a = [1, 1.4, 1.1, .9, .7], b = [.8, .6, .7, .95, 1.2];
  const lower = points(a), upper = points(b).map(p => ({ x: 2 * p.x, y: 3 }));
  const map = createBoundaryArcStationMap({ lower, upper });
  for (const fraction of [0, .001, .25, .5, .9, .999, 1]) {
    const s = map.at(fraction); assert.equal(s[0], 0); assert.equal(s.at(-1), 1);
    assert.ok(ratio(s) <= 1.4 + 1e-13);
    s.forEach((value, i) => { if(i)assert.ok(value > s[i - 1]); });
  }
  for (const [path, f] of [[lower, 0], [upper, 1]]) path.forEach((p, i) => assert.ok(Math.abs(map.at(f)[i] - p.x / path.at(-1).x) < 1e-14));
  assert.ok(Math.abs(map.maximumAdjacentRatio - 1.4) < 1e-14);
});

test('boundary arc control is independent of rigid transforms and scale; curved boundaries use polygon length explicitly', () => {
  const lower = [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: 2 }], upper = points([2, 2]);
  const a = createBoundaryArcStationMap({ lower, upper });
  assert.equal(a.lowerLength, 3); assert.equal(a.upperLength, 4);
  assert.equal(a.at(.5)[1], .5 * (1 / 3 + .5));
  const transform = path => path.map(p => ({ x: 7 - 3 * p.y, y: -2 + 3 * p.x }));
  const b = createBoundaryArcStationMap({ lower: transform(lower), upper: transform(upper) });
  assert.deepEqual(a.at(.3), b.at(.3));
  assert.throws(() => createBoundaryArcStationMap({ lower: [lower[0], lower[0]], upper }), /positive/);
  assert.throws(() => createBoundaryArcStationMap({ lower: lower.slice(1), upper }), /counts/);
  assert.throws(() => a.at(1.1), /fraction/);
});
