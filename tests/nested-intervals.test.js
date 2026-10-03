import test from 'node:test';
import assert from 'node:assert/strict';
import { nestedIntervalMap } from '../src/geometry/nested-intervals.js';

test('local station refinement retains parent stations and uses each requested interval count', () => {
  const counts = [1, 4, 2, 3, 1], r = nestedIntervalMap(counts.length, { subdivisions: counts });
  assert.deepEqual(r.retained, [0, 1, 5, 7, 10, 11]);
  assert.equal(r.uniformFactor, null);
  for (let j = 0; j < counts.length; j++) {
    assert.equal(r.coordinates[r.retained[j]], j);
    for (let k = 0; k <= counts[j]; k++) {
      assert.equal(r.coordinate(r.retained[j] + k, j), k / counts[j]);
      assert.ok(Math.abs(r.coordinates[r.retained[j] + k] - (j + k / counts[j])) < 1e-15);
    }
  }
  assert.equal(r.coordinates.at(-1), counts.length);
  assert.deepEqual(counts, [1, 4, 2, 3, 1]);
});

test('uniform station refinement retains previous global and branch-local arithmetic exactly', () => {
  for (const factor of [1, 2, 3, 4]) {
    const a = nestedIntervalMap(11, { factor }), b = nestedIntervalMap(11, { subdivisions: Array(11).fill(factor) });
    assert.deepEqual(a.coordinates, b.coordinates);
    for (let origin = 0; origin < 11; origin++) for (let i = factor * origin; i <= factor * 11; i++) {
      assert.equal(a.coordinate(i), i / factor);
      assert.equal(a.coordinate(i, origin), (i - factor * origin) / factor);
    }
  }
  assert.equal(nestedIntervalMap(3).uniformFactor, 2);
});

test('ambiguous or malformed station refinement is rejected', () => {
  for (const controls of [{factor:0},{factor:null},{factor:5},{factor:1.5},{subdivisions:null},
    {subdivisions:[1]},{subdivisions:[1,0,2]},{factor:2,subdivisions:[2,2,2]}])
    assert.throws(() => nestedIntervalMap(3, controls), /refinement/);
});
