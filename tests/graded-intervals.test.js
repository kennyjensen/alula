import test from 'node:test';
import assert from 'node:assert/strict';
import { gradedIntervals } from '../src/geometry/graded-intervals.js';

test('graded intervals recover a known geometric series and exactly meet the far boundary', () => {
  const r = gradedIntervals({ extent: 15, firstSpacing: 1, intervals: 4 });
  r.x.forEach((v, i) => assert.ok(Math.abs(v - [0, 1, 3, 7, 15][i]) < 1e-13));
  assert.ok(Math.abs(r.growth - 2) < 1e-14);
  assert.equal(r.x.at(-1), 15);
});

test('graded end blocks preserve physical spacing, monotonic growth and units over the GUI counts', () => {
  for (const intervals of [16, 32, 64, 128]) {
    const r = gradedIntervals({ extent: 4, firstSpacing: .003, intervals });
    const scaled = gradedIntervals({ extent: 12, firstSpacing: .009, intervals });
    const widths = r.x.slice(1).map((x, i) => x - r.x[i]);
    assert.equal(r.x.length, intervals + 1); assert.equal(widths[0], .003);
    r.x.forEach((x, i) => assert.ok(Math.abs(3 * x - scaled.x[i]) < 2e-13));
    widths.slice(1).forEach((h, i) => {
      assert.ok(h >= widths[i]);
      assert.ok(Math.abs(h / widths[i] - r.growth) < 1e-11, 'the final interval must follow the same grading');
    });
  }
  const uniform = gradedIntervals({ extent: 8, firstSpacing: 4, intervals: 4 });
  assert.deepEqual(uniform.x, [0, 2, 4, 6, 8]);
  for (const options of [{ extent: 0 }, { firstSpacing: 0 }, { intervals: 1 }, { intervals: 2.5 }, { intervals: Infinity }])
    assert.throws(() => gradedIntervals({ extent: 4, firstSpacing: .01, intervals: 32, ...options }), /controls/);
});
