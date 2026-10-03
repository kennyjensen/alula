import test from 'node:test';
import assert from 'node:assert/strict';
import { allocateOrderedIntervalCounts } from '../src/numerics/ordered-interval-counts.js';

test('overlapping passage counts use compatible integer prefix sums and preserve fixed ends', () => {
  const input = { minimumCounts: [2, 2, 2, 2], fixedCounts: [{ block: 0, intervals: 2 }, { block: 3, intervals: 2 }],
    requirements: [{ from: 0, to: 2, minimum: 6 }, { from: 1, to: 3, minimum: 7 }] };
  const a = allocateOrderedIntervalCounts(input), before = structuredClone(input);
  assert.deepEqual(a.counts, [2, 4, 3, 2]); assert.equal(a.total, 11);
  let bruteMinimum = Infinity;
  for (let b = 2; b <= 8; b++) for (let c = 2; c <= 8; c++) if (2 + b >= 6 && b + c >= 7) bruteMinimum = Math.min(bruteMinimum, 4 + b + c);
  assert.equal(a.total, bruteMinimum); assert.deepEqual(input, before);
  assert.deepEqual(allocateOrderedIntervalCounts({ ...input, requirements: input.requirements.toReversed() }), a);
});

test('fixed equalities propagate backward and incompatible counts or budgets fail explicitly', () => {
  const a = allocateOrderedIntervalCounts({ minimumCounts: [2, 2, 2], fixedCounts: [{ block: 2, intervals: 2 }],
    requirements: [{ from: 0, to: 3, minimum: 10 }] });
  assert.deepEqual(a.counts, [2, 6, 2]);
  assert.throws(() => allocateOrderedIntervalCounts({ minimumCounts: [2, 2], fixedCounts: [{ block: 0, intervals: 2 }],
    requirements: [{ from: 0, to: 1, minimum: 3 }] }), /conflict/);
  assert.throws(() => allocateOrderedIntervalCounts({ minimumCounts: [4, 5], maximumTotal: 8 }), /budget/);
  assert.throws(() => allocateOrderedIntervalCounts({ minimumCounts: [2, 2], requirements: [{ from: 1, to: 0, minimum: 2 }] }), /Invalid/);
});
