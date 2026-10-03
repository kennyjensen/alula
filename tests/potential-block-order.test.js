import test from 'node:test';
import assert from 'node:assert/strict';
import { orderPotentialBlockEvents } from '../src/geometry/potential-block-order.js';

const event = (body, edge, potential) => ({ body, edge, potential });
// Two lifting elements: a TE has different upper/lower potential values.
// Means would put body 0's TE first (5.1 < 5.3), while the physical gap
// requires body 1's TE first (4.9 < 5.2).
const passages = [
  [event(0, 'LE', 4.8), event(0, 'TE', 5)],
  [event(0, 'LE', 4.8), event(0, 'TE', 5.2), event(1, 'LE', 4.1), event(1, 'TE', 4.9)],
  [event(1, 'LE', 4.1), event(1, 'TE', 5.7)],
];

test('block order respects the facing branches instead of mean trailing-edge potentials', () => {
  const original = structuredClone(passages), result = orderPotentialBlockEvents(passages);
  assert.deepEqual(result.order.map(e => e.id), ['1:LE', '0:LE', '1:TE', '0:TE']);
  for (const row of result.passages) for (let i = 1; i < row.length; i++)
    assert.ok(result.rank[row[i - 1].id] < result.rank[row[i].id]);
  assert.deepEqual(passages, original);
});

test('independent passage potential gauges, units and input enumeration do not change topology', () => {
  const expected = orderPotentialBlockEvents(passages).order;
  const changed = passages.map((row, g) => row.map(e => ({ ...e, potential: 3 * e.potential + 19 * g - 7 })).reverse()).reverse();
  assert.deepEqual(orderPotentialBlockEvents(changed).order, expected);
});

test('conflicting or degenerate block constraints are rejected explicitly', () => {
  assert.throws(() => orderPotentialBlockEvents([
    [event(0, 'TE', 1), event(1, 'TE', 2)],
    [event(1, 'TE', 3), event(0, 'TE', 4)],
  ]), /Conflicting/);
  assert.throws(() => orderPotentialBlockEvents([[event(0, 'LE', 1), event(1, 'LE', 1)]]), /Coincident/);
  assert.throws(() => orderPotentialBlockEvents([[event(0, 'LE', 1), event(0, 'TE', 0)]]), /must increase/);
  assert.throws(() => orderPotentialBlockEvents([[event(0, 'LE', 0), event(0, 'LE', 1)]]), /Duplicate/);
  assert.throws(() => orderPotentialBlockEvents([[event(0, 'LE', 0), event(0, 'TE', NaN)]]), /Invalid/);
});
