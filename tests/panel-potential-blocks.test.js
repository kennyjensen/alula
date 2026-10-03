import test from 'node:test';
import assert from 'node:assert/strict';
import { createPanelPotentialBlocks } from '../src/geometry/panel-potential-blocks.js';

const asymmetric = () => ({
  bodies: [
    { leading: 4.8, trailing: { lower: 5, upper: 5.2 }, inlet: 0, outletIncrement: 6 },
    { leading: 4.1, trailing: { lower: 4.9, upper: 5.7 }, inlet: 0, outletIncrement: 6 },
  ],
  outer: [{ inlet: 0, outlet: 12 }, { inlet: 0, outlet: 12 }],
});
const close = (a, b, tolerance = 2e-13) => assert.ok(Math.abs(a - b) <= tolerance, a + ' != ' + b);
const event = (result, id) => result.events.find(e => e.id === id);

test('facing TE branches determine scalar block order, including propagated events', () => {
  const input = asymmetric(), before = structuredClone(input), blocks = createPanelPotentialBlocks(input);
  assert.deepEqual(blocks.order.map(e => e.id), ['1:LE', '0:LE', '1:TE', '0:TE']);
  assert.deepEqual(blocks.rank, { '1:LE': 1, '0:LE': 2, '1:TE': 3, '0:TE': 4 });
  assert.equal(blocks.inletRank, 0);
  assert.equal(blocks.outletRank, 5);
  assert.deepEqual(blocks.passages[1].map(e => e.potential), [4.1, 4.8, 4.9, 5.2]);
  // The main TE terminates on the lower body's upper wall; its event is
  // absent from the bottom passage rather than assigned a fictitious value.
  assert.deepEqual(blocks.passages[0].map(e => e.id), ['1:LE', '0:LE', '0:TE']);
  assert.deepEqual(blocks.passages[2].map(e => e.id), ['1:LE', '1:TE', '0:TE']);
  assert.deepEqual(input, before);
});

test('a shared upstream cut retains its gauge and a shared wake retains its increment', () => {
  const input = asymmetric(), blocks = createPanelPotentialBlocks(input);
  const upstream = event(blocks, '1:LE');
  assert.deepEqual(upstream.targets.lower, [4.1, 4.1]);
  assert.deepEqual(upstream.targets.upper, [4.1, 4.1]);
  assert.deepEqual(upstream.targets.outer, [4.1, 4.1]);
  assert.deepEqual(upstream.passagePotentials, [4.1, 4.1, 4.1]);
  const wake = event(blocks, '0:TE');
  close(wake.targets.upper[1] - input.bodies[1].trailing.upper,
    wake.targets.lower[1] - input.bodies[1].trailing.lower);
  close(wake.targets.lower[1], 5.2);
  close(wake.targets.upper[1], 6);
  close(wake.targets.outer[1], 6);
  assert.equal(wake.targets.outer[0], 5);
});

test('a wall terminates propagation and leaves all shielded targets unknown', () => {
  const blocks = createPanelPotentialBlocks({
    bodies: [
      { leading: 1, trailing: { lower: 2, upper: 2 }, inlet: -2, outletIncrement: 10 },
      { leading: .5, trailing: { lower: 3, upper: 3 }, inlet: -2, outletIncrement: 10 },
      { leading: 0, trailing: { lower: 4, upper: 4 }, inlet: -2, outletIncrement: 10 },
    ],
    outer: [{ inlet: -2, outlet: 20 }, { inlet: -2, outlet: 20 }],
  });
  for (const id of ['0:LE', '0:TE']) {
    const e = event(blocks, id);
    assert.equal(e.targets.lower[1], id === '0:LE' ? 1 : 2);
    assert.equal(e.targets.upper[1], null);
    assert.equal(e.targets.lower[2], null);
    assert.equal(e.targets.upper[2], null);
    assert.equal(e.targets.outer[1], null);
    assert.deepEqual(e.passagePotentials.slice(2), [null, null]);
  }
  assert.deepEqual(blocks.order.map(e => e.id), ['2:LE', '1:LE', '0:LE', '0:TE', '1:TE', '2:TE']);
});

test('known scalar boundary constraints include endpoints and omit shielded ranks', () => {
  const blocks = createPanelPotentialBlocks(asymmetric());
  assert.deepEqual(blocks.constraints.lower[0].map(c => [c.rank, c.potential]),
    [[0, 0], [1, 4.1], [2, 4.8], [4, 5], [5, 11]]);
  assert.deepEqual(blocks.constraints.upper[1].map(c => c.rank), [0, 1, 3, 4, 5]);
  for (const family of Object.values(blocks.constraints)) for (const row of family) {
    assert.equal(row[0].event, 'inlet');
    assert.equal(row.at(-1).event, 'outlet');
    for (let i = 1; i < row.length; i++) {
      assert.ok(row[i].rank > row[i - 1].rank);
      assert.ok(row[i].potential > row[i - 1].potential);
    }
  }
});

test('a global potential gauge and positive units change preserve topology and targets', () => {
  const input = asymmetric(), scale = 2.5, offset = -17, transform = x => scale * x + offset;
  const converted = {
    bodies: input.bodies.map(b => ({ leading: transform(b.leading),
      trailing: { lower: transform(b.trailing.lower), upper: transform(b.trailing.upper) },
      inlet: transform(b.inlet), outletIncrement: scale * b.outletIncrement })),
    outer: input.outer.map(b => ({ inlet: transform(b.inlet), outlet: transform(b.outlet) })),
  };
  const a = createPanelPotentialBlocks(input), b = createPanelPotentialBlocks(converted);
  assert.deepEqual(b.rank, a.rank);
  for (const ae of a.events) {
    const be = event(b, ae.id);
    for (const side of ['lower', 'upper', 'outer']) for (let j = 0; j < ae.targets[side].length; j++) {
      const target = ae.targets[side][j];
      if (target === null) assert.equal(be.targets[side][j], null);
      else close(be.targets[side][j], transform(target));
    }
  }
});

test('coincident block boundaries, unavailable paths and invalid source order fail explicitly', () => {
  const coincident = asymmetric(); coincident.bodies[0].leading = 4.9;
  assert.throws(() => createPanelPotentialBlocks(coincident), /Coincident/);
  const outsideOuter = asymmetric(); outsideOuter.outer[1].outlet = 5.8;
  assert.throws(() => createPanelPotentialBlocks(outsideOuter), /outside outer/);
  const outsideWake = asymmetric(); outsideWake.bodies[1].outletIncrement = .1;
  assert.throws(() => createPanelPotentialBlocks(outsideWake), /outside body/);
  const corner = asymmetric(); corner.outer[0].outlet = 4.8;
  assert.throws(() => createPanelPotentialBlocks(corner), /Coincident/);
  const reverse = asymmetric(); reverse.bodies[0].trailing.upper = 4.7;
  assert.throws(() => createPanelPotentialBlocks(reverse), /must increase/);
  const nonfinite = asymmetric(); nonfinite.bodies[1].trailing.lower = NaN;
  assert.throws(() => createPanelPotentialBlocks(nonfinite), /Invalid/);
  assert.throws(() => createPanelPotentialBlocks(), /require bodies/);
});

// Cyclic common-index graphs are tested independently in
// potential-block-order.test.js. This constructor delegates their rejection
// to that utility; its monotone transfers do not invent a cycle fixture.

