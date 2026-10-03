import test from 'node:test';
import assert from 'node:assert/strict';
import { createPanelPotentialBlocks } from '../src/geometry/panel-potential-blocks.js';
import { createPassageSpacingComponents } from '../src/geometry/passage-spacing-components.js';
import { distributeFacingDensity } from '../src/geometry/facing-density-stations.js';

const source = () => ({ bodies: [
  { leading: 4.8, trailing: { lower: 5, upper: 5.2 }, inlet: 0, outletIncrement: 6 },
  { leading: 4.1, trailing: { lower: 4.9, upper: 5.7 }, inlet: 0, outletIncrement: 6 },
], outer: [{ inlet: 0, outlet: 12 }, { inlet: 0, outlet: 12 }] });

test('main/flap passages form seven spacing components with cuts shared and walls shielding', () => {
  const blocks = createPanelPotentialBlocks(source()), before = structuredClone(blocks);
  const c = createPassageSpacingComponents(blocks);
  assert.deepEqual(c.map(({ start, end, passages }) => [start, end, passages]), [
    [0, 1, [0, 1, 2]], [1, 2, [0, 1]], [1, 3, [2]], [2, 3, [1]], [2, 4, [0]], [3, 4, [1, 2]], [4, 5, [0, 1, 2]],
  ]);
  assert.deepEqual(c[1].boundaries, [{ kind: 'outer', side: 'lower' }, { kind: 'cut', body: 0, end: 'upstream' }, { kind: 'wall', body: 1, side: 'lower' }]);
  assert.deepEqual(c[2].boundaries, [{ kind: 'wall', body: 1, side: 'upper' }, { kind: 'outer', side: 'upper' }]);
  assert.deepEqual(c[3].boundaries, [{ kind: 'wall', body: 0, side: 'upper' }, { kind: 'wall', body: 1, side: 'lower' }]);
  assert.deepEqual(c[5].boundaries, [{ kind: 'wall', body: 0, side: 'upper' }, { kind: 'cut', body: 1, end: 'wake' }, { kind: 'outer', side: 'upper' }]);
  assert.deepEqual(blocks, before);
  // Each actual passage interval is represented once; shielded global
  // event ranks do not split that passage's interval.
  const got = c.flatMap(q => q.passages.map(g => `${g}:${q.start}:${q.end}`));
  const expected = blocks.passages.flatMap((events, g) => {
    const ranks = [0, ...events.map(e => blocks.rank[e.id]), 5];
    return ranks.slice(1).map((r, k) => `${g}:${ranks[k]}:${r}`);
  });
  assert.deepEqual(got.sort(), expected.sort()); assert.equal(new Set(got).size, got.length);
});

test('an isolated foil and a nested three-element stack retain opposite-side independence', () => {
  for (const n of [1, 3]) {
    const bodies = Array.from({ length: n }, (_, i) => ({ leading: 1 - .2 * i,
      trailing: { lower: 2 + .2 * i, upper: 2 + .2 * i }, inlet: -2, outletIncrement: 10 }));
    const b = createPanelPotentialBlocks({ bodies, outer: [{ inlet: -2, outlet: 20 }, { inlet: -2, outlet: 20 }] });
    const c = createPassageSpacingComponents(b);
    for (let body = 0; body < n; body++) for (const q of c) {
      const sides = q.boundaries.filter(p => p.kind === 'wall' && p.body === body).map(p => p.side);
      assert.ok(sides.length <= 1, 'a solid must not connect its two spacing requests');
    }
    for (const q of c) for (const boundary of q.boundaries.filter(p => p.kind === 'cut')) {
      assert.ok(q.end <= b.rank[`${boundary.body}:LE`] || q.start >= b.rank[`${boundary.body}:TE`]);
    }
    if (n === 1) assert.deepEqual(c.map(q => q.passages), [[0, 1], [0], [1], [0, 1]]);
  }
});

test('potential gauge and units cannot change the components', () => {
  const a = source(), transform = x => 3.7 * x - 19;
  const b = { bodies: a.bodies.map(p => ({ leading: transform(p.leading), inlet: transform(p.inlet),
    trailing: { lower: transform(p.trailing.lower), upper: transform(p.trailing.upper) }, outletIncrement: p.outletIncrement * 3.7 })),
    outer: a.outer.map(p => ({ inlet: transform(p.inlet), outlet: transform(p.outlet) })) };
  assert.deepEqual(createPassageSpacingComponents(createPanelPotentialBlocks(a)), createPassageSpacingComponents(createPanelPotentialBlocks(b)));
});

test('missing cut-bank events and malformed endpoint orders are rejected', () => {
  const a = createPanelPotentialBlocks(source()); a.passages[0].shift();
  assert.throws(() => createPassageSpacingComponents(a), /shared cut|unsplit/);
  const b = createPanelPotentialBlocks(source()); b.rank['0:TE'] = b.rank['0:LE'];
  assert.throws(() => createPassageSpacingComponents(b), /ordered body/);
  assert.throws(() => createPassageSpacingComponents({}), /skeleton/);
});

test('changing flap density affects its facing passage but cannot refine the shielded main upper side', () => {
  const components = createPassageSpacingComponents(createPanelPotentialBlocks(source()));
  const uniform = Array.from({ length: 9 }, (_, i) => i / 8);
  const requests = { '0:upper': uniform, '0:lower': uniform, '1:upper': uniform, '1:lower': uniform };
  const distribute = demand => components.map(c => {
    const profiles = c.boundaries.filter(b => b.kind === 'wall').map(b => demand[`${b.body}:${b.side}`]);
    return profiles.length ? distributeFacingDensity({ profiles, intervals: 32 }).progress : null;
  });
  const before = distribute(requests);
  const after = distribute({ ...requests, '0:upper': [0, .01, .03, .08, .2, .4, .65, .85, 1] });
  const shielded = components.findIndex(c => c.start === 1 && c.end === 3);
  const facing = components.findIndex(c => c.start === 2 && c.end === 3);
  assert.deepEqual(after[shielded], before[shielded]);
  assert.ok(Math.max(...after[facing].map((u, i) => Math.abs(u - before[facing][i]))) > .1);
});
