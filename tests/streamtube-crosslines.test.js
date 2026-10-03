import test from 'node:test';
import assert from 'node:assert/strict';
import { createBodyPotentialMap, createPotentialCrosslines } from '../src/geometry/streamtube-crosslines.js';

const near = (a, b, tolerance = 2e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} differs from ${b}`);

test('potential block maps preserve stagnation, separate TE potentials, wake circulation and invert monotonically', () => {
  const anchors = { start: -4, end: 5, leading: -.2, upperTrailing: 1.3, lowerTrailing: .7, inlet: -3.8, outletIncrement: 3.7 };
  const map = createBodyPotentialMap(anchors);
  for (const side of ['upper', 'lower']) {
    near(map.value(side, -4), -3.8); near(map.value(side, -.2), -.2);
    near(map.value(side, 1), anchors[side + 'Trailing']);
    near(map.value(side, 5), anchors[side + 'Trailing'] + 3.7);
    let previous = -Infinity;
    for (let i = 0; i <= 900; i++) {
      const x = -4 + i / 100, phi = map.value(side, x);
      assert.ok(phi > previous); previous = phi; near(map.inverse(side, phi), x);
      if (i > 0 && i < 900) {
        const h = 1e-6, derivative = (map.value(side, x + h) - map.value(side, x - h)) / (2 * h);
        assert.ok(derivative >= map.minimumDerivative - 2e-9);
      }
    }
  }
  for (const x of [-4, -1, -.2]) near(map.value('upper', x), map.value('lower', x));
  for (const x of [1, 1.1, 3, 5]) near(map.value('upper', x) - map.value('lower', x), .6);
  const transformed = createBodyPotentialMap(Object.fromEntries(Object.entries(anchors).map(([key, value]) => [key, 3 * value + (key === 'outletIncrement' ? 0 : 7)])));
  for (const side of ['upper', 'lower']) for (const x of [-4, -.7, -.2, .5, 1, 2, 5]) near(transformed.value(side, 3 * x + 7), 3 * map.value(side, x) + 7);
  assert.throws(() => createBodyPotentialMap({ ...anchors, inlet: -.21 }), /not monotone/);
  assert.throws(() => createBodyPotentialMap({ ...anchors, lowerTrailing: -.1 }), /not monotone/);
  assert.throws(() => map.inverse('upper', -10), /outside/);
});

test('common spacing integrates an independent uniform-body/geometric-tail metric exactly', () => {
  const row = [0, .25, .5, .75, 1], r = createPotentialCrosslines({ profiles: [row], start: -2, end: 3, growth: .5 });
  // Integral of 1/(.25+.5*d) over either length-two tail is 2*log(5).
  near(r.metricLength, 4 + 4 * Math.log(5));
  assert.deepEqual(r.anchors, [-2, 0, 1, 3]);
  assert.deepEqual(r.blocks.map(b => b.intervals), [4, 4, 4]);
  assert.deepEqual(r.x.slice(4, 9), row);
  // Each tail interval has the same logarithmic metric measure.
  for (let i = 0; i < 4; i++) {
    const a = r.x[i], b = r.x[i + 1];
    near(2 * Math.log((.25 - .5 * a) / (.25 - .5 * b)), .5 * Math.log(5));
  }
  const custom = createPotentialCrosslines({ profiles: [row], start: -2, end: 3, growth: .5, inletIntervals: 7, outletIntervals: 9 });
  assert.deepEqual(custom.blocks.map(b => b.intervals), [7, 4, 9]);
  assert.deepEqual(custom.x.slice(7, 12), row);
  for (let i = 0; i < 7; i++) near(2 * Math.log((.25 - .5 * custom.x[i]) / (.25 - .5 * custom.x[i + 1])), 2 * Math.log(5) / 7);
  assert.throws(() => createPotentialCrosslines({ profiles: [row], start: -2, end: 3, inletIntervals: 1 }), /Invalid/);
});

test('facing-surface spacing merges demands without duplicate-node slivers and is invariant to units and order', () => {
  const profiles = [[0, .1, .3, .7, 1], [0, .1 + 1e-12, .3 + 1e-12, .7 + 1e-12, 1], [.2, .35, .6, .8]];
  const r = createPotentialCrosslines({ profiles, start: -1, end: 2 });
  assert.deepEqual(r.anchors, [-1, 0, .2, .8, 1, 2]);
  assert.ok(Math.min(...r.x.slice(1).map((v, i) => v - r.x[i])) > .01);
  const reversed = createPotentialCrosslines({ profiles: profiles.toReversed(), start: -1, end: 2 });
  const scaled = createPotentialCrosslines({ profiles: profiles.map(row => row.map(x => 2 * x - 3)), start: -5, end: 1 });
  assert.equal(r.x.length, reversed.x.length); assert.equal(r.x.length, scaled.x.length);
  r.x.forEach((x, i) => { near(reversed.x[i], x); near(scaled.x[i], 2 * x - 3); });
  assert.throws(() => createPotentialCrosslines({ profiles: [[0, .1, 1], [1e-12, .2, 1]], start: -1, end: 2 }), /coincident/);
  assert.throws(() => createPotentialCrosslines({ profiles, start: -1, end: 2, maxIntervals: 4 }), /budget/);
  assert.throws(() => createPotentialCrosslines({ profiles: [[0, .2, .1]], start: -1, end: 2 }), /Invalid/);
});
test('explicit block counts retain exact anchors and enforce the total station budget', () => {
  const input = { profiles: [[-1, 0, 1]], start: -2, end: 2, blockIntervals: [6, 11, 9] };
  const r = createPotentialCrosslines(input);
  assert.deepEqual(r.blocks.map(b => b.intervals), [6, 11, 9]);
  assert.equal(r.x[6], -1); assert.equal(r.x[17], 1); assert.equal(r.x.length, 27);
  assert.throws(() => createPotentialCrosslines({ ...input, blockIntervals: [6, 11] }), /per potential block/);
  assert.throws(() => createPotentialCrosslines({ ...input, maxIntervals: 25 }), /budget/);
});
