import test from 'node:test';
import assert from 'node:assert/strict';
import { createPhysicalCutStations } from '../src/geometry/tests/physical-cut-stations.js';

const close = (a, b) => assert.ok(Math.abs(a - b) < 2e-12 * Math.max(1, Math.abs(b)), `${a} != ${b}`);

test('physical station blocks preserve hits/counts and match actual intervals across different counts', () => {
  const r = createPhysicalCutStations({ anchors: [{ index: 4, x: 0 }, { index: 7, x: 7 }, { index: 9, x: 19 }], edgeSpacing: 1 });
  [0, 1, 3, 7, 11, 19].forEach((v, i) => close(r.x[i], v));
  assert.equal(r.firstIndex, 4); assert.equal(r.x[3], 7); assert.equal(r.x[5], 19);
  assert.deepEqual(r.blocks.map(b => b.intervals), [3, 2]);
  close(r.blocks[0].lastSpacing, r.blocks[1].firstSpacing);
});

test('inlet direction reverses exactly and finite distance spacing avoids the stagnation square-root singularity', () => {
  const r = createPhysicalCutStations({ anchors: [{ index: 0, x: -19 }, { index: 2, x: -7 }, { index: 5, x: 0 }], edgeSpacing: 1, bodyEnd: 'last' });
  [-19, -11, -7, -3, -1, 0].forEach((v, i) => close(r.x[i], v));
  const potentials = r.x.map(x => 12 - .5 * x * x);
  // Analytic stagnation field u=-x. Potential is quadratic, although the
  // prescribed near-body physical step is finite and independent of phi.
  close(potentials.at(-1) - potentials.at(-2), .5);
  close(r.x.at(-1) - r.x.at(-2), 1);
});

test('a first step larger than the mean is retained with a decreasing ratio', () => {
  const r = createPhysicalCutStations({ anchors: [{ index: 0, x: 0 }, { index: 3, x: 7 }], edgeSpacing: 4 });
  [0, 4, 6, 7].forEach((v, i) => close(r.x[i], v));
  close(r.blocks[0].growth, .5);
});

test('uniform, translated and very stretched schedules retain strictly positive intervals', () => {
  for (const [origin, length, n, first] of [[2, 8, 4, 2], [0, 100, 128, 1e-7], [-3, 20, 256, .002]]) {
    const r = createPhysicalCutStations({ anchors: [{ index: 3, x: origin }, { index: n + 3, x: origin + length }], edgeSpacing: first });
    assert.equal(r.x[0], origin); assert.equal(r.x.at(-1), origin + length);
    close(r.x[1] - r.x[0], first);
    assert.ok(r.x.every((v, i) => !i || v > r.x[i - 1]));
    const steps = r.x.slice(1).map((v, i) => v - r.x[i]);
    steps.slice(1).forEach((v, i) => close(v / steps[i], r.blocks[0].growth));
  }
});

test('incompatible constraints and unresolved coordinates fail without clipping the request', () => {
  const anchors = [{ index: 0, x: 0 }, { index: 2, x: 1 }];
  assert.throws(() => createPhysicalCutStations({ anchors, edgeSpacing: 1 }), /smaller/);
  assert.throws(() => createPhysicalCutStations({ anchors: [{ index: 0, x: 0 }, { index: 1, x: 1 }], edgeSpacing: .5 }), /One-interval/);
  assert.throws(() => createPhysicalCutStations({ anchors: [{ index: 0, x: 1e10 }, { index: 10, x: 1e10 + 1 }], edgeSpacing: 1e-10 }), /collapse/);
  assert.throws(() => createPhysicalCutStations({ anchors, edgeSpacing: 0 }), /positive/);
  assert.throws(() => createPhysicalCutStations({ anchors: anchors.toReversed(), edgeSpacing: .1 }), /ordered/);
  assert.throws(() => createPhysicalCutStations({ anchors: [{ index: 0, x: 2 ** 40 }, { index: 2, x: 2 ** 40 + 1 }],
    edgeSpacing: 3 * 2 ** -14 }), /accuracy is unresolved/);
  assert.throws(() => createPhysicalCutStations({ anchors: [{ index: 2 ** 53, x: 0 }, { index: 2 ** 53 + 4, x: 7 }],
    edgeSpacing: 1 }), /ordered/);
  assert.throws(() => createPhysicalCutStations({ anchors: [{ index: 0, x: 0 }, { index: 3, x: 7 }, { index: 5, x: 9 }],
    edgeSpacing: 1 }), /smaller/);
});

test('feasible one-interval schedules preserve the input in either direction', () => {
  const anchors = [{ index: 3, x: 2 }, { index: 4, x: 5 }], before = structuredClone(anchors);
  for (const bodyEnd of ['first', 'last']) {
    const r = createPhysicalCutStations({ anchors, edgeSpacing: 3, bodyEnd });
    assert.deepEqual(r.x, [2, 5]); assert.equal(r.blocks[0].achievedFirstSpacing, 3);
  }
  assert.deepEqual(anchors, before);
});

test('a tiny tail and large dynamic range do not lose the fitted remaining length', () => {
  const first = 1 - 2 ** -48;
  const tail = createPhysicalCutStations({ anchors: [{ index: 0, x: 0 }, { index: 2, x: 1 }], edgeSpacing: first });
  assert.equal(tail.x[1], first);
  close(tail.blocks[0].growth / ((1 - first) / first), 1);
  const wide = createPhysicalCutStations({ anchors: [{ index: 0, x: 0 }, { index: 4, x: 1e300 }], edgeSpacing: 1e-320 });
  assert.equal(wide.x[1], 1e-320); assert.equal(wide.x.at(-1), 1e300);
  assert.ok(wide.x.every((v, i) => Number.isFinite(v) && (!i || v > wide.x[i - 1])));
});
