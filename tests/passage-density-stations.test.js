import test from 'node:test';
import assert from 'node:assert/strict';
import { reconcilePassageStations } from '../src/geometry/passage-density-stations.js';
const uniform = n => Array.from({ length: n + 1 }, (_, i) => i / n);
const isolated = () => [
  { start: 0, end: 1, intervals: 16, boundaries: [{ kind: 'cut', body: 0, end: 'upstream', length: 4, progress: uniform(16) }] },
  { start: 1, end: 2, intervals: 16, boundaries: [{ kind: 'wall', body: 0, side: 'upper', length: 1, range: [0, 1], requested: uniform(8) }] },
  { start: 1, end: 2, intervals: 16, boundaries: [{ kind: 'wall', body: 0, side: 'lower', length: 1, range: [0, 1], requested: uniform(8) }] },
  { start: 2, end: 3, intervals: 16, boundaries: [{ kind: 'cut', body: 0, end: 'wake', length: 4, progress: uniform(16) }] },
];

test('isolated symmetric foil fits all physical joins while preserving upper/lower symmetry', () => {
  const c = isolated(), before = structuredClone(c), r = reconcilePassageStations({ components: c });
  assert.deepEqual(c, before); assert.deepEqual(r.fitted[1].progress, r.fitted[2].progress);
  for (const e of r.connections) {
    const left = r.fitted[e.left >> 1].progress, right = r.fitted[e.right >> 1].progress;
    const ratio = right[1] / (1 - left.at(-2));
    assert.ok(Math.abs(ratio / e.ratio - 1) < 1e-13);
  }
  assert.equal(r.endpointFit.exactEquality, true);
  for (const f of r.fitted) { assert.equal(f.progress[0], 0); assert.equal(f.progress.at(-1), 1); assert.equal(f.progress.length, 17); }
});

test('station fit is invariant to component order and length units', () => {
  const c = isolated(), a = reconcilePassageStations({ components: c });
  const altered = c.toReversed().map(q => ({ ...q, boundaries: q.boundaries.map(b => ({ ...b, length: b.length * 3,
    ...(b.kind === 'wall' ? { range: b.range.map(s => 7 + 3 * s), requested: b.requested.map(s => 7 + 3 * s) } : {}) })) }));
  const b = reconcilePassageStations({ components: altered });
  a.fitted.forEach((f, i) => f.progress.forEach((u, j) => assert.ok(Math.abs(u - b.fitted.at(-i - 1).progress[j]) < 2e-14)));
});

test('missing reference cuts and disconnected outlines fail before placement', () => {
  const a = isolated(); a[0].boundaries = [];
  assert.throws(() => reconcilePassageStations({ components: a }), /reference progress/);
  const b = isolated(); b[3].start = 2.1;
  assert.throws(() => reconcilePassageStations({ components: b }), /Disconnected/);
});
