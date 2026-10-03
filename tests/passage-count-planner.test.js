import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { planPassageCounts } from '../src/geometry/passage-count-planner.js';
import { createSpacingEnvelope } from '../src/geometry/spacing-envelope.js';
const uniform = n => Array.from({ length: n + 1 }, (_, i) => i / n);
const symmetric = () => ({ anchors: [0, 1, 2, 3], minimumCounts: [16, 8, 16],
  fixedCounts: [{ block: 0, intervals: 16 }, { block: 2, intervals: 16 }], components: [
    { start: 0, end: 1, intervals: 16, boundaries: [{ kind: 'cut', body: 0, end: 'upstream', length: 4, progress: uniform(16) }] },
    ...['upper', 'lower'].map(side => ({ start: 1, end: 2, intervals: 8,
      boundaries: [{ kind: 'wall', body: 0, side, length: 1, range: [0, 1], requested: uniform(8) }] })),
    { start: 2, end: 3, intervals: 16, boundaries: [{ kind: 'cut', body: 0, end: 'wake', length: 4, progress: uniform(16) }] },
  ] });

test('a symmetric resolved outline needs no extra counts and does not sacrifice wall endpoint density', () => {
  const input = symmetric(), before = structuredClone(input), r = planPassageCounts(input);
  assert.deepEqual(r.counts, [16, 8, 16]); assert.equal(r.history.length, 1);
  assert.deepEqual(r.fitted.fitted[1].progress, r.fitted.fitted[2].progress);
  r.fitted.fitted[1].progress.forEach((s, i) => assert.ok(Math.abs(s - i / 8) < 1e-10));
  assert.deepEqual(input, before);
});

test('saved main/flap boundaries choose counts automatically and satisfy actual post-fit resolution and growth', () => {
  const input = JSON.parse(fs.readFileSync(new URL('./fixtures/passage-counts.json', import.meta.url))), r = planPassageCounts(input);
  assert.equal(r.counts[0], 32); assert.equal(r.counts.at(-1), 32);
  assert.ok(r.counts[2] > input.minimumCounts[2]); assert.ok(r.counts[3] > input.minimumCounts[3]);
  assert.ok(r.history.length > 1);
  assert.equal(r.globalMinimumProven, false);
  r.components.forEach((c, k) => {
    assert.equal(c.intervals, r.indices[input.anchors.indexOf(c.end)] - r.indices[input.anchors.indexOf(c.start)]);
    const p = r.fitted.fitted[k].progress, h = p.slice(1).map((v, i) => v - p[i]);
    assert.ok(h.every(v => v > 0));
    for (let i = 1; i < h.length; i++) assert.ok(Math.max(h[i] / h[i - 1], h[i - 1] / h[i]) < 1.5 + 2e-8);
    const walls = c.boundaries.filter(b => b.kind === 'wall');
    if (walls.length) {
      const e = createSpacingEnvelope({ profiles: walls.map(b => b.requested.map(s => (s - b.range[0]) / b.length)), start: 0, end: 1 });
      for (let i = 1; i < p.length; i++) assert.ok(e.metric(p[i]) - e.metric(p[i - 1]) <= 1 + 1e-10);
    }
  });
});

test('physical units preserve counts and bounded failures cannot report success', () => {
  const input = symmetric(), changed = structuredClone(input);
  changed.components.forEach(c => c.boundaries.forEach(b => { b.length *= 7;
    if (b.kind === 'wall') { b.range = b.range.map(v => 5 + 7 * v); b.requested = b.requested.map(v => 5 + 7 * v); }
  }));
  assert.deepEqual(planPassageCounts(changed).counts, planPassageCounts(input).counts);
  assert.throws(() => planPassageCounts({ ...input, maximumTotal: 30 }), /budget/);
  const saved = JSON.parse(fs.readFileSync(new URL('./fixtures/passage-counts.json', import.meta.url)));
  assert.throws(() => planPassageCounts({ ...saved, maximumPasses: 1 }), /did not satisfy/);
});
