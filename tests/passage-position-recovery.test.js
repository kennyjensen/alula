import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { planPassageCounts } from '../src/geometry/passage-count-planner.js';
import { createSpacingEnvelope } from '../src/geometry/spacing-envelope.js';

const read = name => JSON.parse(fs.readFileSync(new URL(`./fixtures/${name}.json`, import.meta.url)));
const hash = value => createHash('sha256').update(JSON.stringify(value)).digest('hex');

test('successful existing main/flap and THREE station plans retain identical serialized results', () => {
  assert.equal(hash(planPassageCounts(read('passage-counts'))), '36f535880bd8de10c295c67877d57d23207d121fe28323f831ff15080f365ddc');
  assert.equal(hash(planPassageCounts(read('three-passage-counts'))), 'a8c02b10094ef69c91e48319d02c94c59f6dbc60c17668327eb4e6471f0dd4d0');
});

test('30P cove recovery satisfies every original density request and independently bounds position error', () => {
  const input = read('30p-passage-counts'), before = structuredClone(input), r = planPassageCounts(input);
  assert.deepEqual(r.counts, [32, 74, 129, 11, 40, 36, 32]);
  assert.equal(r.placementRecovery.attempted, true);
  assert.equal(r.placementRecovery.originalAttempt.code, 'passage-density-unresolved');
  assert.equal(r.placementRecovery.originalAttempt.history.length, 32);
  assert.deepEqual(input, before);
  for (let k = 0; k < r.components.length; k++) {
    const component = r.components[k], projection = r.fitted.fitted[k].projection;
    const p = projection.positions, h = p.slice(1).map((v, i) => v - p[i]);
    assert.ok(h.every(v => v > 0));
    assert.ok(Math.max(...h.slice(1).map((v, i) => Math.max(v / h[i], h[i] / v))) <= 1.5 + 1e-10);
    const c = projection.forwardCertificate;
    assert.ok(c.physicalPositionErrorBound < 1e-10);
    assert.equal(c.maximumInactiveViolationBound, 0); assert.ok(c.minimumActiveMultiplierBound >= 0);
    const walls = component.boundaries.filter(b => b.kind === 'wall');
    if (walls.length) {
      const e = createSpacingEnvelope({ profiles: walls.map(b => b.requested.map(s => (s - b.range[0]) / b.length)), start: 0, end: 1 });
      assert.ok(Math.max(...p.slice(1).map((v, i) => e.metric(v) - e.metric(p[i]))) <= 1 + 1e-10);
    }
  }
});

test('invalid geometry, count budgets and unresolved growth constraints do not trigger density recovery', () => {
  const input = read('30p-passage-counts');
  for (const changed of [{ ...input, anchors: [0, 0] }, { ...input, maximumTotal: 30 }, { ...input, maximumPasses: 1 }]) {
    assert.throws(() => planPassageCounts(changed), error => !error.placementRecovery);
  }
});
