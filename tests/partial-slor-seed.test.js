// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createPartialSlorSeedRecorder } from '../src/geometry/partial-slor-seed.js';
const stage = 'fixed-request';
const nodes = () => Array.from({ length: 5 }, (_, i) => Array.from({ length: 4 }, (_, j) => ({ x: i, y: j })));
const failure = { origin: 'solver', termination: 'sweep-limit', completedRequestedStages: true, converged: false };
const row = (iteration, residual, merit = residual ** 2) => ({ stage, iteration, residual, merit });

test('retains only a lower-merit same-stage accepted state with nondegraded seed conditioning and independent ownership', () => {
  const seed = nodes(); seed[2][1].x += .3;
  const candidate = nodes(), record = createPartialSlorSeedRecorder({ nodes: seed, stage });
  assert.equal(record.observe(row(0, 3), seed), true);
  assert.equal(record.observe({ ...row(1, 1), stage: 'different-equations' }, candidate), false);
  assert.equal(record.observe(row(1, 1), candidate), true);
  candidate[2][1].x += 10;
  const selected = record.afterRequestedFailure(failure);
  assert.equal(selected.report.iteration, 1); assert.equal(selected.report.smoothingConverged, false);
  assert.equal(selected.report.eligibilityCornerSine, .1);
  assert.deepEqual(selected.nodes, nodes()); selected.nodes[2][1].x -= 10;
  assert.deepEqual(record.afterRequestedFailure(failure).nodes, nodes());
  assert.equal(record.observe(row(1, .1), nodes()), false, 'same or earlier sweep cannot replace a later checkpoint');
  assert.equal(record.observe(row(0, .01), seed), false);
  const folded = nodes(); folded[2][1].y = folded[2][2].y;
  assert.equal(record.observe(row(2, .01), folded), false);
  const movedBoundary = nodes(); movedBoundary[0][0].x += .1;
  assert.equal(record.observe(row(3, .01), movedBoundary), false);
  assert.equal(record.observe(row(4, .01, 100), nodes()), false, 'both measures must decrease');
  assert.equal(record.observe(row(5, NaN), nodes()), false);
  for (const overrides of [{ converged: true }, { completedRequestedStages: false }, { origin: 'observer' },
    { origin: 'exception' }, { termination: 'invalid-grid' }, { termination: 'spacing-nonconvex' }])
    assert.equal(record.afterRequestedFailure({ ...failure, ...overrides }), null);
});

test('a finite measured healthy seed is required; a nearly collapsed seed cannot become a fallback', () => {
  const seed = nodes(), record = createPartialSlorSeedRecorder({ nodes: seed, stage });
  assert.equal(record.afterRequestedFailure(failure), null, 'unobserved residuals cannot qualify');
  assert.equal(record.observe(row(0, 1), seed), true);
  const selected = record.afterRequestedFailure(failure);
  assert.equal(selected.report.iteration, 0); assert.deepEqual(selected.nodes, seed);
  const thin = nodes().map(r => r.map(p => ({ x: p.x + 1e10 * p.y, y: p.y })));
  const unsuitable = createPartialSlorSeedRecorder({ nodes: thin, stage });
  assert.equal(unsuitable.observe(row(0, 1), thin), false);
  assert.equal(unsuitable.afterRequestedFailure(failure), null);
  const invalid = nodes(); invalid[2][1] = { ...invalid[2][2] };
  assert.equal(createPartialSlorSeedRecorder({ nodes: invalid, stage }).afterRequestedFailure(failure), null);
});
