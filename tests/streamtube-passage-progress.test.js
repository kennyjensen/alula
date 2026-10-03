// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { repeatedPassageStepLimit } from '../src/euler/streamtube-iteration-progress.js';
const history = () => Array.from({ length: 6 }, (_, iteration) => ({ iteration, step: .01, residual: .08,
  dissipation: { mucon: 1, mcrit: .99 }, residualDecrease: { beforeSquaredNorm: 10, afterSquaredNorm: 9.99 } }));
test('passage recovery detects persistent small steps using paired residual merits', () => {
  const rows = history();
  assert.ok(repeatedPassageStepLimit(rows, 1e-10));
  assert.equal(repeatedPassageStepLimit(rows.slice(1), 1e-10), null);
  assert.equal(repeatedPassageStepLimit(rows, .1), null);
  for (const patch of [{ step: .5 }, { residualDecrease: { beforeSquaredNorm: 10, afterSquaredNorm: 1 } },
    { residualDecrease: undefined }, { dissipation: { mucon: 1, mcrit: .9 } }]) {
    const changed = history(); Object.assign(changed[3], patch);
    assert.equal(repeatedPassageStepLimit(changed, 1e-10), null);
  }
});
