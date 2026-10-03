// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeEquationControls } from '../src/euler/streamtube-equation-selection.js';

test('public ISMOM controls are explicit, detached and leave omitted defaults empty', () => {
  assert.deepEqual(streamtubeEquationControls(), {});
  for (let ismom = 1; ismom <= 4; ismom++) {
    const expected = { streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5, ismom },
      upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
    const first = streamtubeEquationControls(ismom);
    assert.deepEqual(first, expected);
    first.hybrid.ismom = 0; first.upwind.boundary.kind = 'changed';
    assert.deepEqual(streamtubeEquationControls(ismom), expected);
  }
});

test('public ISMOM rejects ambiguous or out-of-range values', () => {
  for (const value of [null, false, true, '1', 'auto', 0, 5, 1.5, NaN, Infinity, {}, [1]])
    assert.throws(() => streamtubeEquationControls(value), /ISMOM.*integer/);
});
