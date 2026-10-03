import test from 'node:test';
import assert from 'node:assert/strict';
import { fitLogRatios } from '../src/numerics/log-ratio-fit.js';
const close = (a, b) => assert.ok(Math.abs(a - b) < 2e-13 * Math.max(1, Math.abs(b)), `${a} versus ${b}`);

test('a chain matches ratios with the analytic minimum-change geometric mean gauge', () => {
  const r = fitLogRatios({ natural: [1, 1, 1, 7], connections: [{ left: 0, right: 1, ratio: 2 }, { left: 1, right: 2, ratio: 3 }] });
  const scale = 12 ** (-1 / 3); [scale, 2 * scale, 6 * scale, 7].forEach((v, i) => close(r.values[i], v));
  assert.equal(r.exactEquality, true); close(r.stationarityResidual, 0);
});

test('conflicting parallel continuations select the geometric mean ratio, without discarding either edge', () => {
  const r = fitLogRatios({ natural: [1, 1], connections: [{ left: 0, right: 1, ratio: 2 }, { left: 0, right: 1, ratio: 8 }] });
  close(r.values[0], .5); close(r.values[1], 2); close(r.maximumMismatch, 2); close(r.stationarityResidual, 0);
  assert.equal(r.exactEquality, false);
});

test('an inconsistent symmetric cycle has constant fitted values and a nonzero reported mismatch', () => {
  const r = fitLogRatios({ natural: [1, 1, 1], connections: [0, 1, 2].map(left => ({ left, right: (left + 1) % 3, ratio: 2 })) });
  r.values.forEach(v => close(v, 1)); r.residuals.forEach(v => close(v, -Math.log(2)));
  close(r.stationarityResidual, 0); close(r.maximumMismatch, 2);
});

test('connection reversal, permutation and physical units do not change relative fits', () => {
  const natural = [.03, .07, .2], connections = [{ left: 0, right: 1, ratio: 4 }, { left: 0, right: 1, ratio: 2 }, { left: 2, right: 1, ratio: 3 }];
  const a = fitLogRatios({ natural, connections });
  const b = fitLogRatios({ natural: natural.map(v => 1e6 * v), connections: connections.toReversed().map(e => ({ left: e.right, right: e.left, ratio: 1 / e.ratio })) });
  a.values.forEach((v, i) => close(v, b.values[i] / 1e6));
  assert.throws(() => fitLogRatios({ natural: [0], connections: [] }), /positive/);
  assert.throws(() => fitLogRatios({ natural: [1], connections: [{ left: 0, right: 1, ratio: 2 }] }), /valid connections/);
});
