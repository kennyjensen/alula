import test from 'node:test';
import assert from 'node:assert/strict';
import { createVinokurStations } from '../src/numerics/vinokur-stations.js';

// Independent direct published expression, for moderate parameters where
// subtractive cancellation is harmless. Endpoint data are evaluated from it.
function reference(t, beta, A, hyperbolic) {
  const f = hyperbolic ? Math.tanh : Math.tan;
  const u = beta === 0 ? t : .5 * (1 + f(beta * (t - .5)) / f(beta / 2));
  return u / (A + (1 - A) * u);
}
test('discrete calibration recovers the published two-ended family on both branches', () => {
  for (const [beta, A, hyperbolic] of [[3, .7, true], [6, 3, true], [1.9, 2, false], [2.8, .3, false], [0, 2, true]]) {
    const n = 32, length = 2.4, firstSpacing = length * reference(1 / n, beta, A, hyperbolic);
    const lastSpacing = length * (1 - reference(1 - 1 / n, beta, A, hyperbolic));
    const r = createVinokurStations({ length, intervals: n, firstSpacing, lastSpacing });
    r.positions.forEach((s, i) => assert.ok(Math.abs(s - length * reference(i / n, beta, A, hyperbolic)) < 3e-12));
    assert.ok(Math.abs(r.parameter - beta) < 2e-6);
  }
});
test('end cell lengths are matched, with strict ordering, through uniform and near-uniform limits', () => {
  for (const factor of [1, 1 - Number.EPSILON, 1 + Number.EPSILON, .01, .2, 1.9, 9]) {
    const n = 32, a = factor / n, b = .7 / n;
    const r = createVinokurStations({ length: 1, intervals: n, firstSpacing: a, lastSpacing: b });
    assert.ok(Math.abs(r.positions[1] - a) < 2e-14);
    assert.ok(Math.abs(1 - r.positions.at(-2) - b) < 2e-14);
    assert.ok(r.positions.every((s, i) => !i || s > r.positions[i - 1]));
  }
  const r = createVinokurStations({ length: 1, intervals: 32, firstSpacing: 1 / 32, lastSpacing: 1 / 32 });
  r.positions.forEach((s, i) => assert.ok(Math.abs(s - i / 32) < 1e-15));
});
test('physical scaling and reversal preserve the distribution', () => {
  const options = { length: 1, intervals: 41, firstSpacing: .0002, lastSpacing: .009 };
  const a = createVinokurStations(options);
  const b = createVinokurStations({ ...options, firstSpacing: options.lastSpacing, lastSpacing: options.firstSpacing });
  a.positions.forEach((s, i) => assert.ok(Math.abs(s + b.positions[41 - i] - 1) < 1e-14));
  for (const scale of [1e-100, 1e100]) {
    const b = createVinokurStations({ ...options, length: scale, firstSpacing: options.firstSpacing * scale, lastSpacing: options.lastSpacing * scale });
    a.positions.forEach((s, i) => assert.ok(Math.abs(s - b.positions[i] / scale) < 3e-13));
  }
});
test('strong clustering is evaluated without cancellation', () => {
  const r = createVinokurStations({ length: 1, intervals: 64, firstSpacing: 1e-12, lastSpacing: 2e-12 });
  assert.ok(Math.abs(r.firstSpacing / 1e-12 - 1) < 1e-12);
  assert.ok(Math.abs(r.lastSpacing / 2e-12 - 1) < 1e-4);
  assert.ok(r.positions.every((s, i) => !i || s > r.positions[i - 1]));
});
test('infeasible, unresolved, and invalid requests fail explicitly', () => {
  const base = { length: 1, intervals: 16, firstSpacing: .01, lastSpacing: .01 };
  for (const patch of [{ intervals: 2 }, { intervals: 4.5 }, { length: Infinity }, { firstSpacing: 0 },
    { firstSpacing: .6, lastSpacing: .4 }, { lastSpacing: NaN }, { firstSpacing: 1e-300, lastSpacing: 1e-300 }])
    assert.throws(() => createVinokurStations({ ...base, ...patch }), /Vinokur/);
});
