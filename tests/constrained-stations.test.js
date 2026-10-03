import test from 'node:test';
import assert from 'node:assert/strict';
import { fitConstrainedStations } from '../src/numerics/constrained-stations.js';
const uniform = n => Array.from({ length: n + 1 }, (_, i) => i / n);
const close = (a, b) => assert.ok(Math.abs(a - b) < 2e-9, `${a} versus ${b}`);

test('inactive growth bounds recover the analytic equality-constrained minimum', () => {
  const f = fitConstrainedStations({ positions: uniform(5), firstSpacing: .15, lastSpacing: .1, maximumGrowth: 3 });
  [0, .15, .4, .65, .9, 1].forEach((s, i) => close(f.positions[i], s));
  assert.ok(f.primalResidual < 1e-9); assert.ok(f.stationarityResidual < 1e-12); assert.ok(f.complementarity < 1e-9);
});

test('an active growth constraint has an independent unique four-interval solution', () => {
  const f = fitConstrainedStations({ positions: uniform(4), firstSpacing: .1, lastSpacing: .3, maximumGrowth: 2 });
  [0, .1, .3, .7, 1].forEach((s, i) => close(f.positions[i], s));
  assert.ok(f.complementarity < 1e-9); assert.ok(f.stationarityResidual < 1e-12);
});

test('the geometric envelopes prove an insufficient interval count before any iteration', () => {
  assert.throws(() => fitConstrainedStations({ positions: uniform(13), firstSpacing: .041050489985241016,
    lastSpacing: .006935052461773304, maxSweeps: 1 }), /infeasible.*total bounds/);
  assert.throws(() => fitConstrainedStations({ positions: uniform(3), firstSpacing: .8, lastSpacing: .1 }), /infeasible/);
});

test('reflection preserves the optimum; every actual interval meets growth and endpoint constraints', () => {
  const positions = [0, .03, .07, .15, .3, .5, .7, .86, .95, 1];
  const a = fitConstrainedStations({ positions, firstSpacing: .06, lastSpacing: .08 });
  const b = fitConstrainedStations({ positions: positions.toReversed().map(s => 1 - s), firstSpacing: .08, lastSpacing: .06 });
  a.positions.forEach((s, i) => close(s, 1 - b.positions.at(-i - 1)));
  const h = a.positions.slice(1).map((s, i) => s - a.positions[i]);
  assert.ok(h.every(v => v > 0)); close(h[0], .06); close(h.at(-1), .08);
  for (let i = 1; i < h.length; i++) assert.ok(Math.max(h[i] / h[i - 1], h[i - 1] / h[i]) <= 1.5 + 2e-9);
  assert.throws(() => fitConstrainedStations({ positions, firstSpacing: .06, lastSpacing: .08, maxSweeps: 1 }), /did not converge/);
});
