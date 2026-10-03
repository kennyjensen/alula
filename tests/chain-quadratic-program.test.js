import test from 'node:test';
import assert from 'node:assert/strict';
import { fitConstrainedStations } from '../src/numerics/constrained-stations.js';
const close = (a, b, tolerance = 2e-9) => assert.ok(Math.abs(a - b) < tolerance, `${a} versus ${b}`);
const solve = options => fitConstrainedStations({ ...options, method: 'primal-dual' });

test('tridiagonal primal-dual projection matches analytic inactive and active optima', () => {
  const a = solve({ positions: [0, .2, .4, .6, .8, 1], firstSpacing: .15, lastSpacing: .1, maximumGrowth: 3 });
  [0, .15, .4, .65, .9, 1].forEach((s, i) => close(a.positions[i], s));
  const b = solve({ positions: [0, .08, .25, .47, .7, 1], firstSpacing: .12, lastSpacing: .2 });
  const middle = .22 + .05 * .22 ** 2 / (.22 ** 2 + .23 ** 2);
  [0, .12, .3, .3 + middle, .8, 1].forEach((s, i) => close(b.positions[i], s));
  for (const r of [a, b]) {
    assert.ok(r.sweeps < 20); assert.ok(r.primalResidual < 1e-10);
    assert.ok(r.stationarityResidual < 1e-10); assert.ok(r.complementarity < 1e-9);
  }
});

test('a feasible set with no strictly feasible interior still recovers its unique spacing', () => {
  const r = solve({ positions: [0, .25, .5, .75, 1], firstSpacing: .1, lastSpacing: .3, maximumGrowth: 2 });
  [0, .1, .3, .7, 1].forEach((s, i) => close(r.positions[i], s));
  assert.ok(r.sweeps < 30); assert.ok(r.primalResidual < 1e-10);
});

test('reflected data preserve the optimum, and iteration exhaustion remains a failure', () => {
  const positions = [0, .03, .07, .15, .3, .5, .7, .86, .95, 1];
  const a = solve({ positions, firstSpacing: .06, lastSpacing: .08 });
  const b = solve({ positions: positions.toReversed().map(s => 1 - s), firstSpacing: .08, lastSpacing: .06 });
  a.positions.forEach((s, i) => close(s, 1 - b.positions.at(-i - 1)));
  const oracle = fitConstrainedStations({ positions, firstSpacing: .06, lastSpacing: .08 });
  close(a.objective, oracle.objective, 1e-8);
  a.positions.forEach((s, i) => close(s, oracle.positions[i]));
  assert.throws(() => solve({ positions, firstSpacing: .06, lastSpacing: .08, maxSweeps: 1 }), /did not converge/);
});
