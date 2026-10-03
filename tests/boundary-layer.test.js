import test from 'node:test';
import assert from 'node:assert/strict';
import { hkin, hsl, cfl, dil } from '../src/boundary-layer/tests/closures.js';
import * as reference from '../third_party/vibefoil/js/xblsys.js';
import { solveFlatPlate } from '../src/boundary-layer/tests/flat-plate.js';
import { blasius } from './fixtures/blasius.js';

test('laminar closure values agree with the supplied Vibefoil port across branches', () => {
  for (const h of [1.2, 2.0, 2.59, 3.99, 4, 4.35, 5.49, 5.5, 7, 10]) {
    for (const rt of [10, 300, 2000]) {
      for (const [local, original, key] of [[hsl, reference.hsl, 'hs'], [cfl, reference.cfl, 'cf'], [dil, reference.dil, 'di']]) {
        assert.ok(Math.abs(local(h, rt)[key] - original(h, rt, 0)[key]) < 1e-14);
      }
    }
  }
  assert.deepEqual(hkin(2.6, 0.4), reference.hkin(2.6, 0.4));
});
test('closure analytic derivatives match independent central differences', () => {
  for (const hk of [1.5, 2.59, 3.8, 4.2, 4.7, 5.3, 6.5]) {
    for (const [fn, key, derivative] of [[hsl, 'hs', 'hsHk'], [cfl, 'cf', 'cfHk'], [dil, 'di', 'diHk']]) {
      const h = 1e-5;
      const fd = (fn(hk + h, 700)[key] - fn(hk - h, 700)[key]) / (2 * h);
      assert.ok(Math.abs(fn(hk, 700)[derivative] - fd) < 2e-9, `${key}, Hk=${hk}`);
    }
    for (const [fn, key, derivative] of [[cfl, 'cf', 'cfRt'], [dil, 'di', 'diRt']]) {
      const fd = (fn(hk, 700.01)[key] - fn(hk, 699.99)[key]) / 0.02;
      assert.ok(Math.abs(fn(hk, 700)[derivative] - fd) < 1e-11);
    }
  }
});
test('global laminar flat plate recovers Blasius within correlation accuracy', () => {
  const exact = blasius();
  assert.ok(Math.abs(exact.theta - 0.664114672) < 1e-7);
  // Closure correlations are fits, so agreement is not expected at roundoff.
  for (const reynolds of [1e5, 1e6, 1e7]) {
    const result = solveFlatPlate({ reynolds, initialFactor: 1.7 });
    assert.equal(result.converged, true);
    assert.ok(result.history.length > 2);
    for (const p of result.stations) {
      assert.ok(Math.abs(p.theta / (exact.theta * Math.sqrt(p.x / reynolds)) - 1) < 0.01);
      assert.ok(Math.abs(p.h / exact.h - 1) < 0.01);
      assert.ok(Math.abs(p.cf / (exact.cf / Math.sqrt(reynolds * p.x)) - 1) < 0.01);
    }
  }
});
