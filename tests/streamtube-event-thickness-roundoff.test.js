// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { coupledStreamtubeTripEvents } from '../src/euler/streamtube-coupled.js';
import { xfoilBLStepLimit } from '../src/euler/streamtube-coupled-xfoil-update.js';

// The real scalar limiter and event gate run on a local adapter. No Euler
// system, geometry initialization, Jacobian or global solve is needed.
function fixture({ slot = 2, thickness = 1.7, mode = 'automatic', event = true } = {}) {
  const target = { body: 0, side: 'upper', from: 1, to: event ? 2 : 1,
    kind: mode === 'automatic' ? 'natural' : 'forced', fraction: .2 };
  const theta = slot === 2 ? thickness : .25 * thickness;
  const delta = slot === 3 ? thickness : 3 * thickness;
  const current = Float64Array.from([.4, 0, theta, delta, 1, .04, theta, delta, 1]);
  let calls = 0;
  const system = { ne: 1, n: current.length, bl: {
    transitionMode: mode,
    kernel: { parameters: { ncrit: 9 } },
    stations: [{ id: 0, regime: 'laminar' }, { id: 1, regime: 'turbulent' }],
    surfaces: [{ ids: [0, 1], transition: 1 }],
    activeTargets: () => [target],
    snapshotActive: () => [1],
    restoreActive: () => {},
    updateActive: () => { calls++; return { changed: true, changes: [target] }; },
  } };
  return { system, current, candidate: current.slice(), calls: () => calls };
}

function prepare(f, blUpdate = 'xfoil') {
  return coupledStreamtubeTripEvents(f.system, { blUpdate }).prepare(f.candidate, f.current);
}

test('real XFOIL half-thickness proposals survive arithmetic roundoff without changing the step or state', () => {
  let roundedBelow = 0;
  for (const slot of [2, 3]) for (const thickness of [1.7, .03]) {
    const f = fixture({ slot, thickness }), direction = new Float64Array(f.system.n);
    direction[slot] = -123.45 * thickness;
    const originalDirection = direction.slice();
    const limited = xfoilBLStepLimit(f.system, f.current, direction);
    assert.equal(limited.limiter.variable, slot === 2 ? 'theta' : 'delta-star');
    assert.equal(limited.limiter.bound, -.5);
    assert.equal(limited.step, -.5 / limited.limiter.normalizedIncrement);
    f.candidate = f.current.map((value, i) => value + limited.step * direction[i]);
    if (f.candidate[slot] < .5 * f.current[slot]) roundedBelow++;
    const before = f.candidate.slice(), current = f.current.slice();
    assert.equal(prepare(f).changed, true);
    assert.equal(f.calls(), 1);
    assert.deepEqual(f.candidate, before);
    assert.deepEqual(f.current, current);
    assert.deepEqual(direction, originalDirection);
    assert.deepEqual(f.system.bl.snapshotActive(), [1]);
  }
  assert.equal(roundedBelow, 4, 'Every case must exercise the former false rejection.');
});

test('exact equality is accepted and genuine excess is rejected at every tested thickness scale', () => {
  for (const mode of ['automatic', 'fixed-trip']) for (const slot of [2, 3]) {
    for (const thickness of [1e-100, 1e-9, .03, 1.7, 1e100]) {
      const equal = fixture({ mode, slot, thickness });
      equal.candidate[slot] = .5 * equal.current[slot];
      assert.equal(prepare(equal).changed, true);

      const excess = fixture({ mode, slot, thickness });
      excess.candidate[slot] = (.5 - 1e-12) * excess.current[slot];
      const before = excess.candidate.slice(), current = excess.current.slice();
      assert.throws(() => prepare(excess), error => {
        assert.equal(error.code, 'COUPLED_TRANSITION_THICKNESS_LIMIT');
        assert.equal(error.diagnostics.station, 0);
        assert.equal(error.diagnostics.variable, slot === 2 ? 'theta' : 'delta-star');
        assert.equal(error.diagnostics.minimum, .5 * thickness);
        assert(error.diagnostics.minimum - error.diagnostics.candidate > error.diagnostics.roundoff);
        return true;
      });
      assert.equal(excess.calls(), 0);
      assert.deepEqual(excess.candidate, before);
      assert.deepEqual(excess.current, current);
      assert.deepEqual(excess.system.bl.snapshotActive(), [1]);
    }
  }
});

test('legacy Giles gate remains strict and unrelated no-event proposals remain exact no-ops', () => {
  const legacy = fixture(), direction = new Float64Array(legacy.system.n);
  direction[2] = -123.45 * legacy.current[2];
  const limited = xfoilBLStepLimit(legacy.system, legacy.current, direction);
  legacy.candidate = legacy.current.map((value, i) => value + limited.step * direction[i]);
  assert(legacy.candidate[2] < .5 * legacy.current[2]);
  assert.throws(() => prepare(legacy, 'giles'), error => {
    assert.equal(error.code, 'COUPLED_TRANSITION_THICKNESS_LIMIT');
    assert.equal(error.diagnostics.roundoff, 0);
    return true;
  });
  assert.equal(legacy.calls(), 0);

  for (const blUpdate of ['xfoil', 'giles']) {
    const unchanged = fixture({ event: false });
    unchanged.candidate[2] *= .49;
    const before = unchanged.candidate.slice();
    assert.deepEqual(prepare(unchanged, blUpdate), { changed: false });
    assert.equal(unchanged.calls(), 0);
    assert.deepEqual(unchanged.candidate, before);
  }
});

test('nonfinite operands cannot turn a thickness violation into an infinite allowance', () => {
  for (const [current, candidate] of [[1.7, -Infinity], [Infinity, 1.7], [1.7, -.1]]) {
    const f = fixture();
    f.current[2] = current;
    f.candidate[2] = candidate;
    const before = f.candidate.slice(), accepted = f.current.slice();
    assert.throws(() => prepare(f), error => {
      assert.equal(error.code, 'COUPLED_TRANSITION_THICKNESS_LIMIT');
      assert(Number.isFinite(error.diagnostics.roundoff));
      if (!Number.isFinite(current) || !Number.isFinite(candidate))
        assert.equal(error.diagnostics.roundoff, 0);
      return true;
    });
    assert.equal(f.calls(), 0);
    assert.deepEqual(f.candidate, before);
    assert.deepEqual(f.current, accepted);
  }
});
