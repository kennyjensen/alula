// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { coupledStreamtubeTripEvents } from '../src/euler/streamtube-coupled.js';
import { xfoilBLStepLimit } from '../src/euler/streamtube-coupled-xfoil-update.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { selectSurfaceTransition, prepareSurfaceTransition } from '../src/viscous/transition-selection.js';

const target = (extra = {}) => ({ from: 1, to: 2, kind: 'natural', fraction: .2, ...extra });

// A local event adapter; only the scalar limit and event handler are real.
// No Euler system, residual/Jacobian, native march or global update is built.
function synthetic({ mode = 'automatic', targets = [target()], base = 1.7 } = {}) {
  const surfaces = targets.map((t, k) => ({ body: k, side: 'upper', transition: t.from,
    ids: Array.from({ length: 4 }, (_, j) => 4 * k + j) }));
  const stations = surfaces.flatMap(s => s.ids.map((id, j) => ({ id, kind: 'surface',
    regime: j < s.transition ? 'laminar' : 'turbulent' })));
  const current = Float64Array.from([.4, ...stations.flatMap(s => [s.regime === 'laminar' ? 0 : .04, 1, 2, base])]);
  let calls = 0, received;
  const system = { ne: 1, n: current.length, bl: { transitionMode: mode, stations, surfaces,
    kernel: { parameters: { ncrit: 9 } }, activeTargets: () => targets,
    snapshotActive: () => surfaces.map(s => s.transition),
    restoreActive: phase => surfaces.forEach((s, i) => { s.transition = phase[i]; }),
    updateActive: (_, __, options) => { calls++; received = options; return { changed: true, changes: targets }; } } };
  return { system, current, candidate: current.slice(), calls: () => calls, received: () => received };
}

function check(f, { accepts = true, blUpdate = 'xfoil', message = /edge velocity/ } = {}) {
  const current = f.current.slice(), before = f.candidate.slice(), phase = f.system.bl.snapshotActive();
  const prepare = () => coupledStreamtubeTripEvents(f.system, { blUpdate }).prepare(f.candidate, f.current);
  if (accepts) { assert.equal(prepare().changed, true); assert.equal(f.calls(), 1); }
  else { assert.throws(prepare, message); assert.equal(f.calls(), 0); }
  assert.deepEqual(f.current, current); assert.deepEqual(f.candidate, before);
  assert.deepEqual(f.system.bl.snapshotActive(), phase);
}

test('real XFOIL scalar-limited positive and negative .375 updates survive event subtraction roundoff', () => {
  let subtractionOvershoot = false;
  for (const base of [1.7, .63]) for (const sign of [-1, 1]) {
    const f = synthetic({ base }), direction = new Float64Array(f.system.n);
    direction[4] = sign * 1.5;
    const before = direction.slice(), limited = xfoilBLStepLimit(f.system, f.current, direction);
    assert.equal(limited.step, .25); assert.equal(limited.limiter.variable, 'edge-speed');
    assert.equal(limited.step * direction[4], sign * .375);
    f.candidate = f.current.map((value, i) => value + limited.step * direction[i]);
    subtractionOvershoot ||= Math.abs(f.candidate[4] - f.current[4]) > .375;
    check(f); assert.deepEqual(direction, before);
  }
  assert.equal(subtractionOvershoot, true, 'A nonbinary base must exercise the operand-roundoff allowance.');
});

test('natural events allow .3 but reject genuine speed excess above .375 in either direction', () => {
  for (const sign of [-1, 1]) {
    const within = synthetic(); within.candidate[4] += sign * .3; check(within);
    const excess = synthetic(); excess.candidate[4] += sign * (.375 + 1e-10); check(excess, { accepts: false });
  }
});

test('legacy, fixed and prescribed forced events retain the .2 speed cap', () => {
  const cases = [
    { mode: 'fixed-trip', targets: [target({ kind: 'forced' })] },
    { targets: [target()], blUpdate: 'giles' },
    { targets: [target({ kind: 'forced' })] },
    { targets: [target({ kind: undefined })] },
    { targets: [target(), target({ kind: 'forced' })] },
  ];
  for (const c of cases) {
    const large = synthetic(c); large.candidate[4] += .3; check(large, { accepts: false, blUpdate: c.blUpdate });
    const small = synthetic(c); small.candidate[4] += .125; check(small, { blUpdate: c.blUpdate });
  }
});

test('automatic terminal fallback shares the natural-event XFOIL speed bound', () => {
  // The terminal fallback is the result of the full natural-onset scan,
  // unlike an independently prescribed material trip. Both are covered by
  // the original XFOIL global Ue limiter before this event adapter runs.
  for (const targets of [[target({ kind: 'trailing-edge', to: 3, fraction: 1 })],
    [target(), target({ kind: 'trailing-edge', to: 3, fraction: 1 })]]) {
    for (const sign of [-1, 1]) {
      const within = synthetic({ targets }); within.candidate[4] += sign * .375; check(within);
      const excess = synthetic({ targets }); excess.candidate[4] += sign * (.375 + 1e-10);
      check(excess, { accepts: false });
    }
  }
});

test('only pre-transfer changed or reconciled targets determine natural speed eligibility', () => {
  // A forced surface with no event must not veto another natural event.
  const unchangedForced = target({ kind: 'forced', from: 1, to: 1 });
  const unrelated = synthetic({ targets: [target(), unchangedForced] });
  unrelated.candidate[4] += .3; check(unrelated);

  for (const trigger of ['selector', 'packed-n']) {
    const natural = target({ from: 1, to: 1, ...(trigger === 'selector' ? { amplificationReconciliation: true } : {}) });
    const same = synthetic({ targets: [natural] }); same.candidate[4] += .3;
    if (trigger === 'packed-n') same.candidate[1] = 9;
    check(same); assert.deepEqual(same.received().reconcileAmplificationSurfaces, [0]);

    // Same-phase forced reconciliation is itself involved in this event.
    const forced = { ...unchangedForced, ...(trigger === 'selector' ? { amplificationReconciliation: true } : {}) };
    const mixed = synthetic({ targets: [target(), forced] }); mixed.candidate[4] += .3;
    if (trigger === 'packed-n') mixed.candidate[1 + 4 * 4] = 9;
    check(mixed, { accepts: false });
  }
});

test('no-event proposals remain exact identity without invoking conversion or event speed checks', () => {
  for (const mode of ['automatic', 'fixed-trip']) for (const blUpdate of ['xfoil', 'giles']) {
    const f = synthetic({ mode, targets: [target({ from: 1, to: 1 })] }); f.candidate[4] += .7;
    const before = f.candidate.slice(), current = f.current.slice(), phase = f.system.bl.snapshotActive();
    assert.deepEqual(coupledStreamtubeTripEvents(f.system, { blUpdate }).prepare(f.candidate, f.current), { changed: false });
    assert.equal(f.calls(), 0); assert.deepEqual(f.candidate, before); assert.deepEqual(f.current, current);
    assert.deepEqual(f.system.bl.snapshotActive(), phase);
  }
});

test('natural speed eligibility does not change the half-thickness event gate', () => {
  for (const slot of [2, 3]) {
    const f = synthetic(); f.candidate[4] += .3; f.candidate[slot] = f.current[slot] * (.5 - 1e-10);
    check(f, { accepts: false, message: /thickness/ });
    const boundary = synthetic(); boundary.candidate[4] += .3; boundary.candidate[slot] *= .5; check(boundary);
  }
});

const native = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/automatic-transition.json', import.meta.url)));
const profile = native.profiles.find(p => p.parameters.ncrit === 9 && !p.tripS);

test('native-derived profile conversion with .3 speed increment preserves every physical primitive and caller state', () => {
  assert.match(native.provenance.scope, /Executed unmodified original Fortran/);
  const fixtureBefore = structuredClone(profile);
  for (const previousIndex of [1, 8]) {
    const kernel = createIntegralKernel({ ...profile.parameters, exactJacobian: true });
    const ids = profile.states.map((_, id) => id), surface = { body: 0, side: 'upper', ids, transition: previousIndex };
    const candidate = Float64Array.from([.4, ...profile.states.flatMap((s, j) =>
      [j < previousIndex ? 0 : .03, s.theta, s.deltaStar, s.ue])]);
    const current = candidate.slice(); ids.forEach(id => { current[1 + 4 * id + 3] -= .3; });
    const source = current.slice(), before = candidate.slice();
    const decode = packed => profile.states.map((s, id) => ({ ...s,
      aux: packed[4 * id], theta: packed[4 * id + 1], deltaStar: packed[4 * id + 2], ue: packed[4 * id + 3] }));
    const states = decode(candidate.subarray(1));
    const expected = prepareSurfaceTransition(kernel, states, { previousIndex });
    assert.equal(expected.index, profile.expected.index); assert.equal(expected.kind, 'natural');
    assert.ok(Math.abs(expected.s - profile.expected.s) < 5e-5);
    let calls = 0;
    const bl = { transitionMode: 'automatic', kernel, surfaces: [surface], snapshotActive: () => [surface.transition],
      restoreActive: phase => { surface.transition = phase[0]; },
      activeTargets: (_, packed) => {
        const t = selectSurfaceTransition(kernel, decode(packed));
        const up = profile.states[t.index - 1].s, down = profile.states[t.index].s;
        return [{ body: 0, side: 'upper', from: surface.transition, to: t.index, kind: t.kind, fraction: (t.s - up) / (down - up) }];
      },
      updateActive: (packed, _, options) => {
        calls++; assert.deepEqual(options.reconcileAmplificationSurfaces, []);
        const p = prepareSurfaceTransition(kernel, decode(packed), { previousIndex: surface.transition });
        p.auxiliary.forEach((aux, id) => { packed[4 * id] = aux; }); surface.transition = p.index;
        return { changed: p.changed, changes: p.converted };
      } };
    const result = coupledStreamtubeTripEvents({ ne: 1, bl }, { blUpdate: 'xfoil' }).prepare(candidate, current);
    assert.equal(calls, 1); assert.equal(result.changed, true);
    const changedRegimes = Array.from({ length: Math.abs(previousIndex - expected.index) },
      (_, j) => Math.min(previousIndex, expected.index) + j);
    // The newly selected mixed endpoint has a new shear equation even when
    // it was already turbulent. Its event record is additional to the
    // laminar/turbulent flips when the interval moves downstream.
    const changedEquations = previousIndex < expected.index ? [...changedRegimes, expected.index] : changedRegimes;
    assert.deepEqual(result.changes.map(c => c.index), changedEquations);
    const mixed = result.changes.find(c => c.index === expected.index);
    assert.equal(mixed.transitionShearInitialization, true);
    assert(Math.abs(mixed.shearInitialization.residual) <= mixed.shearInitialization.tolerance);
    assert.deepEqual(bl.snapshotActive(), [expected.index]); assert.deepEqual(current, source);
    assert.equal(candidate[0], before[0]);
    for (const id of ids) {
      assert.equal(candidate[1 + 4 * id], expected.auxiliary[id]);
      for (const offset of [1, 2, 3]) assert.equal(candidate[1 + 4 * id + offset], before[1 + 4 * id + offset]);
      if (id < expected.index) assert.ok(candidate[1 + 4 * id] >= 0 && candidate[1 + 4 * id] < 9);
      else assert.ok(candidate[1 + 4 * id] > 0);
    }
  }
  assert.deepEqual(profile, fixtureBefore);
});
