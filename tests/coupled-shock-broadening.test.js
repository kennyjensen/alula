// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { temporaryShockMcrit, maximumAppliedDensityChange, rebaseCoupledMcrit } from '../scripts/validation/coupled-shock-broadening.js';

const serial = v => JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x);
const clone = v => JSON.parse(serial(v));
const original = JSON.parse(fs.readFileSync('docs/coupled-current-profile-preparation/two-element-six-update/initial.json')).checkpoint;
function fixture(mach = .6) {
  const cp = clone(original), f = cp.restart;
  f.input = { ...f.input, mach, streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  f.options.blThermodynamics = 'historical-common-isentrope';
  // Manufacture a thermally admissible profile at the new test gas; the
  // original .2 profile contains wakes too thin for Mach .6. This setup is
  // separate from the threshold-only invariants tested below.
  for (let i = 0; i < f.initialBL.length; i += 4)
    f.initialBL[i + 2] = Math.max(f.initialBL[i + 2], 1.5 * f.initialBL[i + 1]);
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  // This is a new manufactured gas fixture, not a replay of the old .2
  // source. Prepare its phase/auxiliaries once; the rebase helper does neither.
  system.bl.updateActive(system.initial.subarray(system.ne), system.initial.subarray(0, system.ne), { reinitializeAmplification: true });
  f.options.transitionState = system.bl.snapshotActive();
  f.initialBL = Array.from(system.initial.subarray(system.ne));
  const value = system.evaluate(system.initial);
  f.initialEuler = { x: Array.from(system.initial.subarray(0, system.ne)), nodes: clone(value.outer.nodes), undisplacedNodes: clone(value.outer.undisplacedNodes) };
  cp.families = value.families;
  return { cp, maxMach: value.outer.diagnostics.maxMach, residual: value.residual };
}

test('documented broadening formula has exact endpoints, known intermediate values and robust large-d limit', () => {
  for (const target of [.75, .8, .99, 1]) {
    assert.equal(temporaryShockMcrit({ targetMcrit: target, densityChange: 0 }), target);
    let previous = target;
    for (const d of [1e-200, .001, .0375, .075, .15, .3, .6, 1, 10, 1e150, Number.MAX_VALUE]) {
      const actual = temporaryShockMcrit({ targetMcrit: target, densityChange: d });
      assert.ok(Number.isFinite(actual) && actual >= .75 && actual <= previous); previous = actual;
      if (d <= 1) {
        const r = d ** 3 / (.15 * (d ** 2 + (.15 / 4) ** 2));
        assert.ok(Math.abs(actual - (.75 + (target - .75) * Math.exp(-(r ** 2)))) <= 2e-16);
      }
    }
    assert.equal(temporaryShockMcrit({ targetMcrit: target, densityChange: Number.MAX_VALUE }), .75);
  }
  for (const args of [{ targetMcrit: .74, densityChange: 0 }, { targetMcrit: 1.01, densityChange: 0 },
    { targetMcrit: .99, densityChange: -1 }, { targetMcrit: .99, densityChange: Infinity }, { targetMcrit: NaN, densityChange: 1 }])
    assert.throws(() => temporaryShockMcrit(args));
});

test('density norm measures applied physical log-density change only, including tiny increments and decreases', () => {
  const before = [Math.log(2), Math.log(3), 1e20], after = [Math.log(2 * 1.1), Math.log(3 * .75), -1e20];
  assert.ok(Math.abs(maximumAppliedDensityChange(before, after, 2) - .25) < 2e-16);
  assert.equal(maximumAppliedDensityChange([0], [1e-16], 1), Math.expm1(1e-16));
  assert.equal(maximumAppliedDensityChange([0], [-1000], 1), 1);
  assert.equal(maximumAppliedDensityChange([], [], 0), 0);
  assert.throws(() => maximumAppliedDensityChange([0], [NaN], 1), /finite/);
  assert.throws(() => maximumAppliedDensityChange([0], [1000], 1), /representable/);
  assert.throws(() => maximumAppliedDensityChange([0], [0], 2), /prefixes/);
});

test('same-threshold full source replay returns deep-identical checkpoint with all caller arrays detached', t => {
  const f = fixture(), cp = f.cp, before = serial(cp), value = rebaseCoupledMcrit(cp, .99);
  assert.equal(serial(value.checkpoint), before); assert.equal(serial(cp), before);
  assert.ok(value.diagnostics.sameThreshold && value.diagnostics.sourceReplay.packedExact && value.diagnostics.sourceReplay.nodesExact);
  assert.equal(value.diagnostics.n, 349); assert.ok(value.diagnostics.densityCount > 0);
  assert.equal(value.diagnostics.operations.constructors, 1);
  assert.equal(value.diagnostics.operations.nativeMarches, 0); assert.equal(value.diagnostics.operations.globalJacobians, 0);
  value.checkpoint.restart.initialBL[0] += 1;
  assert.equal(serial(cp), before);
  t.diagnostic(JSON.stringify({ maxMach: f.maxMach, sourceFamilies: cp.families, sourceConstruction: 'one supplied constructor/evaluation; no startup solve' }));
});

test('active-filter threshold rebase changes residuals and families while preserving all physical data and replay', () => {
  const { cp } = fixture(), before = serial(cp), broad = rebaseCoupledMcrit(cp, .75);
  assert.ok(broad.diagnostics.targetResidualChanges > 0 && broad.diagnostics.residualMaximumChange > 1e-8);
  assert.notDeepEqual(broad.checkpoint.families, cp.families);
  assert.deepEqual(broad.checkpoint.restart.initialEuler, cp.restart.initialEuler);
  assert.deepEqual(broad.checkpoint.restart.initialBL, cp.restart.initialBL);
  assert.deepEqual(broad.checkpoint.restart.options, cp.restart.options);
  assert.deepEqual(broad.checkpoint.continuation, cp.continuation);
  const restored = clone(broad.checkpoint); restored.families = cp.families; restored.restart.input.upwind.mcrit = .99;
  assert.equal(serial(restored), before);
  const replay = rebaseCoupledMcrit(broad.checkpoint, .75);
  assert.equal(serial(replay.checkpoint), serial(broad.checkpoint));
  assert.deepEqual(replay.residual, broad.residual);
  const targetAgain = rebaseCoupledMcrit(broad.checkpoint, .99);
  assert.equal(serial(targetAgain.checkpoint), before);
  assert.equal(serial(cp), before);
});

test('source mismatch, nonordinary controls and failed thermal domains preserve caller checkpoint', () => {
  const source = fixture().cp;
  for (const [modify, pattern] of [
    [cp => { cp.families.euler += 1; }, /source checkpoint does not replay/],
    [cp => { cp.residual = [0]; }, /standard coupled checkpoint/],
    [cp => { cp.continuation.stepAcceptance = 'residual'; }, /ordinary controls/],
    [cp => { cp.restart.input.streamwiseMode = 'entropy'; }, /hybrid\/upwind/],
    [cp => { cp.restart.initialBL[3] = 10; }, /thermal|enthalpy|BL edge|sonic|inadmissible/i],
  ]) {
    const cp = clone(source); modify(cp); const before = serial(cp);
    assert.throws(() => rebaseCoupledMcrit(cp, .75), pattern); assert.equal(serial(cp), before);
  }
});

test('a target-only filter-domain failure is atomic after successful source replay', () => {
  const cp = fixture().cp;
  // At this fixture's subsonic Mach, the .99 activation is below rounding,
  // while the .75 activation is appreciable. A deliberately excessive test
  // coefficient makes only the target biased transport leave its domain.
  cp.restart.input.upwind.mucon = 1e20;
  const before = serial(cp); let failure;
  try { rebaseCoupledMcrit(cp, .75); } catch (error) { failure = error; }
  assert.ok(failure, 'The manufactured target must fail its strict domain.');
  assert.equal(failure.shockBroadening.stage, 'target threshold');
  assert.equal(failure.shockBroadening.sourceReplay.familiesExact, true);
  assert.equal(failure.shockBroadening.operations.constructors, 2);
  assert.equal(failure.shockBroadening.sourceUnchanged, true);
  assert.equal(serial(cp), before);
});
