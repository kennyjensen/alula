// SPDX-License-Identifier: GPL-2.0-or-later
// Small adapter contracts. Native-profile and full-grid numerical checks are
// saved separately; these tests run no global solve, Jacobian or factorization.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { extendWarmBoundaryIncrements } from '../src/euler/streamtube-displacement.js';
import { assertConvexStreamtubeGrid } from '../src/geometry/streamtube-convex-step.js';

const clone = structuredClone;
const nodes = () => [[[{ x: 0, y: .002 }, { x: 0, y: .4 }, { x: 0, y: 1 }],
  [{ x: 1, y: .002 }, { x: 1, y: .4 }, { x: 1, y: 1 }]]];
let serial = 0;
async function harness(config = {}) {
  const calls = { constructs: 0, profiles: 0, extensions: [], wakeIncrements: [], reconciliations: [] }, key = `__gridProfileTest${serial++}`;
  const prepared = { input: { mach: .2, hybrid: { ismom: 4 }, bodies: [{ element: 0 }] },
    options: { reynolds: 1e6, ncrit: 4, transitionMode: 'automatic', transitionState: [1, 1], tripFractions: [[1, 1]] },
    initialEuler: { x: [Math.log(1.1), 2, .002, .3], nodes: nodes(), undisplacedNodes: nodes() },
    initialBL: [.01, .001, .002, .9], transfer: { targetFamilies: { boundaryLayer: 4 }, targetTransitionState: [1, 1] },
    initialization: { thicknessFactor: 1, transfer: { targetFamilies: { boundaryLayer: 4 } } } };
  globalThis[key] = {
    assertConvexStreamtubeGrid,
    incrementIndependentWakeWidths(input) {
      calls.wakeIncrements.push(clone(input));
      return { nodes: clone(input.nodes), diagnostics: { method: 'tested-wake-increment-spy' } };
    },
    extendWarmBoundaryIncrements(input) { calls.extensions.push(clone(input)); return extendWarmBoundaryIncrements(input); },
    streamtubeMeshSnapshot: ({ nodes }) => ({ nodes: clone(nodes), quality: { valid: true }, initialization: {} }),
    prepareCoupledMrchduProfiles(input) {
      calls.profiles++; assert.deepEqual(input.initialBL, prepared.initialBL);
      assert.equal(input.targetMach, .2); assert.equal(input.bl.transitionMode, 'automatic');
      return { initialBL: [.02, .00105, .0022, .91], transitionState: [2, 2],
        diagnostics: { operations: { translatedMRCHDUCalls: 1 },
          bodies: [{ localConvergenceWarnings: config.warnings ?? [] }] } };
    },
    createCoupledStreamtubeBody(input, options) {
      const number = ++calls.constructs, initial = Float64Array.from([...options.initialEuler.x, ...options.initialBL]);
      let currentNodes = clone(options.initialEuler.nodes), phase = clone(options.transitionState);
      const source = number === 1;
      if (source && config.changeSourceNodes) currentNodes[0][0][1].y += 1e-12;
      const allocation = { groups: [[{ massFlow: 1 }, { massFlow: 3 }]] };
      const conditions = { mach: input.mach, reynolds: options.reynolds, ncrit: options.ncrit, referenceChord: 1 };
      return { initial, n: 8, ne: 4, conditions,
        bl: { transitionMode: options.transitionMode, trips: options.tripFractions,
          kernel: { parameters: { ncrit: options.ncrit, mach: input.mach, reynolds: options.reynolds } },
          snapshotActive: () => phase.slice(),
          updateActive(_state, _euler, controls) { calls.reconciliations.push(clone(controls));
            return { changed: false, changes: [] }; }, thicknesses: x => ({ delta: x[2], wakes: [[x[2]]] }) },
        euler: { conditions: { h0: 63, mach: input.mach, hybrid: clone(input.hybrid) },
          setDisplacement() {},
          decode() { const target = clone(currentNodes); for (const row of target[0]) row[0].y = options.initialBL[2];
            return { nodes: target }; },
          adoptGeometry(state, updated) { currentNodes = clone(updated); const x = state.slice(); x[2] = updated[0][0][0].y; return x; } },
        admissible: () => source || !config.rejectDomain,
        evaluate(state) {
          const rho = !source && config.changeDensity ? 1.2 : Math.exp(state[0]);
          const value = source ? 4 : 1e-6, residual = Float64Array.from([.03, value, .06]);
          if (number === 3 && config.staleReplay) residual[1] += .1;
          return { residual, families: { euler: .03, boundaryLayer: value, edgeMatching: .06 },
            layers: { states: [{ s: .5, aux: state[4], theta: state[5], deltaStar: state[6], ue: state[7] }] },
            outer: { nodes: clone(currentNodes), undisplacedNodes: clone(currentNodes), allocation: clone(allocation),
              captured: [4], stagnation: [.5], strengths: [.3], sections: [[[{ rho }]]] } };
        } };
    },
  };
  const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-grid-profile.js', import.meta.url), 'utf8')
    .replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
  delete globalThis[key];
  return { ...module, prepared, calls };
}

test('native profile correction extends only the wall increment in physical mass coordinates and retains the entire source', async () => {
  const h = await harness(), before = clone(h.prepared), p = h.prepareCoupledGridProfile(h.prepared);
  assert.deepEqual(h.prepared, before); assert.equal(h.calls.profiles, 1); assert.equal(h.calls.extensions.length, 1);
  const middle = p.initialEuler.nodes[0][0][1];
  assert.ok(Math.abs(middle.y - (.4 + .75 * .0002)) < 1e-15, 'Only the new lower-wall increment is extended at captured mass fraction1/4.');
  assert.equal(p.initialEuler.nodes[0][0][0].y, .0022); assert.equal(p.initialEuler.nodes[0][0][2].y, 1);
  assert.deepEqual(h.calls.extensions[0].masses, [[1, 3]]);
  assert.deepEqual(p.options.transitionState, [2, 2]); assert.equal(p.options.ncrit, 4);
  assert.deepEqual(p.input, h.prepared.input); assert.equal(p.initialEuler.x[0], before.initialEuler.x[0]);
  const diagnostics = p.transfer.profilePreparation;
  assert.equal(diagnostics.canonicalReplayExact, true); assert.equal(diagnostics.physicalEulerInvariantsPreserved, true);
  assert.equal(diagnostics.flowSolved, false); assert.equal(diagnostics.converged, false);
  assert.equal(diagnostics.geometry.fullDisplacementReapplied, false); assert.equal(diagnostics.geometry.independentWakeBanksAveraged, false);
  assert.deepEqual(p.initialization.transfer, p.transfer); assert.equal(p.transfer.targetFamilies.boundaryLayer, 1e-6);
  assert.deepEqual(p.transfer.targetTransitionState, [2, 2]); assert.equal(diagnostics.beforeFamilies.boundaryLayer, 4);
  assert.equal(diagnostics.operations.globalNewtonUpdates, 0); assert.equal(diagnostics.operations.globalLinearSolves, 0);
});

test('explicit TE-center preparation passes total width increments and tight amplification reconciliation', async () => {
  const h = await harness();
  Object.assign(h.prepared.input, { wakeGeometry: 'independent-banks', wakeDisplacementMotion: 'te-center' });
  const before = clone(h.prepared), p = h.prepareCoupledGridProfile(h.prepared,
    { wakeWidthIncrement: true, reinitializeAmplification: true });
  assert.deepEqual(h.prepared, before);
  assert.deepEqual(h.calls.reconciliations, [{ reinitializeAmplification: true }]);
  assert.equal(h.calls.wakeIncrements.length, 1);
  assert.deepEqual(h.calls.wakeIncrements[0].beforeWidths, [[.002]]);
  assert.deepEqual(h.calls.wakeIncrements[0].afterWidths, [[.0022]]);
  assert.deepEqual(h.calls.extensions[0].targetNodes, h.calls.wakeIncrements[0].nodes);
  assert.equal(p.transfer.profilePreparation.geometry.wakeWidthIncrement.method, 'tested-wake-increment-spy');
});

test('new profile controls default to the original path and reject unsupported wake charts before construction', async () => {
  const h = await harness(); h.prepareCoupledGridProfile(h.prepared);
  assert.deepEqual(h.calls.reconciliations, [undefined]); assert.deepEqual(h.calls.wakeIncrements, []);
  for (const options of [{ wakeWidthIncrement: null }, { reinitializeAmplification: 'true' }, { wakeWidthIncrement: true }]) {
    const g = await harness();
    assert.throws(() => g.prepareCoupledGridProfile(g.prepared, options), { code: 'COUPLED_GRID_PROFILE_PREPARATION' });
    assert.equal(g.calls.constructs, 0);
  }
});

test('native local warnings reject before geometry changes and preserve explicit warning details', async () => {
  const h = await harness({ warnings: ['station17 did not converge'] }), before = clone(h.prepared);
  assert.throws(() => h.prepareCoupledGridProfile(h.prepared), error => {
    assert.equal(error.code, 'COUPLED_GRID_PROFILE_PREPARATION');
    assert.deepEqual(error.diagnostics.warnings, ['station17 did not converge']);
    assert.equal(error.diagnostics.stage, 'native profile preparation'); return true;
  });
  assert.equal(h.calls.constructs, 1); assert.equal(h.calls.extensions.length, 0); assert.deepEqual(h.prepared, before);
});

test('strict domain, preserved Euler invariants and complete canonical residual replay remain mandatory', async () => {
  for (const [config, message] of [[{ rejectDomain: true }, /physical\/convex domain/],
    [{ changeDensity: true }, /physical Euler density/], [{ staleReplay: true }, /canonical full residual/]]) {
    const h = await harness(config), before = clone(h.prepared);
    assert.throws(() => h.prepareCoupledGridProfile(h.prepared), error => {
      assert.equal(error.code, 'COUPLED_GRID_PROFILE_PREPARATION'); assert.match(error.message, message); return true;
    });
    assert.deepEqual(h.prepared, before);
  }
});

test('incomplete and nonautomatic grid profiles are rejected before numerical construction', async () => {
  const h = await harness();
  for (const p of [undefined, {}, { ...h.prepared, initialBL: undefined },
    { ...h.prepared, options: { ...h.prepared.options, transitionMode: 'fixed-trip' } }])
    assert.throws(() => h.prepareCoupledGridProfile(p), /complete automatic coupled initial guess/);
  assert.equal(h.calls.constructs, 0); assert.equal(h.calls.profiles, 0);
});


test('source-grid replay rejection reports coordinate differences without mutating the caller', async () => {
  const h = await harness({ changeSourceNodes: true }), before = clone(h.prepared);
  assert.throws(() => h.prepareCoupledGridProfile(h.prepared), error => {
    assert.equal(error.diagnostics.failedInvariant, 'source physical nodes');
    const d = error.diagnostics.geometryReplayDifference;
    assert.equal(d.changedCoordinates, 1);
    assert.ok(d.maximumAbsoluteDifference > 0 && d.maximumAbsoluteDifference < 2e-12);
    assert.equal(d.firstDifference.coordinate, 'y');
    assert.equal(d.callerSourceCloned, true);
    return true;
  });
  assert.deepEqual(h.prepared, before);
  assert.equal(h.calls.profiles, 0, 'replay failure occurs before native profile preparation');
});
