import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { extendWarmBoundaryIncrements } from '../src/euler/streamtube-displacement.js';

let serial = 0;
async function harness({ native = 24, noEvent = false, failure, warnings = [] } = {}) {
  const key = `__transitionProfileTrial${serial++}`, calls = [], ne = 2, count = 28;
  let phase = [25], displacement;
  const x = Float64Array.from([2, 3, ...Array.from({ length: count }, (_, i) => [i < 25 ? 8 : .02, .01, .03, .9]).flat()]);
  const baseNodes = [[[{ x: 0, y: 0 }, { x: 0, y: .4 }, { x: 0, y: 1 }],
    [{ x: 1, y: 0 }, { x: 1, y: .4 }, { x: 1, y: 1 }]]];
  const decoded = { nodes: structuredClone(baseNodes), allocation: { groups: [[{ massFlow: 1 }, { massFlow: 3 }]] } };
  const system = { ne, bl: {
    transitionMode: 'automatic', trips: [[1, 1]], scale: 1,
    stations: Array.from({ length: count }, (_, id) => ({ id })),
    snapshotActive: () => phase.slice(), restoreActive: p => { phase = p.slice(); },
    thicknesses: packed => ({ wakes: [[packed[2]]] }),
    activeTargets: () => [{ body: 0, side: 'upper', from: 25, to: noEvent ? 25 : 27, kind: 'natural' }],
    geometry: () => { calls.push('geometry'); return { coordinates: Array.from({ length: count }, (_, i) => ({ s: .1 + .01 * i })) }; },
    updateActive(packed, state, options) {
      calls.push('normalize'); assert.deepEqual(phase, [native]); assert.equal(options.reinitializeAmplification, true);
      assert.equal(packed[4 * 24], 24 < native ? 8 : .02, 'Incoming aux still carries native meanings.');
      if (failure === 'normalize') throw new Error('deliberate normalization failure');
      for (let i = 0; i < count; i++) packed[4 * i] = i < 25 ? 8 : .03;
      phase = [25]; return { changed: native !== 25 };
    },
  }, euler: {
    conditions: { mach: .2 }, layout: { independentWakeBanks: true, wakeDisplacementMotion: 'te-center' },
    setDisplacement: p => { displacement = structuredClone(p); },
    decode() { calls.push('decode'); const r = structuredClone(decoded); for (const row of r.nodes[0]) row[0].y = .01; return r; },
    adoptGeometry() { throw new Error('A trial must never adopt onto the accepted chart.'); },
  } };
  globalThis[key] = {
    prepareCoupledMrchduProfiles(input) {
      calls.push('predict'); assert.deepEqual(phase, [25]); assert.equal(input.targetMach, .2);
      assert.equal(input.states[25].aux, x[ne + 4 * 25]);
      if (failure === 'predict') throw new Error('deliberate prediction failure');
      const initialBL = input.initialBL.slice();
      for (let i = 0; i < count; i++) initialBL.set([i < native ? 8 : .02, .011, .04, .91], 4 * i);
      return { initialBL, transitionState: [native], diagnostics: { bodies: [{ localConvergenceWarnings: warnings }] } };
    },
    incrementIndependentWakeWidths({ nodes, beforeWidths, afterWidths }) {
      calls.push('wake'); assert.deepEqual(beforeWidths, [[.03]]); assert.deepEqual(afterWidths, [[.04]]);
      if (failure === 'wake') throw new Error('deliberate wake preparation failure');
      return { nodes, diagnostics: { widthIncrement: .01, preservesExistingTangentialSeparation: true } };
    },
    extendWarmBoundaryIncrements(input) { calls.push('extend'); assert.strictEqual(input.sourceNodes, baseNodes); return extendWarmBoundaryIncrements(input); },
  };
  const source = fs.readFileSync(new URL('../src/euler/streamtube-transition-profile-trial.js', import.meta.url), 'utf8')
    .replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(source + '\n//# sourceURL=transition-profile-trial-test.js').toString('base64'));
  delete globalThis[key];
  return { ...module, system, x, baseNodes, decoded, calls, displacement: () => displacement };
}

test('no raw interval change is a bit-exact bypass without native, coordinate or phase work', async () => {
  const h = await harness({ noEvent: true }), before = h.x.slice();
  const p = h.prepareCoupledTransitionProfileTrial(h.system, h.x, h);
  assert.strictEqual(p.x, h.x); assert.strictEqual(p.nodes, h.baseNodes); assert.strictEqual(p.decoded, h.decoded);
  assert.equal(p.diagnostics.active, false); assert.deepEqual(h.calls, []); assert.deepEqual(h.x, before);
});

test('native phase on either side of unchanged tight phase is normalized before returning to old metadata', async () => {
  for (const native of [24, 26]) {
    const h = await harness({ native }), before = h.x.slice(), originalNodes = structuredClone(h.baseNodes);
    const p = h.prepareCoupledTransitionProfileTrial(h.system, h.x, h);
    assert.deepEqual(h.calls, ['geometry', 'predict', 'normalize', 'decode', 'wake', 'extend']);
    assert.deepEqual(p.diagnostics.oldPhase, [25]); assert.deepEqual(p.diagnostics.nativePhase, [native]);
    assert.deepEqual(p.diagnostics.tightPhase, [25]); assert.deepEqual(h.system.bl.snapshotActive(), [25]);
    assert.equal(p.x[2 + 4 * 24], 8); assert.equal(p.x[2 + 4 * 25], .03);
    assert.equal(p.x[2 + 4 * 25 + 1], .011); assert.equal(p.x[2 + 4 * 25 + 2], .04);
    assert.deepEqual(p.x.slice(0, 2), before.slice(0, 2)); assert.deepEqual(h.x, before);
    assert.deepEqual(h.baseNodes, originalNodes); assert.deepEqual(h.displacement(), { wakes: [[.03]] });
    assert.equal(p.nodes[0][0][1].y, .40750000000000003, 'Boundary increment extends from the post-DSLIM geometry.');
  }
});

test('local failures preserve the caller state, old phase and installed displacement', async () => {
  for (const failure of ['predict', 'normalize', 'wake']) {
    const h = await harness({ failure }), before = h.x.slice(), nodes = structuredClone(h.baseNodes);
    assert.throws(() => h.prepareCoupledTransitionProfileTrial(h.system, h.x, h), error => {
      assert.match(error.message, /deliberate/); assert.equal(error.eventProfile.initialGuessOnly, true); return true;
    });
    assert.deepEqual(h.x, before); assert.deepEqual(h.baseNodes, nodes);
    assert.deepEqual(h.system.bl.snapshotActive(), [25]); assert.deepEqual(h.displacement(), { wakes: [[.03]] });
  }
});

test('native warnings remain explicit initial-guess diagnostics', async () => {
  const warnings = ['Convergence failed at native station 7'], h = await harness({ warnings });
  const p = h.prepareCoupledTransitionProfileTrial(h.system, h.x, h);
  assert.deepEqual(p.diagnostics.prediction.bodies[0].localConvergenceWarnings, warnings);
  assert.equal(p.diagnostics.initialGuessOnly, true); assert.equal(p.diagnostics.equationsChanged, false);
});
