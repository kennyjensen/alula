import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';
// Controller seams verify work selection and retained-state reporting, not flow convergence.
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-ncrit-recovery.js', import.meta.url), 'utf8');
let serial = 0;
const nodes = n => [[[{ x: 0, y: n }, { x: 0, y: n + 1 }], [{ x: 1, y: n }, { x: 1, y: n + 1 }]]];
function root(ncrit, nx = 128, policy = {}) {
  const families = { euler: 2e-12, boundaryLayer: 3e-12, edgeMatching: 1e-12 };
  const input = { mach: .2, wakeGeometry: 'independent-banks', wakeDisplacementMotion: 'te-center',
    hybrid: { ismom: 4 }, bodies: [{ leadingIndex: 20, trailingIndex: 100 }], outerLower: Array(nx + 1).fill(0) };
  const options = { ncrit, reynolds: 919616.9, transitionMode: 'automatic', tripFractions: [[1, 1]], transitionState: [24, 17],
    hkFloorLinearization: policy.hkFloorLinearization ?? 'native' };
  const checkpoint = { version: 1, families, restart: { input, options,
    initialEuler: { x: [ncrit], nodes: nodes(ncrit), undisplacedNodes: nodes(ncrit) }, initialBL: [1, 2, 3, 4] },
    continuation: { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', linearOrdering: 'station-auto',
      blUpdate: 'xfoil', ...(policy.legacy ? {} : { shearCoordinate: policy.shearCoordinate ?? 'linear' }) } };
  return { converged: true, families, checkpoint, conditions: { mach: .2, ncrit, reynolds: options.reynolds },
    reason: 'residual', history: [{ iteration: 0 }, { iteration: 1 }],
    flow: { nodes: nodes(ncrit), undisplacedNodes: nodes(ncrit), cells: [[[{ interfacePressure: { lower: ncrit, upper: ncrit } }]]] },
    mesh: { quality: { valid: true }, cells: Array(nx * 3).fill([0, 1, 2, 3]),
      initialization: { gridSmoothing: { enabled: true } } } };
}
const controls = { startupAttempt: 1, maxStartupAttempts: 2, maxIterations: 40, thicknessFactor: 1 };
async function harness({ stop, badReplay = false } = {}) {
  const calls = [], key = `__ncritRecovery${++serial}`;
  const stubs = { observableFlow, observableBL, ...shearPolicy,
    createCoupledStreamtubeBody(input, options) {
      const ncrit = options.ncrit;
      return { initial: [ncrit], euler: { layout: { nx: input.outerLower.length - 1, tubes: [3, 3], bodies: input.bodies } },
        bl: { scale: .001, wakes: [], stations: [{ id: 0, kind: 'surface', body: 0, side: 'upper', i: 54 }],
          surfaces: [{ body: 0, side: 'upper', ids: [0] }] },
        evaluate: () => ({ families: { ...root(ncrit).families, ...(badReplay ? { euler: 9 } : {}) },
          layers: { transitions: [{ kind: 'natural', body: 0, side: 'upper', id: 0, s: .13 }] } }) };
    },
    continueCoupledNcrit(value, options) {
      calls.push({ kind: 'continue', source: value, options });
      const failed = stop === options.targetNcrit, next = failed ? value : root(options.targetNcrit, value.checkpoint.restart.input.outerLower.length - 1, {
        legacy: value.checkpoint.continuation.shearCoordinate === undefined,
        shearCoordinate: value.checkpoint.continuation.shearCoordinate,
        hkFloorLinearization: value.checkpoint.restart.options.hkFloorLinearization });
      options.onStage?.({ stage: 'coupled', actualNcrit: next.conditions.ncrit });
      options.onIteration?.({ stage: 'coupled', iteration: 1, actualNcrit: next.conditions.ncrit, targetNcrit: options.targetNcrit });
      options.onIterationCheckpoint?.(next.checkpoint, { retained: failed });
      const iteration = Object.freeze({ iteration: 1, actualNcrit: next.conditions.ncrit, targetNcrit: options.targetNcrit,
        ncritContinuation: Object.freeze({ actualNcrit: next.conditions.ncrit, targetNcrit: options.targetNcrit, stateConverged: failed }) });
      options.onMesh?.(Object.freeze({ nodes: next.flow.nodes, iteration }), failed ? 'solving' : 'initial', 'coupled');
      options.onFlow?.({ checkpoint: next.checkpoint, flow: next.flow });
      return { result: next, diagnostics: { reachedTarget: !failed, actualNcrit: next.conditions.ncrit,
        attempts: [{ iterations: 2 }] } };
    },
    recoverCoupledTransition(value, options) {
      calls.push({ kind: 'refine', source: value, options });
      const next = root(value.conditions.ncrit, options.plan.refinedNx, {
        legacy: value.checkpoint.continuation.shearCoordinate === undefined,
        shearCoordinate: value.checkpoint.continuation.shearCoordinate,
        hkFloorLinearization: value.checkpoint.restart.options.hkFloorLinearization });
      if (stop === 'checkpointless' || stop === 'checkpointless-residual-root') {
        delete next.checkpoint; next.converged = stop === 'checkpointless-residual-root';
        next.reason = 'Coupled ISES initial redistribution rejected: invalid corner';
        next.initialRedistribution = { accepted: false, rejection: 'invalid corner' };
        next.lastRejectedStep = { code: 'streamtube-grid-nonconvex', diagnostics: { group: 1, i: 17, tube: 0 } };
        return next;
      }
      options.onIteration?.({ iteration: 1, stage: undefined });
      options.onIterationCheckpoint?.(next.checkpoint, { stage: 'coupled' });
      options.onMesh?.(Object.freeze({ nodes: next.flow.nodes, iteration: Object.freeze({ iteration: 1, stage: undefined }) }), 'solving', 'coupled');
      options.onFlow?.({ checkpoint: next.checkpoint, flow: next.flow });
      if (stop === 'refine') next.converged = false;
      if (stop === 'prepare') throw new Error('Native profile rejected');
      return next;
    },
    streamtubeMeshSnapshot: frame => ({ initialization: {}, nodes: frame.nodes, iteration: frame.iteration }),
  };
  globalThis[key] = stubs;
  const text = source.replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis.${key};`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(`${text}\n//# sourceURL=ncrit-recovery-${serial}.js`).toString('base64'));
  delete globalThis[key];
  return { ...module, calls };
}
function failed() { const r = root(9, 128, { legacy: true, hkFloorLinearization: 'exact' }); r.converged = false; r.families.boundaryLayer = .4; return r; }

test('only the failed first admissible automatic terminal-trip attempt selects a thin lower-Ncrit startup', async () => {
  const h = await harness(), r = failed(), input = { ncrit: 9, mach: .2 }, before = structuredClone(r);
  const plan = h.coupledNcritStartupRecoveryPlan(input, r, controls);
  assert.equal(plan.sourceNcrit, 4); assert.equal(plan.targetNcrit, 9); assert.equal(plan.thicknessFactor, .25);
  assert.equal(plan.coarseAdvanceNcrit, 5.5); assert.equal(plan.refinedIterations, 20); assert.equal(plan.equationsChanged, false);
  assert.deepEqual(r, before);
  for (const mutate of [x => { x.converged = true; }, x => { delete x.checkpoint; },
    x => { x.mesh.quality.valid = false; }, x => { x.families.euler = NaN; },
    x => { x.checkpoint.restart.options.transitionMode = 'fixed-trip'; },
    x => { x.checkpoint.restart.options.tripFractions = []; },
    x => { x.checkpoint.restart.options.tripFractions = [[.1, 1]]; },
    x => { x.checkpoint.restart.input.wakeDisplacementMotion = 'fixed'; },
    x => { x.automaticRefinement = { kind: 'transition-local' }; }]) {
    const candidate = failed(); mutate(candidate);
    assert.equal(h.coupledNcritStartupRecoveryPlan(input, candidate, controls), null);
  }
  for (const patch of [{ maxIterations: 0 }, { maxStartupAttempts: 1 }, { startupAttempt: 2 },
    { thicknessFactor: 0 }, { thicknessFactor: Infinity }, { coarseInitialization: {} }])
    assert.equal(h.coupledNcritStartupRecoveryPlan(input, failed(), { ...controls, ...patch }), null);
  for (const input of [undefined, { ncrit: 4, mach: .2 }, { ncrit: 9, mach: .3 }, { ncrit: Infinity, mach: .2 }])
    assert.equal(h.coupledNcritStartupRecoveryPlan(input, failed(), controls), null);
});

test('natural-transition windows use true layout intervals and enforce the entire node budget before refinement', async () => {
  const h = await harness(), r = root(5), before = structuredClone(r);
  const plan = h.coupledStartupTransitionPlan(r);
  assert.equal(plan.parentNx, 128); assert.equal(plan.refinedNx, 143);
  assert.deepEqual(plan.surfaces[0].intervals, [51, 52, 53, 54, 55]);
  assert.equal(plan.nodeCount, 144 * 8); assert.deepEqual(r, before);
  assert.throws(() => h.coupledStartupTransitionPlan(r, { maxNodes: 144 * 8 - 1 }), /exceeding/);
  assert.throws(() => h.coupledStartupTransitionPlan(r, { tolerance: NaN }), /strict complete root/);
  const bad = await harness({ badReplay: true });
  assert.throws(() => bad.coupledStartupTransitionPlan(r), /replay exactly/);
  assert.equal(h.calls.length, 0); assert.equal(bad.calls.length, 0);
});

test('startup advances complete physical roots, prepares native local refinement, then retains the existing final continuation policy', async () => {
  const h = await harness(), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
  const source = root(4), before = structuredClone(source), packets = [];
  const output = h.finishCoupledNcritStartup(source, { plan, maxIterations: 40, tolerance: 1e-10, startupAttempt: 2,
    onFlow: frame => packets.push(frame), onIterationCheckpoint: (cp, detail) => {
      assert.equal(detail.actualNcrit, cp.restart.options.ncrit); assert.equal(detail.targetNcrit, 9);
    } });
  assert.deepEqual(h.calls.map(x => [x.kind, x.source.conditions.ncrit]), [['continue', 4], ['refine', 5.5], ['continue', 5.5]]);
  assert.equal(h.calls[1].options.blPredictor, 'xfoil-mrchdu'); assert.equal(h.calls[1].options.maxIterations, 20);
  assert.equal(h.calls[2].options.maximumStep, undefined);
  assert.equal(output.diagnostics.refinement.actualIntervals, 143);
  assert.equal(output.result.checkpoint.restart.input.outerLower.length - 1, 143);
  assert.equal(output.diagnostics.reachedTarget, true); assert.equal(output.diagnostics.actualNcrit, 9);
  assert.equal(output.diagnostics.iterations, 6, 'Count the source, coarse advance, refined solve and final continuation.');
  assert.deepEqual(source, before);
  assert.ok(packets.every(frame => frame.actualNcrit === frame.checkpoint.restart.options.ncrit && frame.targetNcrit === 9));
});

test('failed final continuation keeps the refined intermediate root and reports the requested criterion separately', async () => {
  const h = await harness({ stop: 9 }), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
  const frames = [];
  const output = h.finishCoupledNcritStartup(root(4), { plan, maxIterations: 40, tolerance: 1e-10,
    onFlow: frame => frames.push(frame) });
  assert.equal(output.result.converged, true); assert.equal(output.result.conditions.ncrit, 5.5);
  assert.equal(output.result.checkpoint.restart.input.outerLower.length - 1, 143);
  assert.equal(output.diagnostics.reachedTarget, false); assert.equal(output.diagnostics.stateConverged, true);
  assert.equal(output.diagnostics.actualNcrit, 5.5); assert.equal(output.diagnostics.targetNcrit, 9);
  assert.equal(frames.at(-1).actualNcrit, 5.5); assert.equal(frames.at(-1).targetNcrit, 9);
  assert.equal(h.calls.length, 3);
});

test('startup completion rejects a missing budget or nonroot intermediate before any work', async () => {
  const h = await harness(), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
  for (const [state, patch] of [[failed(), {}], [root(5), {}], [root(4), { maxIterations: 0 }],
    [root(4), { tolerance: undefined }], [root(4), { tolerance: NaN }]])
    assert.throws(() => h.finishCoupledNcritStartup(state, { plan, maxIterations: 40, tolerance: 1e-10, ...patch }), /converged intermediate source/);
  assert.equal(h.calls.length, 0);
});

test('failed refinement retains the matching complete root and restores its pressure/grid topology', async () => {
  for (const stop of ['refine', 'prepare']) {
    const h = await harness({ stop }), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
    const frames = [], meshes = [], checkpoints = [];
    const output = h.finishCoupledNcritStartup(root(4), { plan, maxIterations: 40, tolerance: 1e-10,
      onFlow: x => frames.push(x), onMesh: x => meshes.push(x), onIterationCheckpoint: cp => checkpoints.push(cp) });
    assert.equal(output.diagnostics.reachedTarget, false); assert.equal(output.diagnostics.stateConverged, true);
    assert.equal(output.result.conditions.ncrit, 5.5); assert.equal(output.result.checkpoint.restart.input.outerLower.length, 129);
    assert.equal(h.calls.length, 2); assert.ok(output.diagnostics.failure);
    assert.deepEqual(frames.at(-1).checkpoint, checkpoints.at(-1));
    assert.deepEqual(frames.at(-1).flow.nodes, meshes.at(-1).nodes);
    assert.equal(frames.at(-1).actualNcrit, 5.5); assert.equal(frames.at(-1).targetNcrit, 9);
  }
});

test('startup observer cancellation propagates all thrown values without subsequent work', async () => {
  for (const key of ['onStage', 'onIteration', 'onMesh', 'onFlow', 'onIterationCheckpoint']) for (const value of [undefined, new Error('cancel')]) {
    const h = await harness(), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
    let caught = false;
    try { h.finishCoupledNcritStartup(root(4), { plan, maxIterations: 40, tolerance: 1e-10, [key]: () => { throw value; } }); }
    catch (error) { caught = true; assert.equal(error, value); }
    assert.equal(caught, true); assert.ok(h.calls.length <= 2);
  }
});

test('checkpointless refinement keeps its original rejection and restores only the complete strict source', async () => {
  for (const stop of ['checkpointless', 'checkpointless-residual-root']) {
  const h = await harness({ stop });
  const plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
  const frames = [], checkpoints = [], source = root(4), before = structuredClone(source);
  const output = h.finishCoupledNcritStartup(source, { plan, maxIterations: 40, tolerance: 1e-10,
    onFlow: frame => frames.push(frame), onIterationCheckpoint: cp => checkpoints.push(cp) });
  assert.equal(output.result, h.calls[1].source);
  assert.equal(output.result.converged, true); assert.equal(output.result.conditions.ncrit, 5.5);
  assert.equal(output.result.checkpoint.restart.input.outerLower.length, 129);
  assert.equal(output.diagnostics.actualNcrit, 5.5); assert.equal(output.diagnostics.reachedTarget, false);
  assert.equal(output.diagnostics.refinement.checkpointAvailable, false);
  assert.equal(output.diagnostics.refinement.converged, false);
  assert.equal(output.diagnostics.refinement.actualIntervals, undefined);
  assert.equal(output.diagnostics.failure.code, 'streamtube-grid-nonconvex');
  assert.match(output.diagnostics.failure.message, /initial redistribution rejected: invalid corner/);
  assert.deepEqual(output.diagnostics.failure.diagnostics.initialRedistribution,
    { accepted: false, rejection: 'invalid corner' });
  assert.deepEqual(output.diagnostics.failure.diagnostics.lastRejectedStep.diagnostics, { group: 1, i: 17, tube: 0 });
  assert.equal(output.diagnostics.failure.diagnostics.checkpointAvailable, false);
  assert.deepEqual(frames.at(-1).checkpoint, checkpoints.at(-1));
  assert.deepEqual(checkpoints.at(-1), output.result.checkpoint);
  assert.equal(frames.at(-1).actualNcrit, 5.5); assert.equal(frames.at(-1).targetNcrit, 9);
  assert.equal(h.calls.length, 2); assert.deepEqual(source, before);
  }
});


test('every recovery iteration and rollback mesh reports the requested target and the matching checkpoint condition', async () => {
  for (const stop of [undefined, 5.5, 9, 'refine', 'prepare']) {
    const h = await harness({ stop }), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
    const packets = [], iterations = [];
    let checkpoint;
    h.finishCoupledNcritStartup(root(4), { plan, maxIterations: 40, tolerance: 1e-10, startupAttempt: 2,
      onIteration: h => iterations.push(h),
      onIterationCheckpoint: value => { checkpoint = value; },
      onMesh: (mesh, state, stage) => {
        assert.ok(checkpoint, 'The full checkpoint precedes its mesh.');
        assert.equal(mesh.iteration.actualNcrit, checkpoint.restart.options.ncrit);
        assert.equal(mesh.iteration.ncritContinuation.actualNcrit, checkpoint.restart.options.ncrit);
        assert.equal(mesh.iteration.targetNcrit, 9);
        assert.equal(mesh.iteration.ncritContinuation.targetNcrit, 9);
        assert.equal(mesh.iteration.ncritContinuation.reachedTarget, checkpoint.restart.options.ncrit === 9);
        assert.equal(mesh.iteration.stage, 'coupled');
        assert.equal(stage, 'coupled');
        assert.ok(['initial', 'solving'].includes(state));
        assert.deepEqual(mesh.nodes, checkpoint.restart.initialEuler.nodes);
        packets.push({ mesh, checkpoint });
      } });
    assert.ok(packets.length);
    assert.ok(iterations.length);
    for (const iteration of iterations) {
      assert.equal(iteration.stage, 'coupled');
      assert.equal(iteration.targetNcrit, 9);
      assert.equal(iteration.ncritContinuation.targetNcrit, 9);
      assert.ok(Number.isFinite(iteration.actualNcrit));
    }
    if (stop === undefined) assert.deepEqual(iterations.map(h => h.actualNcrit), [5.5, 5.5, 9]);
    if (stop === 5.5) assert.ok(packets.every(p => p.mesh.iteration.actualNcrit === 4), 'Coarse failure restores N4, never the abandoned N5.5 label.');
  }
});


test('lower-Ncrit recovery starts linear, selects native Hk, and preserves the explicit Hk opt-out', async () => {
  const h = await harness(), source = failed(), before = structuredClone(source);
  const plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, source, controls);
  assert.equal(plan.sourceShearCoordinate, 'linear');
  assert.equal(plan.shearCoordinate, 'linear');
  assert.equal(plan.sourceHkFloorLinearization, 'exact');
  assert.equal(plan.hkFloorLinearization, 'native');
  const opted = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2, coupledNativeHk: false }, source, controls);
  assert.equal(opted.hkFloorLinearization, 'exact');
  assert.equal(opted.shearCoordinate, 'linear');
  const output = h.finishCoupledNcritStartup(root(4, 128, { hkFloorLinearization: 'exact' }),
    { plan: opted, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(output.diagnostics.reachedTarget, true);
  for (const call of h.calls) {
    assert.equal(call.source.checkpoint.restart.options.hkFloorLinearization, 'exact');
    assert.equal(call.source.checkpoint.continuation.shearCoordinate, 'linear');
  }
  assert.deepEqual(source, before);
});

test('startup policy disagreement rejects before work; a legacy plan and checkpoint retain linear exact semantics', async () => {
  const h = await harness(), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
  for (const [source, patch, pattern] of [
    [root(4, 128, { shearCoordinate: 'logarithmic' }), {}, /shear coordinate/],
    [root(4, 128, { hkFloorLinearization: 'exact' }), {}, /Hk-floor/],
    [root(4), { shearCoordinate: null }, /Invalid/],
    [root(4), { shearCoordinate: 'wrong' }, /Invalid/],
  ]) assert.throws(() => h.finishCoupledNcritStartup(source,
    { plan: { ...plan, ...patch }, maxIterations: 40, tolerance: 1e-10 }), pattern);
  assert.equal(h.calls.length, 0);
  const legacy = { ...plan }; delete legacy.shearCoordinate; delete legacy.hkFloorLinearization;
  const source = root(4, 128, { legacy: true, hkFloorLinearization: 'exact' });
  const before = structuredClone(source);
  const output = h.finishCoupledNcritStartup(source, { plan: legacy, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(output.diagnostics.reachedTarget, true);
  assert.equal(output.result.checkpoint.continuation.shearCoordinate, undefined);
  assert.equal(output.result.checkpoint.restart.options.hkFloorLinearization, 'exact');
  for (const call of h.calls) {
    assert.equal(call.source.checkpoint.continuation.shearCoordinate, undefined);
    assert.equal(call.source.checkpoint.restart.options.hkFloorLinearization, 'exact');
  }
  assert.deepEqual(source, before);
});


test('a completed logarithmic recovery carries its explicit coordinate through Ncrit continuation and refinement', async () => {
  const h = await harness();
  const initialPlan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
  assert.equal(initialPlan.shearCoordinate, 'linear');
  const plan = { ...initialPlan, shearCoordinate: 'logarithmic', shearRecovery: { kind: 'coupled-logarithmic-shear-recovery' } };
  const source = root(4, 128, { shearCoordinate: 'logarithmic' }), before = structuredClone(source);
  const output = h.finishCoupledNcritStartup(source, { plan, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(output.diagnostics.reachedTarget, true);
  assert.equal(output.diagnostics.shearCoordinate, 'logarithmic');
  assert.equal(output.result.checkpoint.continuation.shearCoordinate, 'logarithmic');
  assert.equal(h.calls.length, 3);
  for (const call of h.calls) {
    assert.equal(call.source.checkpoint.continuation.shearCoordinate, 'logarithmic');
    assert.equal(call.source.checkpoint.restart.options.hkFloorLinearization, 'native');
    assert.equal(Object.hasOwn(call.source.checkpoint.restart.options, 'shearCoordinate'), false);
  }
  assert.deepEqual(source, before);
});

test('the new linear startup plan accepts a legacy omitted coordinate with matching native Hk', async () => {
  const h = await harness(), plan = h.coupledNcritStartupRecoveryPlan({ ncrit: 9, mach: .2 }, failed(), controls);
  const source = root(4, 128, { legacy: true }), before = structuredClone(source);
  const output = h.finishCoupledNcritStartup(source, { plan, maxIterations: 40, tolerance: 1e-10 });
  assert.equal(output.diagnostics.reachedTarget, true);
  assert.equal(output.diagnostics.shearCoordinate, 'linear');
  assert.equal(output.result.checkpoint.continuation.shearCoordinate, undefined);
  assert.deepEqual(source, before);
});
