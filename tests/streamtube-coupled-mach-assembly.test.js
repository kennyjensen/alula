// SPDX-License-Identifier: GPL-2.0-or-later
// One tiny source and one actual warm target; no public cold mesh/solve.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { coupledMachPlan, solveCoupledStreamtubeMach } from '../src/euler/tests/streamtube-coupled-mach-assembly.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { coupledAssemblyConditions } from '../src/euler/streamtube-coupled-assembly.js';
import { certifyCoupledStreamtubeHybrid } from '../src/euler/tests/streamtube-coupled-hybrid-certification.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const serial = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
let cached, startupError, warm;
function source() {
  if (startupError) throw startupError;
  if (!cached) try {
    const fixture = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 });
    const input = { ...fixture, bodies: fixture.bodies.map((b, i) => ({ ...b, element: 1 - i })), streamwiseMode: 'isentropic' };
    const result = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', transitionMode: 'automatic', tripFractions: [[1, 1], [1, 1]],
      maxIterations: 12, tolerance: 1e-10, stepAcceptance: 'admissible' });
    assert.equal(result.converged, true, result.reason); assert.equal(result.mesh.quality.valid, true);
    assert.ok(result.x.length < 500);
    const sourceCase = { flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .2, alpha: 0,
      elements: fixture.bodies.map(b => ({ points: b.points })).reverse(),
      referenceChord: 2 * result.conditions.referenceChord, reynolds: 2e6, ncrit: 9, transitionMode: 'automatic',
      gridIntervals: 4, gridTubes: 2, gridEllipticSmoothing: false };
    const settings = coupledAssemblyConditions(sourceCase, result.solverInput.bodies, result.conditions.referenceChord);
    assert.equal(settings.options.reynolds, result.conditions.reynolds);
    const parent = { model: 'research-streamtube-euler-bl', converged: true, mesh: { quality: { valid: true } },
      sourceCase, checkpoint: serial(result.checkpoint), families: result.families, ...settings.normalization,
      initialization: { euler: { iterations: 0, gridSmoothing: { enabled: false } },
        boundaryLayer: { method: 'tiny-test-source', thicknessFactor: 1 } },
      solverSettings: { tolerance: 1e-10, streamwiseMode: 'isentropic', edgeMatching: 'section-velocity' } };
    cached = { parent, result };
  } catch (error) { startupError = error; throw error; }
  return structuredClone(cached.parent);
}

test('cheap route planning chooses unchanged cold policy and allows only a Mach change in a compact cache', () => {
  const parent = source(), before = structuredClone(parent), caseData = { ...parent.sourceCase, mach: .4 };
  const plan = coupledMachPlan(caseData, parent);
  assert.equal(plan.route, 'warm-certified'); assert.equal(plan.sourceMach, .2); assert.equal(plan.targetMach, .4);
  assert.deepEqual(plan.settings.elementOrder, [1, 0]); assert.equal(plan.settings.normalization.referenceReynolds, 2e6);
  assert.deepEqual(coupledMachPlan(caseData), { route: 'cold-baseline', targetMach: .4, sourceMach: .2,
    tolerance: 1e-10, epsilonP: 1e-5, upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } });
  assert.equal(coupledMachPlan({ ...caseData, mach: .31 }).sourceMach, .155);
  assert.equal(coupledMachPlan({ ...caseData, mach: .3 }).route, 'cold');
  assert.equal(coupledMachPlan({ ...caseData, mach: .2 }).route, 'cold');
  // A plan never decodes the grid or touches arrays needed only by evaluation.
  const guarded = structuredClone(parent);
  Object.defineProperty(guarded.checkpoint.restart.input, 'outerLower', { get() { throw new Error('Numerical grid accessed'); } });
  assert.equal(coupledMachPlan(caseData, guarded).route, 'warm-certified');
  for (const [key, value] of [['alpha', 1], ['reynolds', 3e6], ['ncrit', 8], ['referenceChord', 1],
    ['transitionMode', 'fixed-trip'], ['materialTrips', [[.2, .3], [.4, .5]]],
    ['gridIntervals', 8], ['gridTubes', 3], ['gridEllipticSmoothing', true], ['maxIterations', 7]])
    assert.throws(() => coupledMachPlan({ ...caseData, [key]: value }, parent), /stale/);
  const geometry = structuredClone(caseData); geometry.elements[0].points[0].x += 1e-8;
  assert.throws(() => coupledMachPlan(geometry, parent), /stale/);
  assert.throws(() => coupledMachPlan(caseData, null), /converged coupled/);
  assert.deepEqual(parent, before);
});

test('normalization, phase, source model, completeness and iteration controls reject before continuation', () => {
  const parent = source(), caseData = { ...parent.sourceCase, mach: .21 };
  const changed = mutate => { const p = structuredClone(parent); mutate(p); return p; };
  for (const [mutate, pattern] of [
    [p => { p.converged = false; }, /converged coupled/],
    [p => { p.mesh.quality.valid = false; }, /valid grid/],
    [p => { delete p.checkpoint.restart.initialEuler.undisplacedNodes; }, /complete continuation/],
    [p => { p.checkpoint.families.euler = .1; }, /residuals differ/],
    [p => { p.referenceReynolds *= 2; }, /normalization/],
    [p => { p.checkpoint.restart.options.reynolds *= 2; }, /reynolds/],
    [p => { p.checkpoint.restart.options.transitionState.pop(); }, /transition map/],
    [p => { p.checkpoint.restart.input.mach = .25; }, /actual Mach/],
    [p => { p.checkpoint.restart.input.streamwiseMode = 'momentum'; }, /isentropic equations/],
    [p => { p.solverSettings.streamwiseMode = 'hybrid'; }, /stale streamwiseMode/],
  ]) assert.throws(() => coupledMachPlan(caseData, changed(mutate)), pattern);
  for (const options of [{ maxStages: 0 }, { maxStages: 257 }, { maxSubdivisions: -1 }, { stageMaxIterations: -1 },
    { blPredictor: 'unknown' }, { maxWakeRecoveries: -1 }, { maxWakeRecoveries: 5 }, { maxWakeRecoveries: 1.5 }, { maxWakeRecoveries: NaN },
    { eulerMaxIterations: -1 }, { maxIterations: -1 }, { maxStartupAttempts: 3 }, { maxBacktracks: 21 }, { onStage: true }])
    assert.throws(() => solveCoupledStreamtubeMach(caseData, { parentResult: parent, ...options }), /iteration|continuation controls/);
  assert.throws(() => coupledMachPlan(caseData, parent, { epsilonP: 0 }), /hybrid\/upwind/);
  assert.throws(() => coupledMachPlan(caseData, parent, { upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'implicit' } } }), /hybrid\/upwind/);
});

test('same-Mach isentropic cache preserves the complete original zero-update result and force-model branch', () => {
  const parent = source(), checkpoint = parent.checkpoint, controls = checkpoint.continuation;
  const expected = solveCoupledStreamtubeIses(undefined, { resume: checkpoint, maxIterations: 0, tolerance: 1e-10,
    iterationGeometry: controls.iterationGeometry, stepAcceptance: controls.stepAcceptance, stagnationLimiter: controls.stagnationLimiter });
  let accepted = 0;
  const result = solveCoupledStreamtubeMach(parent.sourceCase, { parentResult: parent, stageMaxIterations: 0,
    onPrepared: () => { throw new Error('No hybrid preparation at an unchanged isentropic operating point.'); },
    onStage: () => { throw new Error('No certification or target stage at the unchanged operating point.'); },
    onIteration: () => { throw new Error('No new solve at the unchanged operating point.'); },
    onCheckpoint: (value, details) => {
      accepted++; assert.equal(details.kind, 'accepted'); assert.equal(details.reachedTarget, true);
      assert.deepEqual(serial(value), serial(expected.checkpoint));
    } });
  for (const key of ['flow', 'boundaryLayer', 'residual', 'x', 'families', 'conditions', 'checkpoint',
    'initialRedistribution', 'history', 'linearDiagnostics', 'solverInput', 'coupledOptions'])
    assert.deepEqual(serial(result[key]), serial(expected[key]), `Unchanged original ISES ${key}`);
  assert.equal(result.machContinuation.route, 'warm-isentropic'); assert.equal(result.machContinuation.certificationPerformed, false);
  assert.equal(result.solverInput.streamwiseMode, 'isentropic'); assert.equal(result.solverInput.hybrid, undefined);
  assert.equal(result.solverInput.upwind, undefined); assert.equal(result.conditions.blThermodynamics, undefined);
  assert.equal(result.solverSettings.streamwiseMode, 'isentropic'); assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(accepted, 1); assert.equal(result.converged, true);
});

test('one tiny warm public route certifies the source, preserves normalization and reaches the actual requested Mach', t => {
  const parent = source(), before = structuredClone(parent), caseData = { ...parent.sourceCase, mach: .21 };
  let prepared, latest, accepted = 0, meshes = 0; const events = [];
  const result = solveCoupledStreamtubeMach(caseData, { parentResult: parent, stageMaxIterations: 8, maxStages: 4,
    onPrepared: value => {
      prepared = value; assert.equal(value.system.initialization.suppliedBL, true);
      assert.equal(value.system.bl.surfaces.length, 4); assert.equal(value.system.bl.wakes.length, 2);
      assert.equal(value.certification.certified, true); assert.equal(value.certification.maximumSpeedBias, 0);
      assert.equal(coupledMachPlan(caseData, value.parentResult).route, 'warm-hybrid');
      value.checkpoint.restart.initialBL.fill(-99); value.settings.materialTrips[0][0] = -99;
    },
    onStage: event => { events.push(event.stage); assert.equal(event.targetMach, .21); event.targetMach = 99; },
    onIteration: event => { assert.equal(event.actualMach, .21); assert.equal(event.targetMach, .21); events.push('iteration'); },
    onCheckpoint: (checkpoint, details) => {
      latest = checkpoint; assert.equal(checkpoint.version, 1); assert.equal(checkpoint.restart.input.mach, details.actualMach);
      assert.equal(details.targetMach, .21); if (details.kind === 'accepted') accepted++;
      if (details.kind === 'accepted') assert.ok(Math.max(...Object.values(checkpoint.families)) <= 1e-10);
      // The callback may change its copy without corrupting the solver.
      checkpoint.restart.initialBL[0] = -99;
    },
    onMesh: (mesh, phase, state) => {
      meshes++; assert.ok(prepared); assert.ok(latest); assert.ok(mesh.vertices.length);
      assert.equal(state.actualMach, .21); assert.equal(state.targetMach, .21);
      assert.equal(mesh.actualMach, .21); assert.equal(mesh.targetMach, .21);
      assert.equal(state.flow.nodes.length, 3); assert.equal(state.flow.sections.length, cached.result.flow.sections.length);
      assert.ok(state.flow.sections[0][0][0].p > 0);
      assert.deepEqual(mesh.initialization.gridSmoothing, parent.initialization.euler.gridSmoothing);
      assert.equal(phase, state.iteration.iteration === 0 ? 'initial' : 'solving');
      state.flow.nodes[0][0][0].x = 99; mesh.vertices[0].x = 99;
    },
  });
  assert.equal(result.converged, true, result.reason); assert.equal(result.mesh.quality.valid, true);
  assert.equal(result.mach, .21); assert.equal(result.conditions.mach, .21); assert.equal(result.solverInput.mach, .21);
  assert.equal(result.checkpoint.restart.input.mach, .21); assert.equal(result.sourceCase.mach, .21);
  assert.equal(result.actualMach, .21); assert.equal(result.targetMach, .21); assert.equal(result.requestedCase, undefined);
  assert.equal(result.machContinuation.reachedTarget, true); assert.equal(result.physicalAcceptance, false);
  assert.equal(result.solverSettings.streamwiseMode, 'hybrid');
  assert.equal(result.solverSettings.blThermodynamics, 'historical-common-isentrope');
  assert.equal(result.solverSettings.maxWakeRecoveries, 4);
  assert.equal(result.referenceChord, parent.referenceChord); assert.equal(result.kernelReynolds, parent.kernelReynolds);
  assert.deepEqual(result.elementOrder, [1, 0]); assert.deepEqual(result.materialTrips, [[1, 1], [1, 1]]);
  assert.equal(result.boundaryLayer.surfaces.length, 4); assert.equal(result.boundaryLayer.wakes.length, 2);
  assert.equal(result.initialRedistribution.resumed, true); assert.deepEqual(result.initialRedistribution.passages, []);
  assert.equal(accepted, 2); assert.equal(meshes, result.history.length);
  assert.equal(events[0], 'hybrid-certification'); assert.ok(events.includes('coupled-mach'));
  assert.deepEqual(parent, before); assert.deepEqual(caseData, { ...before.sourceCase, mach: .21 });
  warm = result;
  t.diagnostic(JSON.stringify({ unknowns: result.x.length, actualMach: result.mach, targetMach: result.targetMach,
    updates: result.history.length - 1, residual: Math.max(...result.residual.map(Math.abs)),
    boundaryLayers: 4, wakes: 2, acceptedCheckpoints: accepted, meshes, referenceChord: result.referenceChord,
    referenceReynolds: result.referenceReynolds, kernelReynolds: result.kernelReynolds }));
});

test('same-Mach hybrid cache requires no certification, Newton update or LU and rejects changed research controls', () => {
  assert.ok(warm); const parent = serial(warm), before = structuredClone(parent);
  delete parent.flow; delete parent.x; delete parent.residual; delete parent.boundaryLayer;
  const result = solveCoupledStreamtubeMach(parent.sourceCase, { parentResult: parent, stageMaxIterations: 0 });
  assert.equal(result.converged, true); assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(result.continuation.attempts.length, 0); assert.equal(result.machContinuation.certification, undefined);
  assert.deepEqual(serial(result.checkpoint), before.checkpoint);
  assert.throws(() => coupledMachPlan(parent.sourceCase, parent, { epsilonP: 2e-5 }), /unchanged explicit historical hybrid/);
  assert.throws(() => coupledMachPlan(parent.sourceCase, parent, { upwind: {
    mucon: -1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } }), /unchanged explicit historical hybrid/);
});

test('bounded zero-update failure retains accepted Mach and explicitly separates the requested case', () => {
  assert.ok(warm); const parent = serial(warm), caseData = { ...parent.sourceCase, mach: .22 };
  const result = solveCoupledStreamtubeMach(caseData, { parentResult: parent, stageMaxIterations: 0, maxSubdivisions: 0, maxStages: 1 });
  assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
  assert.equal(result.status, 'research-coupled-target-not-reached');
  for (const mach of [result.mach, result.actualMach, result.conditions.mach, result.sourceCase.mach,
    result.checkpoint.restart.input.mach, result.solverInput.mach]) assert.equal(mach, .21);
  assert.equal(result.targetMach, .22); assert.deepEqual(result.requestedCase, caseData);
  assert.equal(result.machContinuation.reachedTarget, false); assert.equal(result.continuation.attempts.length, 1);
  assert.equal(result.continuation.attempts[0].globalLinearSolves, 0);
  assert.deepEqual(serial(result.checkpoint), parent.checkpoint);
});

test('observer cancellation and corrupt full source state never become automatic recovery', () => {
  const parent = source(), caseData = { ...parent.sourceCase, mach: .21 };
  for (const hook of ['onStage', 'onPrepared', 'onCheckpoint']) {
    const cancel = new Error(`cancel ${hook}`);
    assert.throws(() => solveCoupledStreamtubeMach(caseData, { parentResult: parent, [hook]: () => { throw cancel; } }), e => e === cancel);
  }
  let caught = false;
  try { solveCoupledStreamtubeMach(caseData, { parentResult: parent, onStage: () => { throw undefined; } }); }
  catch (e) { caught = true; assert.equal(e, undefined); }
  assert.equal(caught, true);
  const corrupt = structuredClone(parent); corrupt.checkpoint.restart.initialEuler.x[0] += .001;
  let prepared = 0;
  assert.throws(() => solveCoupledStreamtubeMach(caseData, { parentResult: corrupt, onPrepared: () => prepared++ }), /residual does not replay exactly/);
  assert.equal(prepared, 0);
});

test('invalid cold controls stop before mesh creation and report actual baseline Mach without target success', () => {
  const fixture = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const caseData = { flowModel: 'streamtube-grid', quadBoundaryLayers: true, transitionMode: 'automatic',
    elements: fixture.bodies.map(b => ({ points: b.points })), gridIntervals: 1, mach: .4 };
  let meshes = 0, iterations = 0, checkpoints = 0; const stages = [];
  const result = solveCoupledStreamtubeMach(caseData, {
    onStage: s => stages.push(s), onMesh: () => meshes++, onIteration: () => iterations++, onCheckpoint: () => checkpoints++,
  });
  assert.equal(result.status, 'research-coupled-initialization-failed');
  assert.match(result.reason, /Invalid initial streamtube topology controls/);
  assert.equal(result.availableFlow, false); assert.equal(result.converged, false); assert.equal(result.stateConverged, false);
  assert.equal(result.sourceCase.mach, .2); assert.equal(result.mach, .2); assert.equal(result.actualMach, .2);
  assert.equal(result.targetMach, .4); assert.deepEqual(result.requestedCase, caseData);
  assert.deepEqual(result.solverSettings, { maxIterations: 40, eulerMaxIterations: 20, maxStartupAttempts: 2, tolerance: 1e-10 });
  assert.equal(result.machContinuation.route, 'cold-baseline'); assert.equal(result.machContinuation.reachedTarget, false);
  assert.equal(meshes + iterations + checkpoints, 0); assert.equal(stages.length, 1);
  assert.equal(stages[0].actualMach, .2); assert.equal(stages[0].targetMach, .4);
  const low = solveCoupledStreamtubeMach({ ...caseData, mach: .25, transitionMode: 'fixed-trip' });
  assert.equal(low.actualMach, .25); assert.equal(low.targetMach, .25); assert.equal(low.converged, false);
  assert.equal(low.solverSettings.maxIterations, 20); assert.equal(low.machContinuation.route, 'cold');
  let caught = false;
  try { solveCoupledStreamtubeMach(caseData, { onStage: () => { throw undefined; } }); }
  catch (error) { caught = true; assert.equal(error, undefined); }
  assert.equal(caught, true, 'The old cold adapter cannot swallow an observer cancellation.');
});

test('a typed nonzero-bias mismatch permits exactly one cold route only after complete baseline replay', () => {
  const parent = source(), before = structuredClone(parent), caseData = { ...parent.sourceCase, mach: .21 };
  const upwind = { mucon: 1, mcrit: 0, boundary: { kind: 'unfiltered-first-two' } };
  assert.throws(() => certifyCoupledStreamtubeHybrid(parent.checkpoint, { epsilonP: 1e-5, upwind }), error => {
    assert.equal(error.code, 'coupled-hybrid-certification-speed-bias');
    assert.ok(Number.isFinite(error.diagnostics.bias)); assert.notEqual(error.diagnostics.bias, 0);
    return true;
  });
  const stages = []; let meshes = 0, iterations = 0;
  // This fixture's explicit four-interval metadata intentionally fails the
  // public cold topology guard (minimum eight), before any mesh or flow.
  const result = solveCoupledStreamtubeMach(caseData, { parentResult: parent, upwind,
    onStage: event => stages.push(event), onMesh: () => meshes++, onIteration: () => iterations++ });
  assert.deepEqual(stages.map(s => s.stage), ['hybrid-certification', 'coupled-cold-recovery', 'euler']);
  assert.equal(stages[1].actualMach, .105); assert.equal(stages[1].targetMach, .21);
  assert.equal(stages[1].recovery.originalSourceReplayed, true);
  assert.equal(stages[1].recovery.originalSourceConverged, true);
  assert.equal(result.machContinuation.coldRecovery.coldAttempts, 1);
  assert.equal(result.machContinuation.coldBaselineUsed, true);
  assert.equal(result.machContinuation.failureStage, 'cold-initialization');
  assert.equal(result.mach, .2); assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
  assert.deepEqual(serial(result.checkpoint), parent.checkpoint);
  assert.equal(result.linearDiagnostics.solves, 0); assert.equal(meshes + iterations, 0);
  assert.deepEqual(parent, before);
  for (const mutate of [p => { p.checkpoint.continuation.fractions[0][1] = 0; },
    p => { p.checkpoint.continuation.lastRedistributedStagnation[0] = NaN; },
    p => { p.checkpoint.restart.initialEuler.x[0] += .001; }]) {
    const corrupt = structuredClone(parent), seen = []; mutate(corrupt);
    assert.throws(() => solveCoupledStreamtubeMach(caseData, { parentResult: corrupt, upwind,
      onStage: event => seen.push(event.stage) }), /inlet fractions|continuation state|residual does not replay exactly/);
    assert.equal(seen.includes('coupled-cold-recovery'), false, 'Corrupt source cannot authorize a cold replacement.');
  }
});

test('public continuation forwards bounded defaults and explicit overrides while preserving recovery stage and actual Mach', async () => {
  assert.ok(warm); const parent = serial(warm), before = structuredClone(parent), calls = [];
  // Replace only the automatic numerical boundary. Exercise real planning,
  // source replay, normalization and callback adaptation without a new solve.
  const key = '__publicWakePolicyTest';
  globalThis[key] = (targetMach, options) => {
    calls.push({ targetMach, stageMaxIterations:options.stageMaxIterations, maxStages:options.maxStages, maxWakeRecoveries:options.maxWakeRecoveries, blPredictor:options.blPredictor });
    options.onStage({ mach:targetMach });
    const event = { mach:parent.mach, targetMach, stage:'coupled-wake-grid', fraction:0, wakeRecovery:1 };
    options.onStage(event); options.onIteration({ ...event, iteration:0 });
    options.onMesh?.({ ...event, mesh:{ vertices:[{x:1,y:2}],initialization:{} },iteration:{iteration:0} }, { flow:{ marker:'detached' } });
    options.onCheckpoint(options.initialCheckpoint, { ...event, kind:'accepted', reachedTarget:false });
    return { ...structuredClone(parent), converged:false, stateConverged:true };
  };
  const proxy = 'data:text/javascript;base64,'+Buffer.from(`export const solveCoupledStreamtubeAutomatic=globalThis[${JSON.stringify(key)}];`).toString('base64');
  try {
    const code = fs.readFileSync('src/euler/tests/streamtube-coupled-mach-assembly.js','utf8').replace(/from '([^']+)'/g,(_,specifier) =>
      `from '${specifier==='./streamtube-coupled-automatic.js' ? proxy : pathToFileURL(path.resolve('src/euler/tests',specifier)).href}'`);
    const adapter = await import('data:text/javascript;base64,'+Buffer.from(code).toString('base64'));
    for (const controls of [{},{stageMaxIterations:7,maxStages:5,maxWakeRecoveries:0}]) {
      const stages=[], recoveryCallbacks=[];
      const result = adapter.solveCoupledStreamtubeMach({...parent.sourceCase,mach:.22},{parentResult:parent,...controls,
        onStage:event=>stages.push(event),
        onIteration:event=>recoveryCallbacks.push(event),
        onMesh:(mesh,phase,event)=>{assert.equal(mesh.actualMach,.21);recoveryCallbacks.push(event);},
        onCheckpoint:(cp,event)=>{if(event.stage==='coupled-wake-grid') {recoveryCallbacks.push(event);cp.restart.initialBL[0]=-99;}},
      });
      assert.deepEqual(stages.map(s=>s.stage),['coupled-mach','coupled-wake-grid']);
      assert.equal(stages[1].actualMach,.21); assert.equal(stages[1].targetMach,.22);
      assert.equal(recoveryCallbacks.length,3);
      for(const event of recoveryCallbacks) {
        assert.equal(event.stage,'coupled-wake-grid');assert.equal(event.actualMach,.21);assert.equal(event.targetMach,.22);
      }
      assert.equal(result.solverSettings.stageMaxIterations,controls.stageMaxIterations??20);
      assert.equal(result.solverSettings.maxStages,controls.maxStages??48);
      assert.equal(result.solverSettings.maxWakeRecoveries,controls.maxWakeRecoveries??4);
      assert.equal(result.converged,false);assert.equal(result.actualMach,.21);assert.equal(result.targetMach,.22);
    }
    assert.deepEqual(calls,[{targetMach:.22,stageMaxIterations:20,maxStages:48,maxWakeRecoveries:4,blPredictor:'preserve'},
      {targetMach:.22,stageMaxIterations:7,maxStages:5,maxWakeRecoveries:0,blPredictor:'preserve'}]);
    // Selecting the global XFOIL update changes preparation policy only;
    // this accepted source has the same physical equations and state.
    const nativeParent=structuredClone(parent);
    nativeParent.checkpoint.continuation.blUpdate='xfoil';
    const selected=adapter.solveCoupledStreamtubeMach({...nativeParent.sourceCase,mach:.22},{parentResult:nativeParent});
    assert.equal(calls.at(-1).blPredictor,'xfoil-mrchdu');
    assert.equal(selected.solverSettings.blPredictor,'xfoil-mrchdu');
    adapter.solveCoupledStreamtubeMach({...nativeParent.sourceCase,mach:.22},{parentResult:nativeParent,blPredictor:'preserve'});
    assert.equal(calls.at(-1).blPredictor,'preserve');
    assert.deepEqual(parent,before);
  } finally {delete globalThis[key];}
});

test('public warm alpha continuation retains requested Mach, normalization and actual incidence', () => {
  const original = source();
  const parent = solveCoupledStreamtubeMach({ ...original.sourceCase, mach: .21 }, { parentResult: original });
  assert.equal(parent.converged, true, parent.reason);
  const request = { ...parent.sourceCase, alpha: .1 }, events = [];
  const result = solveCoupledStreamtubeMach(request, { parentResult: parent,
    onStage: e => events.push(e) });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.actualAlpha, .1); assert.equal(result.alpha, .1);
  assert.equal(result.checkpoint.restart.input.alpha, .1);
  assert.equal(result.actualMach, .21);
  assert.deepEqual(result.sourceCase, request);
  assert.equal(result.referenceChord, parent.referenceChord);
  assert.equal(result.referenceReynolds, parent.referenceReynolds);
  assert.ok(events.some(e => e.stage === 'coupled-alpha' && e.targetAlpha === .1));
});

test('public alpha live CL/CD/Cm update and equal final coefficients at the actual angle', async () => {
  const { quadCoupledTransonicCoefficients: loads, quadCoupledTransonicCoefficientsFromResult: finalLoads }
    = await import('../src/ui/quad-coupled-transonic-coefficients.js');
  const original = source();
  const parent = solveCoupledStreamtubeMach({ ...original.sourceCase, mach: .21 }, { parentResult: original });
  const frames = [];
  const result = solveCoupledStreamtubeMach({ ...parent.sourceCase, alpha: .1 }, { parentResult: parent,
    onFlow: frame => frames.push({ iteration: frame.iteration.iteration, alpha: frame.checkpoint.restart.input.alpha,
      loads: loads({ ...frame, ...frame.normalization }) }) });
  assert.equal(result.converged, true, result.reason);
  assert.ok(frames.length > 2);
  assert.ok(frames.every(f => f.alpha === .1 && ['cl', 'cd', 'cm'].every(k => Number.isFinite(f.loads[k]))));
  assert.ok(frames.some(f => f.loads.cl !== frames[0].loads.cl));
  const expected = finalLoads(result);
  for (const k of ['cl', 'cd', 'cm']) assert.equal(frames.at(-1).loads[k], expected[k]);
});

test('public warm combined continuation checks stale controls and updates both operating coordinates', () => {
  const original = source();
  const parent = solveCoupledStreamtubeMach({ ...original.sourceCase, mach: .21 }, { parentResult: original });
  const request = { ...parent.sourceCase, mach: .22, alpha: .1 }, frames = [];
  const result = solveCoupledStreamtubeMach(request, { parentResult: parent, onFlow: f => frames.push(f) });
  assert.equal(result.converged, true, result.reason);
  assert.deepEqual(result.sourceCase, request);
  assert.equal(result.actualMach, .22); assert.equal(result.actualAlpha, .1);
  assert.equal(result.machContinuation.reachedTarget, true);
  assert.ok(frames.some(f => f.stage === 'coupled-operating-point' && f.stepMethod === 'combined'));
  assert.equal(result.referenceReynolds, parent.referenceReynolds);
  assert.throws(() => solveCoupledStreamtubeMach({ ...request, reynolds: request.reynolds * 2 },
    { parentResult: parent }), /stale/);
});
