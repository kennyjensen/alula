// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { coupledWakeRecoveryPlan, recoverCoupledWakeCorrespondence } from '../src/euler/tests/streamtube-coupled-wake-recovery.js';
import { correctCoupledWakeCoordinates, continueCoupledWakeCoordinates } from '../src/euler/streamtube-coupled-startup.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('wake-coordinate correction rejects changing-law residual comparisons', () => {
  assert.throws(() => correctCoupledWakeCoordinates({ continuation: {
    dissipationEnhancement: { method: 'mses', targetMcrit: .99, previousDensityChange: .1 },
  } }), /fixed dissipation law/);
});

test('app wake correction requires the specific merit failure and remaining original budget', () => {
  const endpoint = () => ({ converged: false, history: [{ iteration: 4 }],
    mesh: { quality: { valid: true, invalidCells: [] } },
    checkpoint: { version: 1, restart: { input: { wakeGeometry: 'independent-banks' } },
      continuation: { stepAcceptance: 'event-armijo' } },
    lastRejectedStep: { code: 'COUPLED_RESIDUAL_DECREASE', wakeCoordinateRepair: { method: 'paired-wake-trial' } } });
  const mutations = [
    x => { x.converged = true; },
    x => { x.wakeCorrection = { accepted: false }; },
    x => { x.history[0].iteration = 8; },
    x => { x.history = []; },
    x => { x.lastRejectedStep.code = 'KLU_RESIDUAL_LIMIT'; },
    x => { delete x.lastRejectedStep.wakeCoordinateRepair; },
    x => { x.mesh.quality.valid = false; },
    x => { x.mesh.quality.invalidCells = [1]; },
    x => { x.checkpoint.restart.input.wakeGeometry = 'shared'; },
    x => { x.checkpoint.continuation.dissipationEnhancement = { method: 'mses' }; },
    x => { x.checkpoint.continuation.stepAcceptance = 'admissible'; },
  ];
  for (const mutate of mutations) {
    const source = endpoint(); mutate(source); const before = structuredClone(source);
    const result = continueCoupledWakeCoordinates(source, { maxIterations: 8,
      onCorrection: () => assert.fail('An ineligible endpoint must not start a correction.') });
    assert.equal(result, source); assert.deepEqual(source, before);
  }
  for (const maxIterations of [-1, 1.5, NaN, undefined])
    assert.throws(() => continueCoupledWakeCoordinates(endpoint(), { maxIterations }), /budget/);
});

test('app late-wake correction converges the retained RAE64x9 state within one shared budget', t => {
  const cp = JSON.parse(fs.readFileSync(new URL(
    '../docs/solver-reliability/rae-inviscid-audit/fixtures/rae64x9-before-late-wake-merit-rejection.json', import.meta.url)));
  const original = structuredClone(cp);
  const before = solveCoupledStreamtubeIses(undefined, { ...cp.continuation, resume: cp,
    maxIterations: 1, tolerance: 1e-10 });
  assert.equal(before.converged, false);
  assert.equal(before.lastRejectedStep.code, 'COUPLED_RESIDUAL_DECREASE');
  const frozen = structuredClone(before.checkpoint), phases = [], frames = [];
  let latest;
  const result = continueCoupledWakeCoordinates(before, { maxIterations: 8,
    onCorrection: h => phases.push(h.phase),
    onIteration: h => frames.push(h),
    onCheckpoint: value => { latest = value; },
    onMesh: state => {
      assert.ok(latest, 'Only committed flow may be published, with its matching checkpoint.');
      assert.deepEqual(state.nodes, latest.restart.initialEuler.nodes);
      assert.deepEqual(state.coupledFamilies, latest.families);
    } });
  assert.equal(result.converged, true, result.reason);
  assert.deepEqual(phases, ['start', 'trial', 'accepted']);
  assert.ok(result.wakeCorrection.afterSquaredNorm < result.wakeCorrection.beforeSquaredNorm);
  assert.equal(result.wakeCorrection.attemptedCorrectorIterations, 1);
  assert.ok(result.history.at(-1).iteration <= 8);
  assert.equal(result.history.length, result.history.at(-1).iteration + 1);
  assert.deepEqual(frames.map(h => h.iteration), result.history.slice(1).map(h => h.iteration));
  assert.ok(Object.values(result.families).every(r => r <= 1e-10));
  assert.equal(result.mesh.quality.valid, true);
  assert.ok(result.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert.deepEqual(result.checkpoint.restart.input, cp.restart.input);
  for (const key of ['reynolds', 'ncrit', 'tripFractions', 'transitionMode'])
    assert.deepEqual(result.checkpoint.restart.options[key], cp.restart.options[key]);
  assert.deepEqual(before.checkpoint, frozen); assert.deepEqual(cp, original);
  t.diagnostic(JSON.stringify({ families: result.families, updates: result.history.at(-1).iteration,
    linearSolves: result.linearDiagnostics.solves }));
});

for (const [label, name, fraction] of [
  ['RAE64x7', 'rae64x7-localized-before-wake-contact', .25],
  ['RAE64x9', 'rae64x9-surface-hk-before-wake-contact', .125],
]) test(`${label} damped wake correspondence reduces merit without changing BL or gas unknowns`, () => {
  const cp = JSON.parse(fs.readFileSync(new URL(
    `../docs/solver-reliability/rae-inviscid-audit/fixtures/${name}.json`, import.meta.url))), before = structuredClone(cp);
  const recovered = correctCoupledWakeCoordinates(cp, { maxCoordinateBacktracks: 12, maxCorrectorIterations: 0 });
  assert.equal(recovered.report.accepted, true, recovered.report.reason);
  assert.equal(recovered.report.coordinateFraction, fraction);
  assert.equal(recovered.report.correctorIterations, 0);
  assert.ok(recovered.report.afterSquaredNorm < recovered.report.beforeSquaredNorm);
  const f = cp.restart, next = recovered.checkpoint.restart;
  assert.deepEqual(next.input, f.input);
  assert.deepEqual(Array.from(next.initialBL), f.initialBL);
  assert.deepEqual(next.options.transitionState, f.options.transitionState);
  assert.deepEqual(recovered.checkpoint.continuation, cp.continuation);
  const source = createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: f.initialBL });
  const candidate = createCoupledStreamtubeBody(next.input, { ...next.options,
    initialEuler: next.initialEuler, initialBL: next.initialBL });
  const value = candidate.admissibleValue(candidate.initial, { requireConvex: true });
  assert.ok(value);
  assert.deepEqual(value.families, recovered.checkpoint.families);
  assert.deepEqual(candidate.initial.slice(0, candidate.euler.layout.densityCount),
    source.initial.slice(0, source.euler.layout.densityCount));
  assert.deepEqual(cp, before);
});

for (const [label, fixture, rejectedBudget, minimumUpdates, maximumUpdates, continuationBudget] of [
  ['RAE64x7', 'rae64x7-localized-before-wake-contact', 1, 2, 5, 0],
  ['RAE64x9 late', 'rae64x9-before-late-wake-merit-rejection', 0, 1, 1, 8],
]) test(`bounded ${label} wake correction commits only after reducing the original coupled merit`, t => {
  const cp = JSON.parse(fs.readFileSync(new URL(
    `../docs/solver-reliability/rae-inviscid-audit/fixtures/${fixture}.json`, import.meta.url)));
  const before = structuredClone(cp), f = cp.restart;
  const source = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const original = source.evaluate(source.initial);
  const rejected = correctCoupledWakeCoordinates(cp, { maxCorrectorIterations: rejectedBudget });
  assert.equal(rejected.report.accepted, false);
  assert.equal(rejected.checkpoint, undefined, 'An unconverged correction must not leak an adoptable state.');
  assert.deepEqual(cp, before);
  const recovered = correctCoupledWakeCoordinates(cp);
  assert.equal(recovered.report.accepted, true, recovered.report.reason);
  assert.ok(recovered.report.uncommittedGuess.squaredNorm > recovered.report.beforeSquaredNorm);
  assert.ok(recovered.report.afterSquaredNorm < recovered.report.beforeSquaredNorm);
  assert.ok(recovered.report.correctorIterations >= minimumUpdates
    && recovered.report.correctorIterations <= maximumUpdates,
  'The correction must honor the original merit gate and bounded work allowance.');
  assert.equal(recovered.report.equationsChanged, false);
  assert.equal(recovered.report.physicalConditionsChanged, false);
  const next = recovered.checkpoint.restart;
  assert.deepEqual(next.input, f.input);
  for (const key of ['reynolds', 'ncrit', 'tripFractions', 'transitionMode', 'edgeMatching', 'blThermodynamics'])
    assert.deepEqual(next.options[key], f.options[key]);
  for (const key of ['fractions', 'iterationGeometry', 'stepAcceptance', 'stagnationLimiter',
    'blUpdate', 'linearOrdering', 'hkProjectionRecovery', 'dissipationEnhancement'])
    assert.deepEqual(recovered.checkpoint.continuation[key], cp.continuation[key]);
  const candidate = createCoupledStreamtubeBody(next.input, { ...next.options,
    initialEuler: next.initialEuler, initialBL: next.initialBL });
  const value = candidate.admissibleValue(candidate.initial, { requireConvex: true });
  assert.ok(value);
  assert.deepEqual(value.families, recovered.checkpoint.families);
  assert.deepEqual(source.evaluate(source.initial).residual, original.residual);
  assert.deepEqual(cp, before);
  if (continuationBudget) {
    const result = solveCoupledStreamtubeIses(undefined, { resume: recovered.checkpoint,
      ...recovered.checkpoint.continuation, maxIterations: continuationBudget, tolerance: 1e-10 });
    assert.equal(result.converged, true, result.reason);
    assert.ok(Object.values(result.families).every(r => r <= 1e-10));
    assert.equal(result.mesh.quality.valid, true);
    assert.ok(result.linearDiagnostics.maxRelativeResidual <= 1e-10);
    assert.ok(result.history.slice(-3).every(h => h.step === 1 && h.backtracks === 0));
    assert.deepEqual(result.checkpoint.restart.input, cp.restart.input);
    for (const key of ['reynolds', 'ncrit', 'tripFractions', 'transitionMode'])
      assert.deepEqual(result.checkpoint.restart.options[key], cp.restart.options[key]);
    assert.deepEqual(result.checkpoint.continuation.fractions, cp.continuation.fractions);
    assert.equal(result.checkpoint.continuation.wakeCoordinateRecovery, undefined);
    assert.deepEqual(cp, before);
    t.diagnostic(JSON.stringify({ families: result.families, updates: result.history.length - 1 }));
  }
  t.diagnostic(JSON.stringify({ correctorIterations: recovered.report.correctorIterations,
    beforeSquaredNorm: recovered.report.beforeSquaredNorm,
    guessSquaredNorm: recovered.report.uncommittedGuess.squaredNorm,
    afterSquaredNorm: recovered.report.afterSquaredNorm }));
});

test('wake recovery is limited to finite failed targets with folded TE-to-first-wake cells',()=>{
  // Five intervals and two tubes in each passage. Body TE is station3:
  // lower wake-bank cell7, upper wake-bank cell16; adjacent cell17 is interior.
  const checkpoint={restart:{input:{wakeGeometry:'independent-banks',outerLower:Array(6),bodies:[{trailingIndex:3}]}}};
  const trial={residual:Float64Array.of(1e-12,2e-12,1e-12),families:{euler:1e-12,boundaryLayer:2e-12,edgeMatching:1e-12},
    mesh:{quality:{valid:false,invalidCells:[16]}},flow:{nodes:Array.from({length:2},()=>Array.from({length:6},()=>Array(3)))}};
  const before=structuredClone(trial);
  assert.deepEqual(coupledWakeRecoveryPlan(trial,checkpoint,1e-10)?.bodies,[0]);
  assert.equal(coupledWakeRecoveryPlan(trial,checkpoint,1e-10).residualConverged,true);
  assert.deepEqual(trial,before);
  trial.mesh.quality.invalidCells=[7,16];assert.ok(coupledWakeRecoveryPlan(trial,checkpoint,1e-10));
  for(const ids of [[17],[6],[16,17],[]]){
    trial.mesh.quality.invalidCells=ids;assert.equal(coupledWakeRecoveryPlan(trial,checkpoint,1e-10),null);
  }
  trial.mesh.quality.invalidCells=[16];
  trial.families.boundaryLayer=1.5; trial.residual[1]=-1.5;
  const finiteFailure=coupledWakeRecoveryPlan(trial,checkpoint,1e-10);
  assert.equal(finiteFailure.residualConverged,false);
  assert.deepEqual(finiteFailure.targetFamilies,trial.families);
  finiteFailure.targetFamilies.boundaryLayer=99;
  assert.equal(trial.families.boundaryLayer,1.5);
  for(const value of [NaN,-1,Infinity]){
    trial.families.boundaryLayer=value;assert.equal(coupledWakeRecoveryPlan(trial,checkpoint,1e-10),null);
  }
  trial.families.boundaryLayer=1.5;
  for(const residual of [[],undefined,[NaN],[Infinity],new DataView(new ArrayBuffer(8))]){
    trial.residual=residual;assert.equal(coupledWakeRecoveryPlan(trial,checkpoint,1e-10),null);
  }
  trial.residual=[-1.5];assert.equal(coupledWakeRecoveryPlan(trial,checkpoint,NaN),null);
  trial.mesh.quality.invalidCells=[16,17];assert.equal(coupledWakeRecoveryPlan(trial,checkpoint,1e-10),null);
  trial.mesh.quality.invalidCells=[16];trial.mesh.quality.valid=true;
  assert.equal(coupledWakeRecoveryPlan(trial,checkpoint,1e-10),null);
});

test('one-time wake alignment reconverges an independent-bank coupled state and preserves its source',()=>{
  const input={...intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2,mach:.2}),
    wakeGeometry:'independent-banks',wakeOutlet:'banks',streamwiseMode:'hybrid',hybrid:{epsilonP:1e-5},
    upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}};
  const original=solveCoupledStreamtubeIses(input,{edgeMatching:'section-velocity',blThermodynamics:'historical-common-isentrope',
    transitionMode:'automatic',blUpdate:'xfoil',maxIterations:12,tolerance:1e-10,stepAcceptance:'admissible'});
  assert.equal(original.converged,true,original.reason);
  const cp=structuredClone(original.checkpoint),before=structuredClone(cp);
  let iterations=0;
  const recovered=recoverCoupledWakeCorrespondence(cp,{onIteration:()=>{iterations++;}});
  assert.equal(recovered.result.converged,true,recovered.result.reason);
  assert.equal(recovered.result.conditions.mach,.2);
  assert.equal(recovered.result.mesh.quality.valid,true);
  assert.ok(Math.max(...Object.values(recovered.result.families))<=1e-10);
  assert.equal(recovered.result.boundaryLayer.surfaces.length,4);
  assert.equal(recovered.result.boundaryLayer.wakes.length,2);
  assert.equal(recovered.result.checkpoint.continuation.blUpdate,'xfoil');
  assert.equal(recovered.diagnostics.equationsChanged,false);
  assert.equal(recovered.diagnostics.tangentialGapConstrainedDuringNewton,false);
  assert.ok(iterations>0);assert.deepEqual(cp,before);
  assert.deepEqual(recovered.result.checkpoint.restart.input,cp.restart.input);
  assert.deepEqual(recovered.result.checkpoint.restart.options.tripFractions,cp.restart.options.tripFractions);
});
