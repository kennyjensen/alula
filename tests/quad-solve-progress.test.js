import test from 'node:test';
import assert from 'node:assert/strict';
import { createQuadSolveProgress, quadStartingGridLabel } from '../src/ui/quad-solve-progress.js';

test('targeted Hk recovery is visible during direct solving and persists until another method starts', () => {
  const progress = createQuadSolveProgress({ ncrit: 4, mach: .74 });
  progress.update({ stage: 'coupled', startupAttempt: 1, startupStrategy: 'direct-requested-conditions',
    shearCoordinate: 'logarithmic', hkFloorLinearization: 'exact' }, { stageChange: true });
  const event = progress.update({ stage: 'coupled', iteration: 10, step: 0, accepted: false,
    hkProjectionRecovery: { from: 'exact', to: 'native', equationsChanged: false }, hkFloorLinearization: 'native' });
  assert.match(event.current, /Step 2\/2.*Hk thickness-projection recovery.*native Hk derivative/);
  assert.equal(event.history.length, 2);
  const next = progress.update({ stage: 'coupled', iteration: 11, step: 1 });
  assert.equal(next.current, event.current);
  assert.equal(next.history.length, 2);
  const fresh = progress.update({ stage: 'euler' }, { stageChange: true });
  assert.doesNotMatch(fresh.current, /Hk thickness-projection recovery|native Hk/);
});

test('ordinary, lower-Ncrit and logarithmic methods remain distinct across metadata-light updates', () => {
  const progress = createQuadSolveProgress({ ncrit: 9, mach: .2 });
  const stage = event => progress.update(event, { stageChange: true });
  stage({ stage: 'euler', startupAttempt: 0 });
  stage({ stage: 'boundary-layer-initialization', startupAttempt: 1, shearCoordinate: 'linear', hkFloorLinearization: 'exact' });
  let view = stage({ stage: 'coupled', startupAttempt: 1 });
  assert.match(view.current, /Ordinary Euler\/BL startup.*linear shear.*exact Hk/);
  stage({ stage: 'boundary-layer-initialization', startupAttempt: 2, actualNcrit: 4, targetNcrit: 9,
    shearCoordinate: 'linear', hkFloorLinearization: 'native' });
  view = stage({ stage: 'coupled', startupAttempt: 2 });
  assert.match(view.current, /Lower-Ncrit recovery.*linear shear.*native Hk.*Ncrit 4 → 9/);
  assert.match(progress.update({ stage: 'coupled', iteration: 40 }).current, /Lower-Ncrit recovery.*linear shear/);
  view = stage({ stage: 'coupled', startupAttempt: 2, shearRecovery: { shearCoordinate: 'logarithmic' } });
  assert.match(view.current, /Logarithmic shear recovery.*native Hk/);
  const count = view.history.length;
  view = progress.update({ stage: 'coupled', iteration: 41 });
  assert.match(view.current, /Logarithmic shear recovery.*native Hk/);
  assert.equal(view.history.length, count);
  assert.equal(view.history.filter(label => label.startsWith('Logarithmic')).length, 1);
  assert.ok(view.history.some(label => label.startsWith('Ordinary')));
  assert.ok(view.history.some(label => label.startsWith('Lower-Ncrit')));
});

test('continued condition/refinement phases preserve coordinate but replace the recovery label', () => {
  const progress = createQuadSolveProgress({ ncrit: 9, mach: .3 }, { checkpoint: {
    restart: { input: { mach: .2 }, options: { ncrit: 4, hkFloorLinearization: 'native' } },
    continuation: { shearCoordinate: 'logarithmic' },
  } });
  const stage = event => progress.update(event, { stageChange: true });
  let view = stage({ stage: 'coupled-refinement', startupAttempt: 1 });
  assert.match(view.current, /Grid refinement.*logarithmic shear.*native Hk/);
  view = stage({ stage: 'coupled', startupAttempt: 1, refinement: true });
  assert.match(view.current, /Grid refinement/);
  view = stage({ stage: 'coupled-ncrit-startup', startupAttempt: 1, actualNcrit: 5, targetNcrit: 9 });
  assert.match(view.current, /Ncrit continuation.*logarithmic shear.*Ncrit 5 → 9/);
  view = stage({ stage: 'coupled-mach', startupAttempt: 1, mach: .25, targetMach: .3 });
  assert.match(view.current, /Mach continuation.*logarithmic shear.*Mach 0.250 → 0.300/);
  assert.equal(progress.update({ stage: 'coupled-mach', iteration: 2 }).current, view.current);
});

test('fresh initialization clears prior log/native metadata and a new request clears history', () => {
  const progress = createQuadSolveProgress({ ncrit: 9, mach: .2 });
  progress.update({ stage: 'coupled', startupAttempt: 2, actualNcrit: 4,
    shearRecovery: { shearCoordinate: 'logarithmic' }, hkFloorLinearization: 'native' });
  let view = progress.update({ stage: 'boundary-layer-initialization', startupAttempt: 3 }, { stageChange: true });
  assert.doesNotMatch(view.current, /logarithmic|native/);
  view = progress.update({ stage: 'coupled', startupAttempt: 3, shearCoordinate: 'linear', hkFloorLinearization: 'exact' });
  assert.match(view.current, /Lower-Ncrit recovery.*linear shear.*exact Hk/);
  assert.doesNotMatch(view.current, /native/);
  const fresh = createQuadSolveProgress({ ncrit: 7, mach: .15 }).update({ stage: 'euler' }, { stageChange: true });
  assert.equal(fresh.history.length, 1);
  assert.doesNotMatch(fresh.current, /recovery|native|logarithmic/);
});

test('coupled iteration zero is described as initialization before Newton', () => {
  assert.equal(quadStartingGridLabel({ iteration: { stage: 'coupled', iteration: 0 } }, 'initial'),
    'Starting Euler/BL grid · before the first Newton update');
  assert.equal(quadStartingGridLabel({}, 'smoothing'), 'SLOR initialization · Euler has not started');
  assert.equal(quadStartingGridLabel({ initialization: { flowSolved: true } }, 'initial'),
    'Euler converged on the starting grid · 0 updates');
});

test('selected BL initialization is described from its receipt and retained in method history', () => {
  const progress = createQuadSolveProgress({ ncrit: 9, mach: .2 });
  progress.update({ stage: 'boundary-layer-initialization', startupAttempt: 1 }, { stageChange: true });
  let view = progress.update({ stage: 'boundary-layer-initialization', startupAttempt: 1,
    boundaryLayerInitialization: { method: 'MRCHUE surfaces and ISET-style wake guess', wakeInitialization: 'iset-linear-shape', thicknessFactor: 1 },
    shearCoordinate: 'linear', hkFloorLinearization: 'exact' }, { stageChange: true });
  assert.match(view.current, /ISET wake guess/);
  assert.equal(view.history.length, 1, 'The selected initializer fills in its existing stage receipt.');
  view = progress.update({ stage: 'coupled', startupAttempt: 1 }, { stageChange: true });
  assert.match(view.current, /Ordinary Euler\/BL startup/);
  assert.match(view.history[0], /ISET wake guess/);
  view = progress.update({ stage: 'boundary-layer-initialization', startupAttempt: 2,
    boundaryLayerInitialization: { method: 'uniform initial thickness backtracking', thicknessFactor: .5 } }, { stageChange: true });
  assert.doesNotMatch(view.current, /ISET/);
  assert.match(view.current, /thickness × 0.5/);
});

test('target-Mach grid strategy and temporary MCRIT are visible without duplicating iteration methods', () => {
  const p = createQuadSolveProgress({ mach: .74, ncrit: 4 });
  p.update({ stage: 'coupled-mach', mach: .72, gridLevel: 32, gridStrategy: 'target-mach-before-refinement' }, { stageChange: true });
  let view = p.update({ dissipation: { mcrit: .76, targetMcrit: .99 } });
  assert.match(view.current, /MSES dissipation.*0.760 → 0.990/);
  assert.match(view.current, /coarse-to-fine continuation.*grid 32/);
  view = p.update({ dissipation: { mcrit: .99, targetMcrit: .99 } });
  assert.equal(view.history.length, 1);
  assert.match(view.current, /0.990 → 0.990/);
  view = p.update({ stage: 'euler' }, { stageChange: true });
  assert.doesNotMatch(view.current, /MSES dissipation|coarse-to-fine continuation/);
});

test('automatic damping failure is visibly distinguished from a solver failure', () => {
  const progress = createQuadSolveProgress({ mach: .2, ncrit: 4 });
  const view = progress.update({ stage: 'coupled', iteration: 7, progress: { action: 'ordinary-newton',
    recoveryOutcome: { kind: 'ordinary-policy-restored', reason: 'Automatic merit search exhausted' } } });
  assert.match(view.current, /ordinary Newton resumed after damping trial/);
  assert.doesNotMatch(view.current, /stalled stage stopped/);
});

test('Mach seed is explicitly labelled as preparation for alpha continuation and clears on the alpha walk', () => {
 const p=createQuadSolveProgress({mach:.74,alpha:2.68});
 assert.match(p.update({stage:'coupled-mach',alphaSeed:true,actualAlpha:.84,targetAlpha:2.68,mach:.6,targetMach:.74}, {stageChange:true}).current,
  /lower-incidence seed for alpha continuation/);
 assert.doesNotMatch(p.update({stage:'coupled-operating-point',actualAlpha:.94,targetAlpha:2.68,mach:.74,targetMach:.74,stepMethod:'alpha'}, {stageChange:true}).current,
  /lower-incidence seed/);
});

test('direct solve labels distinguish the two requested-condition stages', () => {
  const progress = createQuadSolveProgress({ mach: .74, alpha: 2.68, ncrit: 4 });
  const common = { startupStrategy: 'direct-requested-conditions', mach: .74, targetMach: .74,
    actualAlpha: 2.68, targetAlpha: 2.68, gridLevel: 128 };
  const first = progress.update({ ...common, stage: 'euler' }, { stageChange: true });
  assert.match(first.current, /Step 1\/2: inviscid at requested conditions/);
  progress.update({ ...common, stage: 'boundary-layer-initialization' }, { stageChange: true });
  const last = progress.update({ ...common, stage: 'coupled' }, { stageChange: true });
  assert.match(last.current, /Step 2\/2: viscous at requested conditions/);
  assert.doesNotMatch(last.current, /continuation|recovery/);
});

test('wake correction progress identifies tentative work and clears when ordinary Newton resumes', () => {
  const progress = createQuadSolveProgress({ mach: .74, alpha: 2.68, ncrit: 4 });
  assert.match(progress.update({ stage: 'coupled', wakeCorrection: { phase: 'start' } }, { stageChange: true }).current,
    /checking wake-grid correction/);
  assert.match(progress.update({ stage: 'coupled', wakeCorrection: { phase: 'trial', iteration: 12 } }, { stageChange: true }).current,
    /checking wake-grid correction · trial 12/);
  assert.match(progress.update({ stage: 'coupled', wakeCorrection: { phase: 'rejected' } }, { stageChange: true }).current,
    /rejected · retained previous flow/);
  assert.doesNotMatch(progress.update({ stage: 'coupled', iteration: 13 }).current, /wake-grid/);
});


test('resolution recovery reports adaptive increases and clears when the previous flow is retained', () => {
  const p = createQuadSolveProgress({ mach: .74, ncrit: 9 });
  const event = { stage: 'euler', resolutionRecovery: { requestedGridIntervals: 8, finalNominalGridIntervals: 16 } };
  assert.match(p.update(event, { stageChange: true }).current, /adaptive refinement 8 → 16/);
  assert.match(p.update({ stage: 'coupled', iteration: 1, resolutionRecovery: 'startup' }).current, /adaptive refinement 8 → 16/);
  assert.doesNotMatch(p.update({ stage: 'coupled', retained: true }).current, /adaptive refinement/);
  assert.match(p.update({ ...event, resolutionRecovery: { requestedGridIntervals: 32, finalNominalGridIntervals: 32 } },
    { stageChange: true }).current, /coarse-to-fine grid recovery/);
  assert.doesNotMatch(p.update({ stage: 'euler' }, { stageChange: true }).current, /grid recovery/);
});
