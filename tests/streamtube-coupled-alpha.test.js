// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { solveCoupledStreamtubeAlpha } from '../src/euler/tests/streamtube-coupled-alpha.js';
import { initializeCoupledStreamtubeFromFlow } from '../src/euler/tests/streamtube-coupled-flow-restart.js';
import { createQuadSolveProgress } from '../src/ui/quad-solve-progress.js';
let cp;
function source() {
  if (!cp) {
    const r = solveCoupledStreamtubeIses({ ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 }),
      alpha: 0, streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
      upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } }, {
      edgeMatching: 'section-velocity', blThermodynamics: 'historical-common-isentrope', transitionMode: 'automatic',
      maxIterations: 12, tolerance: 1e-10, stepAcceptance: 'admissible', iterationRecovery: false });
    assert.equal(r.converged, true, r.reason); cp = r.checkpoint;
  }
  return structuredClone(cp);
}
test('alpha transfer preserves gas, geometry, BL thickness and source; changes incidence residual', () => {
  const checkpoint = source(), before = structuredClone(checkpoint);
  const r = initializeCoupledStreamtubeFromFlow(.2, checkpoint, { targetAlpha: .1 });
  assert.equal(r.checkpoint.restart.input.alpha, .1);
  assert.deepEqual(r.checkpoint.restart.initialEuler, checkpoint.restart.initialEuler);
  assert.deepEqual(r.checkpoint.restart.initialBL, checkpoint.restart.initialBL);
  assert.notDeepEqual(r.checkpoint.families, checkpoint.families);
  assert.deepEqual(checkpoint, before);
});
test('alpha continuation converges in both directions with bounded accepted steps', () => {
  for (const target of [.2, -.2]) {
    const checkpoint = source(), accepted = [];
    const r = solveCoupledStreamtubeAlpha(target, { initialCheckpoint: checkpoint,
      onCheckpoint: (cp, e) => accepted.push(e.actualAlpha) });
    assert.equal(r.converged, true, JSON.stringify(r.alphaContinuation));
    assert.equal(r.actualAlpha, target); assert.equal(r.checkpoint.restart.input.alpha, target);
    assert.equal(r.checkpoint.restart.input.mach, .2);
    assert.ok(accepted.length >= 2);
    let previous = 0; for (const alpha of accepted) { assert.ok(Math.abs(alpha - previous) <= .1 + 1e-14); previous = alpha; }
  }
});
test('failed alpha stages retain the exact last root and cancellation propagates', () => {
  const checkpoint = source();
  const r = solveCoupledStreamtubeAlpha(.2, { initialCheckpoint: checkpoint, stageMaxIterations: 0, maxSubdivisions: 1 });
  assert.equal(r.converged, false); assert.equal(r.actualAlpha, 0); assert.equal(r.stateConverged, true);
  assert.deepEqual(r.checkpoint, checkpoint);
  assert.deepEqual(r.alphaContinuation.attempts.map(a => a.alpha), [.1, .05]);
  const cancel = new Error('cancel');
  assert.throws(() => solveCoupledStreamtubeAlpha(.2, { initialCheckpoint: checkpoint, onIteration: () => { throw cancel; } }), e => e === cancel);
});
test('progress shows alpha continuation and actual versus requested incidence', () => {
  const p = createQuadSolveProgress({ alpha: 2.68, mach: .74 });
  assert.match(p.update({ stage: 'coupled-alpha', actualAlpha: .1, targetAlpha: 2.68, mach: .74 }, { stageChange: true }).current,
    /Alpha continuation.*α 0.10° → 2.68°/);
});

test('minimum useful alpha step bounds retries without relaxing residual acceptance', () => {
  const checkpoint = source();
  const r = solveCoupledStreamtubeAlpha(.2, { initialCheckpoint: checkpoint, stageMaxIterations: 0,
    maxSubdivisions: 20, maxTransitionRecoveries: 0 });
  assert.equal(r.alphaContinuation.stopReason, 'minimum alpha step');
  assert.equal(r.alphaContinuation.attempts.length, 10);
  assert.equal(r.converged, false); assert.equal(r.actualAlpha, 0);
  assert.deepEqual(r.checkpoint, checkpoint);
  assert.ok(r.alphaContinuation.attempts.every(a => a.increment >= 1e-4));
  assert.ok(r.alphaContinuation.lastAttempt.families);
});

test('transition refinement is bounded and failed refinement cannot replace an accepted alpha root', async () => {
  const fs = await import('node:fs');
  const url = new URL('../src/euler/tests/streamtube-coupled-alpha.js', import.meta.url);
  let recoveries = 0;
  const key = `alpha-recovery-${Math.random()}`;
  globalThis[key] = {
    transitionRecoveryPlan: () => ({ reason: 'repeated-adjacent-natural-transition', parentNx: 10, refinedNx: 20 }),
    recoverCoupledTransition: (trial, options) => {
      recoveries++; options.onStage({ stage: 'transition-refinement' });
      options.onFlow?.({ checkpoint: trial.checkpoint, flow: {}, bl: { refined: true },
        normalization: options.normalization, iteration: { iteration: 1 } });
      return { ...trial, converged: false, reason: 'controlled refinement failure' };
    },
  };
  const code = fs.readFileSync(url, 'utf8')
    .replace("import { transitionRecoveryPlan, recoverCoupledTransition } from '../streamtube-transition-recovery.js';",
      `const { transitionRecoveryPlan, recoverCoupledTransition } = globalThis[${JSON.stringify(key)}];`)
    .replace(/from '(\.[^']+)'/g, (_, relative) => `from '${new URL(relative, url).href}'`);
  try {
    const { solveCoupledStreamtubeAlpha: solve } = await import(`data:text/javascript;base64,${Buffer.from(code).toString('base64')}`);
    const checkpoint = source(), stages = [], meshes = [], frames = [];
    const r = solve(.2, { initialCheckpoint: checkpoint, stageMaxIterations: 0, maxSubdivisions: 3,
      onStage: e => stages.push(e), onMesh: (mesh, phase) => meshes.push(phase),
      onFlow: frame => frames.push(frame) });
    assert.equal(recoveries, 2); assert.equal(r.alphaContinuation.transitionRecoveries.length, 2);
    assert.ok(r.alphaContinuation.transitionRecoveries.every(r => !r.accepted));
    assert.deepEqual(r.checkpoint, checkpoint); assert.equal(r.actualAlpha, 0);
    assert.equal(stages.filter(e => e.transitionRecovery).length, 2);
    assert.equal(meshes.filter(p => p === 'retained').length, 2);
    assert.equal(frames.filter(f => f.transitionRecovery && f.bl.refined).length, 2);
    assert.ok(frames.filter(f => f.transitionRecovery).every(f => f.actualAlpha === f.checkpoint.restart.input.alpha));
    const cancel = new Error('cancel recovery');
    assert.throws(() => solve(.2, { initialCheckpoint: checkpoint, stageMaxIterations: 0,
      onStage: e => { if (e.transitionRecovery) throw cancel; } }), e => e === cancel);
  } finally { delete globalThis[key]; }
});

test('alpha iterations publish detached matching flow/checkpoint frames even without mesh observers', () => {
  const initialCheckpoint = source(), frames = [];
  const plain = v => JSON.parse(JSON.stringify(v));
  const baseline = solveCoupledStreamtubeAlpha(.1, { initialCheckpoint });
  const result = solveCoupledStreamtubeAlpha(.1, { initialCheckpoint,
    normalization: { solverLength: 1, referenceChord: 1 },
    onFlow: frame => { frames.push(structuredClone(frame)); frame.checkpoint.restart.initialBL[0] = -99; } });
  assert.equal(result.converged, true);
  assert.deepEqual(plain(result.checkpoint), plain(baseline.checkpoint));
  assert.ok(frames.length > 2);
  assert.equal(frames[0].iteration.iteration, 0);
  for (const frame of frames) {
    assert.equal(frame.checkpoint.restart.input.alpha, .1);
    assert.equal(frame.actualAlpha, .1); assert.equal(frame.stage, 'coupled-alpha');
    assert.deepEqual(plain(frame.flow.nodes), plain(frame.checkpoint.restart.initialEuler.nodes));
    assert.equal(frame.bl.stations.length * 4, frame.checkpoint.restart.initialBL.length);
    assert.equal(frame.flow.sections.at(-1).length, frame.checkpoint.restart.input.weights.length);
  }
  assert.deepEqual(plain(frames.at(-1).checkpoint), plain(result.checkpoint));
  const failedFrames = [];
  const failed = solveCoupledStreamtubeAlpha(.1, { initialCheckpoint, stageMaxIterations: 0, maxSubdivisions: 0,
    onFlow: frame => failedFrames.push(frame) });
  assert.equal(failed.converged, false);
  assert.equal(failedFrames.at(-1).retained, true);
  assert.equal(failedFrames.at(-1).actualAlpha, 0);
  assert.deepEqual(plain(failedFrames.at(-1).checkpoint), plain(initialCheckpoint));
  const cancel = new Error('cancel flow');
  assert.throws(() => solveCoupledStreamtubeAlpha(.1, { initialCheckpoint, onFlow: () => { throw cancel; } }), e => e === cancel);
});

test('combined Mach and alpha steps converge with exact target and unchanged source', () => {
  const checkpoint = source(), before = structuredClone(checkpoint), events = [], frames = [];
  const r = solveCoupledStreamtubeAlpha(.2, { initialCheckpoint: checkpoint, targetMach: .22,
    maxMachStep: .01, onStage: e => events.push(e), onFlow: e => frames.push(e) });
  assert.equal(r.converged, true, JSON.stringify(r.operatingPointContinuation));
  assert.equal(r.checkpoint.restart.input.mach, .22);
  assert.equal(r.checkpoint.restart.input.alpha, .2);
  assert.ok(events.some(e => e.stepMethod === 'combined'));
  assert.ok(frames.length > 0);
  for (const f of frames) {
    assert.equal(f.checkpoint.restart.input.mach, f.mach);
    assert.equal(f.checkpoint.restart.input.alpha, f.actualAlpha);
  }
  assert.deepEqual(checkpoint, before);
  assert.ok(Object.values(r.families).every(v => v <= 1e-10));
});

test('combined rejected directions restart from the same root and bound retries', () => {
  const checkpoint = source();
  const r = solveCoupledStreamtubeAlpha(.2, { initialCheckpoint: checkpoint, targetMach: .22,
    stageMaxIterations: 0, maxSubdivisions: 0 });
  assert.equal(r.converged, false);
  assert.deepEqual(r.checkpoint, checkpoint);
  assert.deepEqual(r.operatingPointContinuation.attempts.map(a => a.stepMethod), ['combined', 'alpha', 'mach']);
  const p = createQuadSolveProgress({ alpha: .2, mach: .22 });
  assert.match(p.update({ stage: 'coupled-operating-point', actualAlpha: .1, targetAlpha: .2,
    mach: .21, targetMach: .22, stepMethod: 'combined' }, { stageChange: true }).current,
    /Mach–alpha continuation.*α 0.10°.*Mach 0.210.*combined/);
});
