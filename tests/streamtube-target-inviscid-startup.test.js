// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createQuadSolveProgress } from '../src/ui/quad-solve-progress.js';

async function adapter(solve, alphaSolve) {
  const url = new URL('../src/euler/tests/streamtube-coupled-mach-assembly.js', import.meta.url);
  const key = `target-start-${Math.random()}`;
  globalThis[key] = solve;
  let code = fs.readFileSync(url, 'utf8').replace(
    'coupledCheckpointHkPolicy, coupledAssemblyConditions, solveCoupledStreamtubeAssembly',
    'coupledCheckpointHkPolicy, coupledAssemblyConditions');
  if (alphaSolve) {
    globalThis[`${key}-alpha`] = alphaSolve;
    code = code.replace("import { solveCoupledStreamtubeAlpha } from './streamtube-coupled-alpha.js';",
      `const solveCoupledStreamtubeAlpha = globalThis[${JSON.stringify(`${key}-alpha`)}];`);
  }
  code = `const solveCoupledStreamtubeAssembly = globalThis[${JSON.stringify(key)}];\n` + code;
  code = code.replace(/from '(\.[^']+)'/g, (_, relative) => `from '${new URL(relative, url).href}'`);
  try { const module = await import(`data:text/javascript;base64,${Buffer.from(code).toString('base64')}`);
    return (c, o) => module.solveCoupledStreamtubeMach(c, { alphaContinuation: false, targetGridSequencing: false, ...o }); }
  finally { delete globalThis[key]; delete globalThis[`${key}-alpha`]; }
}
const input = { flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .74, alpha: 2.68,
  elements: [{ points: [{x:1,y:0},{x:0,y:0},{x:1,y:0}] }], reynolds: 2700000, ncrit: 4, transitionMode: 'automatic', gridIntervals: 32, gridTubes: 7 };
const result = (c, converged) => ({ converged, mach: c.mach, conditions: { mach: c.mach },
  mesh: { quality: { valid: true } }, reason: converged ? 'converged' : 'iteration limit',
  initialization: { euler: { iterations: 20 } } });

test('cold transonic target performs one inviscid/viscous attempt without changing operating conditions', async () => {
  const calls = [], events = [];
  const solve = await adapter((c, o) => {
    calls.push({ c, o });
    o.onStage({ stage: 'euler' }); o.onStage({ stage: 'boundary-layer-initialization' });
    o.onIteration({ stage: 'coupled', iteration: 1 });
    return result(c, true);
  });
  const r = solve(input, { onStage: e => events.push(e), onIteration: e => events.push(e) });
  assert.equal(calls.length, 1); assert.deepEqual(calls[0].c, input);
  assert.equal(calls[0].o.coarseStartup, false); assert.equal(calls[0].o.maxStartupAttempts, 1);
  assert.equal(r.converged, true); assert.equal(r.actualMach, .74);
  assert.equal(r.machContinuation.route, 'cold-target-inviscid');
  assert.equal(r.machContinuation.coldBaselineUsed, false);
  assert.ok(events.every(e => e.startupStrategy === 'target-inviscid-then-viscous' && e.actualMach === .74));
});

test('failed target attempt falls back once and reports failure without accepting provisional target flow', async () => {
  for (const throws of [false, true]) {
    const calls = [], events = [];
    const solve = await adapter((c, o) => {
      calls.push(c.mach);
      if (c.mach === .74 && throws) throw new Error('Euler state inadmissible');
      return result(c, false);
    });
    const r = solve(input, { onStage: e => events.push(e) });
    assert.deepEqual(calls, [.74, .2]); assert.equal(r.converged, false);
    assert.equal(r.actualMach, .2); assert.equal(r.targetMach, .74);
    assert.equal(r.machContinuation.targetInviscidStartup.converged, false);
    assert.equal(events[0].startupStrategy, 'low-mach-fallback');
  }
});

test('opt-out, subsonic and zero-iteration requests preserve old routing; cancellation never triggers fallback', async () => {
  const calls = [];
  const solve = await adapter((c, o) => { calls.push(c.mach); o.onStage({ stage: 'euler' }); return result(c, false); });
  solve(input, { targetInviscidStartup: false });
  solve({ ...input, mach: .2 }); solve(input, { maxIterations: 0 });
  assert.deepEqual(calls, [.2, .2, .2]); calls.length = 0;
  const cancel = new Error('cancel');
  assert.throws(() => solve(input, { onStage: () => { throw cancel; } }), e => e === cancel);
  assert.deepEqual(calls, [.74]);
  assert.throws(() => solve(input, { targetInviscidStartup: 1 }), /Invalid target inviscid/);
});

test('progress distinguishes target Euler startup, viscosity introduction and fallback', () => {
  const p = createQuadSolveProgress(input);
  assert.match(p.update({ stage: 'euler', startupStrategy: 'target-inviscid-then-viscous' }, { stageChange: true }).current,
    /target-condition inviscid startup/);
  assert.match(p.update({ stage: 'boundary-layer-initialization', startupStrategy: 'target-inviscid-then-viscous' }, { stageChange: true }).current,
    /viscous solve at requested conditions/);
  assert.match(p.update({ stage: 'coupled-cold-recovery', startupStrategy: 'low-mach-fallback' }, { stageChange: true }).current, /low-Mach fallback/);
  assert.doesNotMatch(p.update({ stage: 'euler' }, { stageChange: true }).current, /target-condition/);
});


test('public transonic route solves a nearby lower incidence before approaching the requested alpha at fixed Mach', async () => {
  const calls = [], stages = [];
  const solve = await adapter((c, o) => {
    calls.push(['seed', c.alpha, c.mach]); o.onStage({ stage: 'euler' });
    return { ...result(c, true), checkpoint: { restart: { input: { alpha: c.alpha, mach: .74, streamwiseMode: 'hybrid' }, options: {} }, continuation: {} } };
  }, (target, o) => {
    calls.push(['alpha', target, o.initialCheckpoint.restart.input.mach]);
    assert.equal(o.maxAlphaStep, .1);
    assert.equal(o.targetMach, .74);
    const checkpoint = structuredClone(o.initialCheckpoint); checkpoint.restart.input.alpha = target;
    return { checkpoint, mach: .74, alpha: target, actualAlpha: target, converged: true, alphaContinuation: { reachedTarget: true } };
  });
  const r = solve(input, { alphaContinuation: true, onStage: e => stages.push(e) });
  assert.deepEqual(calls, [['seed', input.alpha - 1, .74], ['alpha', 2.68, .74]]);
  assert.equal(r.sourceCase.alpha, 2.68); assert.equal(r.converged, true);
  assert.ok(stages.every(e => e.actualAlpha === input.alpha - 1 && e.targetAlpha === 2.68));
});

test('failed lower-incidence startup tries another lower angle and never starts final-alpha Mach-only fallback', async () => {
 const calls=[];
 const solve=await adapter((c)=>{calls.push(c.alpha);return result(c,false);});
 const r=solve(input,{alphaContinuation:true});
 assert.equal(r.converged,false);
 assert.ok(calls.includes(input.alpha-1));assert.ok(calls.includes((input.alpha-1)/2));
 assert.ok(calls.every(a=>a!==input.alpha));
 assert.match(r.reason,/Lower-incidence startup/);
 assert.equal(r.alphaContinuation.targetAlpha,input.alpha);
});
