import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { initializeCoupledStreamtubeFromFlow } from '../src/euler/tests/streamtube-coupled-flow-restart.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
const serial=v=>JSON.parse(JSON.stringify(v,(_,x)=>ArrayBuffer.isView(x)?Array.from(x):x));
const savedPath='docs/rae2822/current-rae2822-automatic-slor-mach074/checkpoints/000-accepted.json';
const source=JSON.parse(fs.readFileSync(savedPath)).checkpoint;
const originalSource=fs.readFileSync('docs/rae2822/mach-target-phase/streamtube-coupled-flow-restart.before.js.txt','utf8');
const old=await import('data:text/javascript;base64,'+Buffer.from(originalSource.replace(/from '([^']+)'/g,(_,p)=>`from '${pathToFileURL(path.resolve('src/euler',p)).href}'`)).toString('base64'));
const extract=p=>serial({initial:p.initial,residual:p.value.residual,checkpoint:p.checkpoint,diagnostics:p.diagnostics,nodes:p.value.outer.nodes,undisplacedNodes:p.value.outer.undisplacedNodes});
const difference=(a,b,p='')=>{if(a===b)return null;if(!a||!b||typeof a!=='object'||typeof b!=='object')return {path:p,a,b};for(const k of new Set([...Object.keys(a),...Object.keys(b)])){const d=difference(a[k],b[k],p+'/'+k);if(d)return d;}return null;};
const equal=(a,b,label)=>assert.equal(difference(serial(a),serial(b)),null,label);

test('same-Mach automatic replay keeps the exact archived result and source state',()=>{
 const cp=structuredClone(source),before=serial(cp);
 const p=initializeCoupledStreamtubeFromFlow(.2,cp),original=old.initializeCoupledStreamtubeFromFlow(.2,cp);
 equal(extract(p),extract(original),'same-Mach old/new');
 equal(p.checkpoint,cp,'complete checkpoint');equal(cp,before,'source isolated');
 assert.equal(p.diagnostics.targetPhaseInitialization,undefined);
});

test('a changed-Mach automatic phase uses only native auxiliary transfer and produces a replayable complete checkpoint',t=>{
 const cp=structuredClone(source),before=serial(cp),p=initializeCoupledStreamtubeFromFlow(.29,cp);
 const phase=p.diagnostics.targetPhaseInitialization;
 assert.ok(phase,'This frozen profile must change its active interval at the selected target.');
 assert.notDeepEqual(phase.after,phase.before);
 assert.ok(phase.auxiliaryChanges.length>0);assert.equal(p.diagnostics.packedBLPreserved,false);
 assert.equal(p.diagnostics.transitionMapPreserved,false);
 const r=p.checkpoint.restart;
 for(let i=0;i<r.initialBL.length;i++)if(i%4!==0)assert.equal(r.initialBL[i],cp.restart.initialBL[i]);
 equal(r.initialEuler,cp.restart.initialEuler,'all saved Euler/geometry arrays');
 equal(r.options.tripFractions,cp.restart.options.tripFractions,'material trip limits');
 equal(p.checkpoint.continuation,cp.continuation,'maintenance and update policy');
 equal(p.system.initial,p.initial,'prepared system uses prepared auxiliaries');
 assert.ok(p.value.residual.every(Number.isFinite));
 const {iterationGeometry,stepAcceptance,stagnationLimiter}=p.checkpoint.continuation;
 const replay=solveCoupledStreamtubeIses(undefined,{resume:p.checkpoint,iterationGeometry,stepAcceptance,stagnationLimiter,maxIterations:0,tolerance:1e-10});
 equal(replay.residual,p.value.residual,'complete target residual replay');
 equal(replay.checkpoint,p.checkpoint,'complete target checkpoint replay');
 assert.equal(replay.linearDiagnostics.solves,0);assert.equal(replay.initialRedistribution.resumed,true);
 equal(cp,before,'accepted source unchanged');
 t.diagnostic(JSON.stringify({sourceMach:.2,targetMach:.29,phases:{before:phase.before,after:phase.after,auxiliaryChanges:phase.auxiliaryChanges.length},maximumResidual:p.diagnostics.targetResidual,zeroReplayLinearSolves:0}));
});

test('old fixed-trip Mach transfer arithmetic remains exactly unchanged',t=>{
 const input={...intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2,mach:.2}),streamwiseMode:'hybrid',hybrid:{epsilonP:1e-5},upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}};
 const root=solveCoupledStreamtubeIses(input,{edgeMatching:'section-velocity',blThermodynamics:'historical-common-isentrope',transitionMode:'fixed-trip',maxIterations:12,tolerance:1e-10,stepAcceptance:'admissible'});
 assert.equal(root.converged,true,root.reason);assert.ok(root.x.length<500);
 const cp=serial(root.checkpoint),before=serial(cp),p=initializeCoupledStreamtubeFromFlow(.21,cp),original=old.initializeCoupledStreamtubeFromFlow(.21,cp);
 equal(extract(p),extract(original),'fixed-trip old/new');equal(cp,before,'fixed-trip source isolated');
 assert.equal(p.diagnostics.targetPhaseInitialization,undefined);
 t.diagnostic(JSON.stringify({unknowns:root.x.length,sourceUpdates:root.history.length-1}));
});

test('phase preparation does not bypass the actual frozen raw-Hk failure or invalid source geometry/thermal states',()=>{
 const evidence=JSON.parse(fs.readFileSync('docs/rae2822/mach-target-phase/frozen-target.json'));
 assert.equal(evidence.passed,true);assert.equal(evidence.targetMach,.335);assert.equal(evidence.admissible,false);
 assert.equal(evidence.failure.code,'BL_EDGE_STATE_DOMAIN');assert.equal(evidence.failure.diagnostics.condition,'raw-hk');
 assert.ok(evidence.failure.diagnostics.rawHk<1);
 assert.deepEqual(evidence.failure.targetPhaseInitialization.before,[17,20]);assert.deepEqual(evidence.failure.targetPhaseInitialization.after,[40,21]);
 const folded=structuredClone(source);folded.restart.initialEuler.nodes[0][2][1].x=folded.restart.initialEuler.nodes[0][4][1].x;
 const foldedBefore=serial(folded);assert.throws(()=>initializeCoupledStreamtubeFromFlow(.29,folded),/valid convex physical grid/);equal(folded,foldedBefore,'folded source unchanged');
 const thermal=structuredClone(source);thermal.restart.initialBL[3]=100;const thermalBefore=serial(thermal);
 assert.throws(()=>initializeCoupledStreamtubeFromFlow(.29,thermal),/thermal|BL|coupled/i);equal(thermal,thermalBefore,'thermal source unchanged');
});
