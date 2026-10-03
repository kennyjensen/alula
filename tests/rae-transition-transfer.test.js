// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {initializeCoupledStreamtubeFromFlow} from '../src/euler/tests/streamtube-coupled-flow-restart.js';
import {solveCoupledStreamtubeIses} from '../src/euler/streamtube-coupled-ises.js';
import {createCoupledStreamtubeBody} from '../src/euler/streamtube-coupled.js';
const read=()=>JSON.parse(fs.readFileSync('docs/solver-reliability/rae-small-alpha-audit/walk-root.json')).checkpoint;

test('RAE natural-transition crossing at 1.75 to 1.85 converges without replaying downstream shear',()=>{
 const cp=read(),before=structuredClone(cp),seed=initializeCoupledStreamtubeFromFlow(.74,cp,{targetAlpha:1.85});
 const r=solveCoupledStreamtubeIses(undefined,{resume:seed.checkpoint,...seed.checkpoint.continuation,
  maxIterations:15,tolerance:1e-10,iterationRecovery:false});
 assert.equal(r.converged,true,r.reason);
 assert.ok(Object.values(r.families).every(v=>v<=1e-10));
 assert.ok(r.history.some(h=>h.changes?.some(c=>c.downstreamShearPreserved)));
 assert.equal(r.checkpoint.restart.input.mach,.74);
 assert.equal(r.checkpoint.restart.input.alpha,1.85);
 assert.equal(r.checkpoint.restart.options.ncrit,4);
 assert.equal(r.checkpoint.restart.options.transitionMode,'automatic');
 assert.deepEqual(cp,before);
});

test('an infinitesimal natural-transition boundary crossing leaves existing downstream shear continuous',()=>{
 const cp=read(),f=cp.restart;
 const evidence=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-transition-event/continuity.json'));
 for(const epsilon of [1e-4,1e-6,1e-8]){
  const sides=[];
  for(const sign of [-1,1]){
   const s=createCoupledStreamtubeBody(f.input,{...f.options,ncrit:evidence.criterionAtBoundary+sign*epsilon,initialEuler:f.initialEuler,initialBL:f.initialBL});
   const x=s.initial.subarray(s.ne),before=x.slice();
   s.bl.updateActive(x,s.initial.subarray(0,s.ne));
   assert.ok(x.every((v,i)=>i%4===0||v===before[i]));
   sides.push({interval:s.bl.snapshotActive()[0],shear:s.bl.surfaces[0].ids.map(id=>x[4*id])});
  }
  assert.deepEqual(sides.map(s=>s.interval),[37,38]);
  assert.deepEqual(sides[0].shear.slice(39),sides[1].shear.slice(39));
 }
});
