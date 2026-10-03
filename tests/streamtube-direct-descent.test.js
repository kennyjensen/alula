// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {solveStreamtubeIses} from '../src/euler/streamtube-ises-update.js';
import {solveCoupledStreamtubeIses} from '../src/euler/streamtube-coupled-ises.js';
import {intrinsicBodyFixture} from './fixtures/intrinsic-body.js';

test('Euler Armijo checks maintained residuals and preserves checkpoint policy',()=>{
 const r=solveStreamtubeIses(intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2}),{
  stepAcceptance:'armijo',maxIterations:20,tolerance:1e-10,retainCheckpoint:true,retainBestCheckpoint:true});
 assert.equal(r.converged,true,r.reason);
 assert.ok(r.history.slice(1).every(h=>h.residualDecrease &&
  (h.residualDecrease.converged || h.residualDecrease.afterSquaredNorm<=h.residualDecrease.allowedSquaredNorm)));
 const replay=solveStreamtubeIses(undefined,{resume:r.checkpoint,stepAcceptance:'armijo',maxIterations:0,tolerance:1e-10});
 assert.equal(replay.converged,true,replay.reason);
 assert.throws(()=>solveStreamtubeIses(undefined,{resume:r.checkpoint,stepAcceptance:'admissible',maxIterations:0}),/controls/);
});

test('RAE transition event escapes a fixed-topology merit barrier without changing operating conditions',()=>{
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/direct-startup/viscous-armijo.json')).checkpoint;
 cp.continuation.stepAcceptance='event-armijo';
 const before=structuredClone(cp);
 const r=solveCoupledStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:15,tolerance:1e-10});
 assert.equal(r.converged,true,r.reason);
 assert.ok(Object.values(r.families).every(v=>v<=1e-10));
 assert.ok(r.history.some(h=>h.activeChange && h.transitionEventAcceptance));
 for(const h of r.history.slice(1)){
  if(h.transitionEventAcceptance){assert.equal(h.transitionEventAcceptance.method,'admissible-transition-event');assert.equal(h.activeChange,true);}
  else assert.ok(h.residualDecrease && (h.residualDecrease.converged || h.residualDecrease.afterSquaredNorm<=h.residualDecrease.allowedSquaredNorm));
 }
 assert.equal(r.checkpoint.restart.input.mach,.74);
 assert.equal(r.checkpoint.restart.input.alpha,2.68);
 assert.equal(r.checkpoint.restart.options.ncrit,4);
 assert.equal(r.checkpoint.restart.options.transitionMode,'automatic');
 assert.equal(r.x.length,4493);
 assert.deepEqual(cp,before);
});
