import test from 'node:test';
import assert from 'node:assert/strict';
import {createStreamtubeBodySystem} from '../src/euler/streamtube-body.js';
import {intrinsicBodyFixture} from './fixtures/intrinsic-body.js';
import {sampleStreamtubePanelVelocity,initializeStreamtubePanelVelocity} from '../src/euler/streamtube-panel-topology.js';

test('panel transfer preserves geometry, captured mass, globals and total enthalpy while reproducing section speeds',()=>{
 const system=createStreamtubeBodySystem({...intrinsicBodyFixture({elements:1,bodySegments:4,tubes:2,alpha:0}),mach:.2});
 const x=system.initial.slice(),before=system.decode(x),sample={method:'uncorrected-panel-velocity',speeds:Array(system.layout.densityCount).fill(.9),maximumDirectionMismatch:0};
 const r=initializeStreamtubePanelVelocity(system,x,sample),after=system.decode(r.initial);
 assert.deepEqual(system.initial,x);assert.deepEqual(before.nodes,after.nodes);assert.deepEqual(before.allocation,after.allocation);
 assert.deepEqual(r.initial.slice(system.layout.densityCount),x.slice(system.layout.densityCount));
 for(const row of r.flow.sections)for(const group of row)for(const s of group){
  assert.ok(Math.abs(s.q-.9)<1e-14);
  assert.ok(Math.abs(s.enthalpy+.5*s.q*s.q-system.conditions.h0)<1e-13);
  assert.ok(s.rho>0&&s.p>0);
 }
 assert.equal(r.diagnostics.pgCorrection,false);assert.ok(r.diagnostics.maximumSpeedError<1e-14);
});

test('panel transfer rejects reversed/nonfinite velocities and impossible gas states without clipping',()=>{
 const s=createStreamtubeBodySystem({...intrinsicBodyFixture({elements:1,bodySegments:4,tubes:2,alpha:0}),mach:.2}),x=s.initial.slice();
 assert.throws(()=>sampleStreamtubePanelVelocity(s,x,()=>({u:NaN,v:0})),/forward and finite/);
 const sample={method:'uncorrected-panel-velocity',speeds:Array(s.layout.densityCount).fill(100)};
 assert.throws(()=>initializeStreamtubePanelVelocity(s,x,sample),/enthalpy/);
 sample.speeds[0]=-1;assert.throws(()=>initializeStreamtubePanelVelocity(s,x,sample),/Invalid panel velocity/);
 assert.deepEqual(s.initial,x);
});

test('panel sampler projects velocity onto the final section direction without a PG speed multiplier',()=>{
 const nodes=[Array.from({length:3},(_,i)=>[{x:i,y:0},{x:i,y:1}])];
 const system={layout:{nx:2,tubes:[1],densityCount:2,densityIndex:i=>i},
  decode:()=>({nodes,allocation:{groups:[[{massFlow:1}]]}})};
 const sample=sampleStreamtubePanelVelocity(system,[],()=>({u:2,v:1}));
 assert.deepEqual(sample.speeds,[2,2]);
 assert.ok(Math.abs(sample.maximumDirectionMismatch-1/Math.sqrt(5))<1e-15);
 assert.throws(()=>sampleStreamtubePanelVelocity(system,[],()=>({u:-1,v:0})),/forward and finite/);
});
