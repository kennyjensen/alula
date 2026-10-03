import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';
import{createCoupledStreamtubeBody,coupledStreamtubeTripEvents}from'../src/euler/streamtube-coupled.js';
const base=new URL('../docs/nlr-finite-base/quad-coupled-control/',import.meta.url),read=p=>JSON.parse(fs.readFileSync(new URL(p,base)));
const accepted=read('checkpoints/iteration-008.json').checkpoint,failed=read('rejected-candidate-qualified-state.json'),reference=read('bl-reconciliation-local.json');
const serial=v=>JSON.stringify(v,(_,x)=>ArrayBuffer.isView(x)?Array.from(x):x);
const make=r=>createCoupledStreamtubeBody(r.input,{...r.options,initialEuler:{...r.initialEuler,x:Float64Array.from(r.initialEuler.x)},initialBL:Float64Array.from(r.initialBL)});

test('valid frozen simultaneous state preserves packed N, complete residual and transition exactly',()=>{
 const system=make(accepted.restart),x=system.initial.slice(),before=x.slice(),phase=system.bl.snapshotActive();
 const value=system.evaluate(x),event=coupledStreamtubeTripEvents(system).prepare(x,before);
 assert.deepEqual(event,{changed:false});assert.ok(x.every((v,k)=>v===before[k]));assert.deepEqual(system.bl.snapshotActive(),phase);
 assert.equal(serial(value.residual),serial(read('result.json').residual));assert.equal(serial(system.residual(x)),serial(value.residual));
 assert.equal(x[system.ne+4*122],8.99999403709203,'valid Newton N must not be remarched');
});

test('actual rejected NLR trial reconciles only affected surface amplification and passes unchanged equations',t=>{
 const system=make(failed.restart),x=system.initial.slice(),before=x.slice(),phase=system.bl.snapshotActive();
 assert.equal(x[system.ne+4*122],9.000014120161401);
 assert.throws(()=>system.evaluate(x),/Amplification left the active laminar interval/);
 const targets=system.bl.activeTargets(x.subarray(0,system.ne),x.subarray(system.ne));assert.ok(targets.every(s=>s.from===s.to));
 const event=coupledStreamtubeTripEvents(system).prepare(x,before);
 assert.equal(event.changed,true);assert.equal(event.changes.length,1);
 assert.deepEqual([event.changes[0].body,event.changes[0].side,event.changes[0].kind,event.changes[0].auxiliaryOnly],[1,'upper','amplification-reconciliation',true]);
 assert.equal(event.changes[0].from,event.changes[0].to);assert.deepEqual(system.bl.snapshotActive(),phase);
 assert.equal(serial(x.subarray(system.ne)),serial(reference.correctedRestart.initialBL));
 const affected=new Set(system.bl.surfaces[2].ids.slice(0,phase[2]).map(id=>system.ne+4*id));
 assert.ok(x.every((v,k)=>affected.has(k)||v===before[k]),'geometry, theta, delta-star, Ue, other surfaces, wakes and turbulent shear are exact');
 const value=system.evaluate(x);assert.ok(value.residual.every(Number.isFinite));assert.equal(system.admissible(x),true);
 const maximum=Math.max(...system.bl.surfaces[2].ids.slice(0,phase[2]).map(id=>Math.abs(value.residual[system.ne+4*id])));
 assert.ok(maximum<1e-10);assert.ok(x[system.ne+4*122]<system.bl.kernel.parameters.ncrit);
 t.diagnostic(JSON.stringify({beforeN:before[system.ne+4*122],afterN:x[system.ne+4*122],maxNativeAmplificationResidual:maximum,families:value.families,changedAuxiliarySlots:event.changes[0].stations.length}));
});

test('failed selective native preparation is atomic and invalid surface controls remain rejected',()=>{
 const system=make(failed.restart),x=system.initial.slice(),before=x.slice(),phase=system.bl.snapshotActive();
 assert.throws(()=>system.bl.updateActive(x.subarray(system.ne),x.subarray(0,system.ne),{reconcileAmplificationSurfaces:[4]}),/Invalid amplification/);
 const original=system.bl.kernel.station;system.bl.kernel.station=()=>{throw new Error('Manufactured native domain failure');};
 assert.throws(()=>coupledStreamtubeTripEvents(system).prepare(x,before),/Manufactured native domain failure/);
 system.bl.kernel.station=original;assert.ok(x.every((v,k)=>v===before[k]));assert.deepEqual(system.bl.snapshotActive(),phase);
});
