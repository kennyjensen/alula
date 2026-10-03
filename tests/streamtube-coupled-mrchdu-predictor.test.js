import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { prepareCoupledMrchduProfiles } from '../src/euler/streamtube-coupled-mrchdu-predictor.js';
import { initializeCoupledStreamtubeFromFlow } from '../src/euler/tests/streamtube-coupled-flow-restart.js';
const serial=v=>JSON.stringify(v,(_,x)=>ArrayBuffer.isView(x)?Array.from(x):x);
const fixture=JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-body.json',import.meta.url)));
const cases=fixture.cases.filter(c=>c.name.startsWith('default-two-element-root'));
function bodyProfiles(){
  const scale=1/Math.sqrt(cases[0].parameters.reynolds),states=[],initialBL=[],surfaces=[],wakes=[],expected=[];
  const put=p=>p.map(a=>{const id=states.length;initialBL.push(a.aux,a.theta/scale,a.deltaStar/scale,a.ue);
    states.push({...a,theta:scale*initialBL[4*id+1],deltaStar:scale*initialBL[4*id+2]});return id;});
  cases.forEach((c,body)=>{
    const wake={body,ids:put(c.input.wake)};wakes.unshift(wake);
    expected.push({ids:wake.ids,states:c.native.wake});
    for(const side of [1,0]){
      const surface={body,side:['upper','lower'][side],ids:put(c.input.surfaces[side]),transition:c.input.phases[side]};
      surfaces.unshift(surface);expected.push({ids:surface.ids,states:c.native.surfaces[side].states});
    }
  });
  return {bl:{transitionMode:'automatic',surfaces,wakes,scale,trips:[[1,1],[1,1]],kernel:{parameters:cases[0].parameters}},
    states,initialBL,targetMach:cases[0].parameters.mach,expected};
}
test('mixed predictor preserves native physical profiles through reordered body/station metadata and packed scaling',()=>{
  const input=bodyProfiles(),before=serial(input),r=prepareCoupledMrchduProfiles(input);
  assert.equal(serial(input),before);let maximum=0;
  for(const p of input.expected)p.states.forEach((a,k)=>{
    const id=p.ids[k],actual={aux:r.initialBL[4*id],theta:input.bl.scale*r.initialBL[4*id+1],deltaStar:input.bl.scale*r.initialBL[4*id+2],ue:r.initialBL[4*id+3]};
    for(const name of ['aux','theta','deltaStar','ue']){
      const error=Math.abs(actual[name]-a[name])/(name==='aux'?Math.max(.01,Math.abs(a[name])):Math.abs(a[name]));
      maximum=Math.max(maximum,error);assert.ok(error<1e-10,`Native ${name} discrepancy at station ${id}: ${error}`);
    }
  });
  r.transitionState.forEach((phase,k)=>{const s=input.bl.surfaces[k];assert.equal(phase,cases[s.body].native.surfaces[s.side==='upper'?0:1].transition);});
  assert.ok(r.diagnostics.bodies.every(b=>b.localConvergenceWarnings.length===0));
  assert.ok(r.diagnostics.bodies.every(b=>b.inputShapeRecovery===undefined));assert.equal(r.diagnostics.operations.translatedMRCHDUCalls,2);
  assert.equal(r.diagnostics.operations.globalLinearSolves,0);assert.ok(maximum<1e-10);
});
test('mixed predictor rejects unsupported trips, mismatched profiles and thermal failure without source mutation',()=>{
  for(const [change,pattern] of [
    [x=>{x.bl.trips[0][0]=.5;},/terminal material trips/],
    [x=>{x.bl.transitionMode='fixed-trip';},/automatic/],
    [x=>{x.states[0].theta*=2;},/packed physical BL/],
    [x=>{x.targetMach=.99;x.states.forEach((p,i)=>{p.ue=10;x.initialBL[4*i+3]=10;});},/thermal|BL edge|enthalpy/],
  ]){const x=bodyProfiles();change(x);const before=serial(x);assert.throws(()=>prepareCoupledMrchduProfiles(x),pattern);assert.equal(serial(x),before);}
});
const cp=JSON.parse(fs.readFileSync('docs/rae2822/current-rae2822-automatic-slor-mach074/checkpoints/000-accepted.json')).checkpoint;
const originalSource=fs.readFileSync('docs/rae2822/target-mrchdu-integration/before/streamtube-coupled-flow-restart.js','utf8');
const old=await import('data:text/javascript;base64,'+Buffer.from(originalSource.replace(/from '([^']+)'/g,(_,p)=>`from '${pathToFileURL(path.resolve('src/euler',p)).href}'`)).toString('base64'));
const extract=p=>({initial:p.initial,residual:p.value.residual,checkpoint:p.checkpoint,diagnostics:p.diagnostics,nodes:p.value.outer.nodes});
test('same-Mach opt-in and omitted predictor retain exact archived replay and source state',()=>{
  const source=structuredClone(cp),before=serial(source),baseline=old.initializeCoupledStreamtubeFromFlow(.2,source);
  const current=initializeCoupledStreamtubeFromFlow(.2,source),selected=initializeCoupledStreamtubeFromFlow(.2,source,{blPredictor:'xfoil-mrchdu'});
  const preserve=initializeCoupledStreamtubeFromFlow(.2,source,{blPredictor:'preserve'});
  assert.ok(serial(extract(current))===serial(extract(baseline)),'Omitted same-Mach result changed.');
  assert.ok(serial(extract(selected))===serial(extract(baseline)),'Opt-in same-Mach result changed.');assert.equal(serial(source),before);
  assert.ok(serial(extract(preserve))===serial(extract(baseline)),'Explicit preserve same-Mach result changed.');
});
test('omitted changed-Mach preparation is exactly archived and invalid predictor controls fail explicitly',()=>{
  const source=structuredClone(cp),before=serial(source);
  const current=initializeCoupledStreamtubeFromFlow(.29,source),baseline=old.initializeCoupledStreamtubeFromFlow(.29,source);
  const preserve=initializeCoupledStreamtubeFromFlow(.29,source,{blPredictor:'preserve'});
  assert.ok(serial(extract(current))===serial(extract(baseline)),'Omitted changed-Mach result changed.');assert.equal(serial(source),before);
  assert.ok(serial(extract(preserve))===serial(extract(baseline)),'Explicit preserve changed-Mach result changed.');
  assert.throws(()=>initializeCoupledStreamtubeFromFlow(.21,source,{blPredictor:'unknown'}),/Unknown coupled Mach BL predictor/);
  const trip=structuredClone(source);trip.restart.options.tripFractions[0][0]=.5;const tripBefore=serial(trip);
  assert.throws(()=>initializeCoupledStreamtubeFromFlow(.21,trip,{blPredictor:'xfoil-mrchdu'}),/terminal material trips/);assert.equal(serial(trip),tripBefore);
  // A deliberately incomplete fixed-trip checkpoint must reach the shared
  // completeness guard, not an automatic-only predictor gate.
  const fixed={version:1,restart:{options:{transitionMode:'fixed-trip',tripFractions:[[.1,.1]]}},continuation:{},families:{}};
  assert.throws(()=>initializeCoupledStreamtubeFromFlow(.21,fixed,{blPredictor:'preserve'}),/complete finite Euler\/BL/);
});
