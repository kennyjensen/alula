// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { coupledStreamtubeTripEvents } from '../src/euler/streamtube-coupled.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { selectSurfaceTransition, prepareSurfaceTransition } from '../src/viscous/transition-selection.js';

function synthetic({ mode='automatic',kind='natural',from=1,to=4,fraction=.8 }={}) {
 const target={body:0,side:'upper',from,to,fraction,kind}, x=Float64Array.from([.4,...Array.from({length:6},()=>[0,1,2,1]).flat()]);
 let calls=0;
 const expected={changed:true,changes:[target],sentinel:'unchanged transfer result'};
 const system={ne:1,bl:{transitionMode:mode,kernel:{parameters:{ncrit:9}},surfaces:[{ids:[0,1,2,3,4,5],transition:from}],
  activeTargets:()=>[target],snapshotActive:()=>[from],restoreActive:()=>{},updateActive:()=>{calls++;return expected;}}};
 return {system,x,target,expected,calls:()=>calls};
}

test('XFOIL automatic natural selection can cross several intervals in both directions and a full new interval',()=>{
 for(const [from,to,fraction] of [[1,4,.8],[4,1,.1],[1,2,.9],[2,1,.05]]){
  const f=synthetic({from,to,fraction}),current=f.x.slice();
  assert.throws(()=>coupledStreamtubeTripEvents(f.system).prepare(f.x,current),/crosses too much/);
  assert.equal(f.calls(),0);
  assert.equal(coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,current),f.expected);
  assert.equal(f.calls(),1);assert.deepEqual(f.x,current);
 }
});

test('fixed trips, interior forced automatic targets and legacy drivers retain the original migration gate',()=>{
 for(const [mode,kind,blUpdate] of [['fixed-trip','natural','xfoil'],['fixed-trip',undefined,'giles'],
  ['automatic','forced','xfoil'],['automatic','trailing-edge','giles'],['automatic','natural','giles']]){
  const f=synthetic({mode,kind}),before=f.x.slice();
  assert.throws(()=>coupledStreamtubeTripEvents(f.system,{blUpdate}).prepare(f.x,before),error=>{
   assert.match(error.message,/crosses too much/);assert.equal(error.code,'COUPLED_TRANSITION_MIGRATION_LIMIT');
   assert.deepEqual(error.diagnostics.targets,[f.target]);return true;
  });
  assert.equal(f.calls(),0);assert.deepEqual(f.x,before);
 }
 const f=synthetic({mode:'fixed-trip',kind:'forced',from:1,to:2,fraction:.2});
 assert.equal(coupledStreamtubeTripEvents(f.system).prepare(f.x,f.x.slice()),f.expected);
 assert.throws(()=>coupledStreamtubeTripEvents(f.system,{blUpdate:'unknown'}),/Unknown BL update policy/);
});

test('automatic XFOIL terminal reselection crosses the full surface with its existing global speed bound',()=>{
 for(const from of [1,4]){
  const f=synthetic({kind:'trailing-edge',from,to:5,fraction:1}),current=f.x.slice();f.x[4]+=.375;
  assert.equal(coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,current),f.expected);
  assert.equal(f.calls(),1);
 }
 for(const [slot,value,message] of [[2,.49,/thickness/],[3,.99,/thickness/],[4,1.376,/edge velocity/]]){
  const f=synthetic({kind:'trailing-edge',from:1,to:5,fraction:1}),current=f.x.slice();f.x[slot]=value;
  assert.throws(()=>coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,current),message);assert.equal(f.calls(),0);
 }
});

test('natural index freedom retains the half-thickness and global XFOIL edge-speed bounds',()=>{
 for(const [slot,value,message] of [[2,.49,/thickness/],[3,.99,/thickness/],[4,1.376,/edge velocity/]]){
  const f=synthetic(),current=f.x.slice();f.x[slot]=value;const before=f.x.slice();
  assert.throws(()=>coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,current),message);
  assert.equal(f.calls(),0);assert.deepEqual(f.x,before);
 }
});

test('unchanged automatic phase is an exact no-op under both policies',()=>{
 for(const blUpdate of ['giles','xfoil']){
  const f=synthetic({from:2,to:2}),before=f.x.slice();
  assert.deepEqual(coupledStreamtubeTripEvents(f.system,{blUpdate}).prepare(f.x,before),{changed:false});
  assert.equal(f.calls(),0);assert.deepEqual(f.x,before);
 }
});

const native=JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/automatic-transition.json',import.meta.url)));
const profile=native.profiles.find(p=>p.parameters.ncrit===9&&!p.tripS);
// A local profile adapter isolates event selection and real native-derived
// auxiliary conversion, with no Euler system, Jacobian or Newton solve.
function physical(previousIndex,{ncrit=profile.parameters.ncrit,tripS}={}){
 const controls=tripS===undefined?{}:{tripS};
 const kernel=createIntegralKernel({...profile.parameters,ncrit,exactJacobian:true}),ids=profile.states.map((_,i)=>i);
 const surface={body:0,side:'upper',ids,transition:previousIndex};
 const x=Float64Array.from([.4,...profile.states.flatMap((s,i)=>[i<previousIndex?0:.03,s.theta,s.deltaStar,s.ue])]);
 const decode=packed=>profile.states.map((s,j)=>({...s,aux:packed[4*j],theta:packed[4*j+1],deltaStar:packed[4*j+2],ue:packed[4*j+3]}));
 const target=packed=>{
  const t=selectSurfaceTransition(kernel,decode(packed),controls),down=profile.states[t.index].s,up=profile.states[t.index-1].s;
  return {body:0,side:'upper',from:surface.transition,to:t.index,kind:t.kind==='forced'&&tripS===profile.states.at(-1).s?'trailing-edge':t.kind,fraction:(t.s-up)/(down-up)};
 };
 const bl={transitionMode:'automatic',kernel,surfaces:[surface],snapshotActive:()=>[surface.transition],restoreActive:a=>{surface.transition=a[0];},
  activeTargets:(_,packed)=>[target(packed)],updateActive:packed=>{
   const t=target(packed),p=prepareSurfaceTransition(kernel,decode(packed),{previousIndex:surface.transition,...controls});
   p.auxiliary.forEach((value,j)=>{packed[4*j]=value;});surface.transition=p.index;
   return {changed:p.changed,changes:[{...t,converted:p.converted}]};
  }};
 return {system:{ne:1,bl},x,decode};
}

test('real native-derived natural selection converts every crossed N/shear slot without changing physical primitives',()=>{
 for(const previous of [1,8]){
  const f=physical(previous),before=f.x.slice(),packed=f.x.subarray(1),expected=prepareSurfaceTransition(f.system.bl.kernel,f.decode(packed),{previousIndex:previous});
  assert.equal(expected.index,profile.expected.index);assert.equal(expected.kind,'natural');
  assert.ok(Math.abs(previous-expected.index)>1);
  const result=coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,before);
  assert.equal(result.changed,true);assert.equal(result.changes[0].converted.length,Math.abs(previous-expected.index)+Number(previous<expected.index));
  assert.deepEqual(f.system.bl.snapshotActive(),[expected.index]);
  for(let j=0;j<profile.states.length;j++){
   assert.equal(packed[4*j],expected.auxiliary[j]);
   for(let k=1;k<4;k++)assert.equal(packed[4*j+k],before[1+4*j+k]);
   if(j<expected.index)assert.ok(packed[4*j]>=0&&packed[4*j]<9);else assert.ok(packed[4*j]>0);
  }
  assert.equal(f.x[0],before[0]);
  const stable=f.x.slice();assert.deepEqual(coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,stable),{changed:false});
  assert.deepEqual(f.x,stable);
 }
});

test('real automatic terminal selection converts the turbulent prefix back to N without changing physical primitives',()=>{
 const tripS=profile.states.at(-1).s,previous=1,f=physical(previous,{ncrit:100,tripS}),before=f.x.slice(),packed=f.x.subarray(1);
 const expected=prepareSurfaceTransition(f.system.bl.kernel,f.decode(packed),{previousIndex:previous,tripS});
 assert.equal(expected.kind,'forced');assert.equal(expected.index,profile.states.length-1);assert.ok(expected.index-previous>1);
 const target=f.system.bl.activeTargets(null,packed)[0];assert.equal(target.kind,'trailing-edge');assert.equal(target.fraction,1);
 const result=coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,before);
 assert.equal(result.changed,true);assert.equal(result.changes[0].converted.length,expected.index-previous+1);assert.deepEqual(f.system.bl.snapshotActive(),[expected.index]);
 assert.equal(result.changes[0].converted.at(-1).terminalShearInitialization,true);
 for(let j=0;j<profile.states.length;j++){
  assert.equal(packed[4*j],expected.auxiliary[j]);for(let k=1;k<4;k++)assert.equal(packed[4*j+k],before[1+4*j+k]);
  if(j<expected.index)assert.ok(packed[4*j]>=0&&packed[4*j]<100);else assert.ok(packed[4*j]>0);
 }
 assert.equal(f.x[0],before[0]);const stable=f.x.slice();assert.deepEqual(coupledStreamtubeTripEvents(f.system,{blUpdate:'xfoil'}).prepare(f.x,stable),{changed:false});assert.deepEqual(f.x,stable);
});
