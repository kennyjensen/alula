import fs from 'node:fs';
import { createHash } from 'node:crypto';
import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { createStreamtubeBoundaryLayers } from '../src/euler/streamtube-boundary-layers.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

import { refineCoupledStreamtubeBody as refine } from '../src/euler/streamtube-coupled-refinement.js';
const baselineBytes=fs.readFileSync(new URL('../docs/rae2822/transonic-refinement-audit/finite-base-draft/production-before/streamtube-coupled-refinement.js.txt',import.meta.url));
assert.equal(createHash('sha256').update(baselineBytes).digest('hex'),'eab6845ef5339f0f2f71227915d94e6b7e4c79f229106bc5669ab638ecbaca97');
const baselineText=baselineBytes.toString().replace(/from '([^']+)'/g,(_,p)=>`from '${new URL(p,new URL('../src/euler/streamtube-coupled-refinement.js',import.meta.url))}'`);
const { refineCoupledStreamtubeBody:currentRefine }=await import('data:text/javascript;base64,'+Buffer.from(baselineText).toString('base64'));

const model='historical-common-isentrope',gamma=1.4;
const plain=x=>JSON.parse(JSON.stringify(x,(_,v)=>ArrayBuffer.isView(v)?Array.from(v):v));
const close=(a,b,tol=2e-12)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
function fixture({edgeMach=.6,elements=1,historical=true,jump=false}={}){
 const mach=historical?.5:.2;
 const input={...intrinsicBodyFixture({bodySegments:4,tubes:1,elements,mach}),normalStencil:'body-stations',stagnationMotion:'walls-only',
  streamwiseMode:historical?'hybrid':'isentropic',...(historical?{hybrid:{epsilonP:1e-5},upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}}: {})};
 const euler=createStreamtubeBodySystem(input),x=euler.initial.slice(),base=euler.evaluate(x);
 const speed=Math.sqrt((gamma-1)*edgeMach**2*euler.conditions.h0/(1+.5*(gamma-1)*edgeMach**2));
 const first=Math.max(...euler.layout.bodies.map(b=>b.trailingIndex))+1, last=euler.layout.nx-2;
 for(let i=first;i<=last;i++) for(let g=0;g<euler.layout.tubes.length;g++){
  const s=base.sections[i][g][0]; const q=jump&&i===last?speed*.8:speed;
  x[euler.layout.densityIndex(i,g,0)]=Math.log(s.rho*s.q/q);
 }
 const options={reynolds:1e6,ncrit:9,transitionMode:'automatic',tripFractions:input.bodies.map(()=>[1,1]),edgeMatching:'section-velocity',...(historical?{blThermodynamics:model}: {})};
 const bl=createStreamtubeBoundaryLayers(euler,x,{...options,...(historical?{allowSupersonicEdge:true}: {})});
 const phase=bl.surfaces.map(s=>s.ids.length-1); bl.restoreActive(phase);
 const initialBL=new Float64Array(4*bl.stations.length);
 for(const s of bl.stations){const theta=s.kind==='wake'?8e-7:4e-7,delta=theta*(2.2*(1+.113*edgeMach**2)+.29*edgeMach**2);
  initialBL.set([['transition','leading-transition','turbulent','trailing-edge','wake'].includes(s.regime)?.03:0,theta/bl.scale,delta/bl.scale,speed],4*s.id);}
 const parent=createCoupledStreamtubeBody(input,{...options,transitionState:phase,initialEuler:{x,nodes:base.nodes},initialBL});
 assert.ok(parent.admissible(parent.initial),'Manufactured parent must be admissible, not claimed converged.');
 return{input,parent,first,last};
}
function checks(f,r,streamwiseFactor=1){
 assert.deepEqual(r.system.conditions,f.parent.conditions);
 assert.deepEqual(r.system.euler.conditions,f.parent.euler.conditions);
 assert.equal(r.system.bl.transitionMode,'automatic');assert.deepEqual(r.system.bl.trips,f.parent.bl.trips);
 assert.ok(r.system.admissible(r.system.initial));assert.equal(r.diagnostics.quality.valid,true);
 const value=r.system.evaluate(r.system.initial),old=f.parent.evaluate(f.parent.initial);
 assert.ok(value.residual.every(Number.isFinite));
 const restart={input:r.input,options:{...r.options,transitionState:r.system.bl.snapshotActive()},initialEuler:r.initialEuler,initialBL:r.initialBL};
 const replay=createCoupledStreamtubeBody(restart.input,{...restart.options,initialEuler:restart.initialEuler,initialBL:restart.initialBL});
 assert.deepEqual(replay.initial,r.system.initial);assert.deepEqual(replay.evaluate(replay.initial).residual,value.residual);
 assert.deepEqual(replay.bl.snapshotActive(),r.system.bl.snapshotActive());
 for(let g=0;g<old.outer.nodes.length;g++)for(let i=0;i<old.outer.nodes[g].length;i++) for(let j=0;j<old.outer.nodes[g][i].length;j++){
  const p=old.outer.nodes[g][i][j],q=value.outer.nodes[g][streamwiseFactor*i][2*j];close(p.x,q.x);close(p.y,q.y);
 }
 return{old,value};
}
for(const edgeMach of [.6,1.2])test(`complete historical refiner retains physical density and ${edgeMach<1?'sub':'super'}sonic interior states`,()=>{
 const f=fixture({edgeMach}),before=plain({input:f.input,x:f.parent.initial}),r=refine(f.input,f.parent,{streamwiseFactor:1,normalFactor:2});
 const {old,value}=checks(f,r);assert.equal(r.options.blThermodynamics,model);assert.ok(r.diagnostics.physicalDensityTransfer);
 for(let i=0;i<f.parent.euler.layout.nx;i++)for(let g=0;g<f.parent.euler.layout.tubes.length;g++)for(let j=0;j<2;j++){
  assert.equal(r.system.initial[r.system.euler.layout.densityIndex(i,g,j)],f.parent.initial[f.parent.euler.layout.densityIndex(i,g,0)]);
  if(i>=f.first&&i<=f.last){
    close(value.outer.sections[i][g][j].q,old.outer.sections[i][g][0].q,2e-10);
    assert.equal(value.outer.sections[i][g][j].machSquared>1,edgeMach>1);
  }
 }
 assert.deepEqual(plain({input:f.input,x:f.parent.initial}),before);
});
test('two-body historical entropy variation survives full normal refinement with every phase/model flag intact',()=>{
 const f=fixture({edgeMach:1.2,elements:2,jump:true}),r=refine(f.input,f.parent,{streamwiseFactor:1,normalFactor:2}),{old,value}=checks(f,r);
 assert.equal(r.system.bl.surfaces.length,4);assert.equal(r.system.bl.wakes.length,2);assert.deepEqual(r.system.bl.snapshotActive(),f.parent.bl.snapshotActive());
 let span=0;
 for(let i=0;i<f.parent.euler.layout.nx;i++)for(let g=0;g<f.parent.euler.layout.tubes.length;g++){
  const entropy=s=>Math.log(s.enthalpy)/(gamma-1)-Math.log(s.rho),original=entropy(old.outer.sections[i][g][0]);span=Math.max(span,Math.abs(original-entropy(old.outer.sections[0][g][0])));
  if(i>=f.first&&i<=f.last)for(let j=0;j<2;j++)close(entropy(value.outer.sections[i][g][j]),original,2e-10);
 }
 assert.ok(span>.01,'This must contain a real physical entropy variation, not merely a uniform isentrope.');
});
test('omitted historical mode follows the exact pre-existing subsonic refiner path',()=>{
 const f=fixture({historical:false}),controls={streamwiseFactor:1,normalFactor:2};
 const old=currentRefine(f.input,f.parent,controls),next=refine(f.input,f.parent,controls);
 for(const key of ['input','options','initialEuler','initialBL','diagnostics'])assert.deepEqual(plain(next[key]),plain(old[key]));
 assert.deepEqual(next.system.initial,old.system.initial);assert.deepEqual(next.system.evaluate(next.system.initial).residual,old.system.evaluate(old.system.initial).residual);
});

test('simultaneous streamwise/normal historical refinement preserves geometry, BL primitives, phases and canonical child state',()=>{
 const f=fixture({edgeMach:1.2}),r=refine(f.input,f.parent,{streamwiseFactor:2,normalFactor:2}),{old,value}=checks(f,r,2);
 assert.equal(r.system.euler.layout.nx,2*f.parent.euler.layout.nx);
 assert.deepEqual(r.options.transitionState,r.system.bl.snapshotActive());
 assert.ok(r.diagnostics.physicalDensityTransfer.minimumStaticEnthalpy>0);
 for(const s of f.parent.bl.stations){
  const child=r.system.bl.stations.find(t=>t.kind===s.kind&&t.body===s.body&&t.side===s.side&&t.i===2*s.i);
  assert.ok(child);
  for(const key of ['theta','deltaStar','ue'])assert.equal(value.layers.states[child.id][key],old.layers.states[s.id][key]);
 }
 assert.equal(r.system.conditions.blThermodynamics,model);
});

for(const mode of ['hybrid','momentum'])test(`${mode} upwind Euler density transfer is independent of the subsonic BL thermodynamic adapter`,()=>{
 const f=fixture({historical:false}),old=f.parent.evaluate(f.parent.initial),input={...f.input,streamwiseMode:'hybrid',hybrid:{epsilonP:1e-5},upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}};
 if(mode==='momentum'){input.streamwiseMode='momentum';delete input.hybrid;}
 const parent=createCoupledStreamtubeBody(input,{reynolds:1e6,ncrit:9,edgeMatching:'section-velocity',transitionMode:'automatic',tripFractions:f.parent.bl.trips,
  transitionState:f.parent.bl.snapshotActive(),initialEuler:{x:f.parent.initial.slice(0,f.parent.ne),nodes:old.outer.nodes,undisplacedNodes:old.outer.undisplacedNodes},initialBL:f.parent.initial.slice(f.parent.ne)});
 const r=refine(input,parent,{streamwiseFactor:1,normalFactor:2});checks({input,parent},r);
 assert.equal(r.options.blThermodynamics,undefined);assert.ok(r.diagnostics.physicalDensityTransfer);
 for(let i=0;i<parent.euler.layout.nx;i++)for(let g=0;g<parent.euler.layout.tubes.length;g++)for(let j=0;j<2;j++)
  assert.equal(r.system.initial[r.system.euler.layout.densityIndex(i,g,j)],parent.initial[parent.euler.layout.densityIndex(i,g,0)]);
});
