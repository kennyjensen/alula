import test from 'node:test';
import assert from 'node:assert/strict';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody as refine } from '../src/euler/streamtube-coupled-refinement.js';
const close=(a,b,t=2e-12,label='')=>assert.ok(Math.abs(a-b)<=t,`${label}: ${a} != ${b}`);
const plain=v=>JSON.parse(JSON.stringify(v,(_,v)=>ArrayBuffer.isView(v)?[...v]:v));
function fixture({automatic=false}={}){
 const f=twoActiveFiniteBaseWakes(),old=f.system.evaluate(f.x);
 const input={...f.input,streamwiseMode:'hybrid',hybrid:{epsilonP:1e-5},upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}};
 const initialBL=f.x.slice(f.system.ne);
 if(automatic)for(const s of f.system.bl.stations)initialBL[4*s.id]=s.kind==='surface'&&s.i<input.bodies[s.body].trailingIndex?0:.03;
 const options={reynolds:1e6,ncrit:9,edgeMatching:'section-velocity',blThermodynamics:'historical-common-isentrope',
  ...(automatic?{transitionMode:'automatic',tripFractions:[[1,1],[1,1]],transitionState:f.system.bl.surfaces.map(s=>s.ids.length-1)}:{})};
 const system=createCoupledStreamtubeBody(input,{...options,initialEuler:{x:f.x.slice(0,f.system.ne),nodes:old.outer.nodes,undisplacedNodes:old.outer.undisplacedNodes},initialBL});
 assert.ok(system.admissible(system.initial),'Finite two-tail parent is an admissible manufactured off-root state.');
 return{input,system};
}
function independentGap(base,arc,L){
 if(arc===0)return base.width/L;
 const a=base.upperDerivative,b=base.lowerDerivative,c=(a.x*b.y-a.y*b.x)/(Math.hypot(a.x,a.y)*Math.hypot(b.x,b.y));
 const k=Math.max(-1.2,Math.min(1.2,c/Math.sqrt(1-c*c))),g=base.width,ell=2.5*g,t=arc/ell;
 return arc>=ell?0:(g+k*arc-(3*g+2*k*ell)*t*t+(2*g+k*ell)*t*t*t)/L;
}
function check(f,r,factor){
 const old=f.system.evaluate(f.system.initial),v=r.system.evaluate(r.system.initial),L=r.system.euler.conditions.lengthScale;
 assert.deepEqual(r.system.conditions,f.system.conditions);assert.deepEqual(r.system.euler.conditions,f.system.euler.conditions);
 assert.ok(r.system.admissible(r.system.initial));assert.equal(r.diagnostics.quality.valid,true);
 assert.equal(r.system.bl.surfaces.length,4);assert.equal(r.system.bl.wakes.length,2);
 assert.deepEqual(r.input.bodies.map(b=>b.points),f.input.bodies.map(b=>b.points));
 let maxNaiveError=0,maximumRetainedTotalDeltaChange=0,active=[];
 for(const wake of r.system.bl.wakes){
  const original=f.system.bl.wakes.find(w=>w.body===wake.body),base=r.system.euler.baseGeometry[wake.body];
  let previous=base.center,arc=0,count=0;
  for(let k=0;k<wake.ids.length;k++){
   const id=wake.ids[k],s=v.layers.states[id],station=r.system.bl.stations[id];
   if(k){const a=v.outer.nodes[wake.body][station.i].at(-1),b=v.outer.nodes[wake.body+1][station.i][0],p={x:(a.x+b.x)/2,y:(a.y+b.y)/2};arc+=Math.hypot(p.x-previous.x,p.y-previous.y);previous=p;}
   const gap=independentGap(base,arc,L);close(s.wakeGap,gap,2e-15,'independent expanded cubic');if(arc>0&&gap>0)count++;
   const u=k/factor,left=Math.min(original.ids.length-2,Math.floor(u)),t=u-left;
   const a=old.layers.states[original.ids[left]],b=old.layers.states[original.ids[left+1]],mix=(a,b)=>t===0?a:t===1?b:a+t*(b-a);
   const fluid=mix(a.deltaStar-a.wakeGap,b.deltaStar-b.wakeGap);
   close(s.deltaStar-s.wakeGap,fluid,2e-15,'fluid displacement');
   for(const key of ['theta','ue','aux'])close(s[key],mix(a[key],b[key]),2e-15,key);
   maxNaiveError=Math.max(maxNaiveError,Math.abs(s.deltaStar-mix(a.deltaStar,b.deltaStar)));
   if(Number.isInteger(u)) maximumRetainedTotalDeltaChange=Math.max(maximumRetainedTotalDeltaChange,Math.abs(s.deltaStar-mix(a.deltaStar,b.deltaStar)));
   assert.ok(s.deltaStar-s.wakeGap>s.theta); assert.ok(r.system.bl.kernel.station(s,'wake').rho>0);
  }
  active.push(count);
 }
 assert.ok(active.every(n=>n>=3),'Both distinct tails must be simultaneously active.');
 if(factor>1)assert.ok(maxNaiveError>1e-5,'A naive total-delta interpolation must measurably fail this fixture.');
 const replay=createCoupledStreamtubeBody(r.input,{...r.options,initialEuler:r.initialEuler,initialBL:r.initialBL});
 assert.deepEqual(replay.initial,r.system.initial);assert.deepEqual(replay.evaluate(replay.initial).residual,v.residual);
 assert.deepEqual(replay.bl.snapshotActive(),r.system.bl.snapshotActive());
 return{active,maxNaiveError,maximumRetainedTotalDeltaChange,quality:r.diagnostics.quality,families:v.families};
}
for(const automatic of [false,true])test(`two active cubic tails: ${automatic?'automatic':'fixed'} full streamwise/normal hybrid refinement restores final-arc gaps once`,t=>{
 const f=fixture({automatic}),before=plain({input:f.input,state:f.system.initial}),r=refine(f.input,f.system,{streamwiseFactor:2,normalFactor:2});
 const result=check(f,r,2);
 if(automatic) {
  assert.deepEqual(f.system.bl.snapshotActive(),f.system.bl.surfaces.map(s=>s.ids.length-1));
  assert.deepEqual(r.system.bl.snapshotActive(),r.system.bl.surfaces.map(s=>s.ids.length-1));
  assert.deepEqual(r.options.transitionState,r.system.bl.snapshotActive());
 }
 assert.deepEqual(plain({input:f.input,state:f.system.initial}),before);
 assert.equal(r.system.bl.transitionMode,automatic?'automatic':'fixed-trip');assert.deepEqual(r.system.bl.trips,f.system.bl.trips);
 t.diagnostic(JSON.stringify(result));
});
test('normal-only finite refinement preserves all original wake primitives and geometry',()=>{
 const f=fixture(),r=refine(f.input,f.system,{streamwiseFactor:1,normalFactor:2});check(f,r,1);
 const old=f.system.evaluate(f.system.initial),v=r.system.evaluate(r.system.initial);
 for(let g=0;g<old.outer.nodes.length;g++)for(let i=0;i<old.outer.nodes[g].length;i++)for(let j=0;j<old.outer.nodes[g][i].length;j++){
  close(v.outer.nodes[g][i][2*j].x,old.outer.nodes[g][i][j].x);close(v.outer.nodes[g][i][2*j].y,old.outer.nodes[g][i][j].y);
 }
 assert.deepEqual(r.system.bl.snapshotActive(),f.system.bl.snapshotActive());
});

test('the finite-base fluid-gap correction also applies to the unchanged legacy BL thermodynamic model',()=>{
 const f=twoActiveFiniteBaseWakes(),r=refine(f.input,f.system,{streamwiseFactor:2,normalFactor:2});
 check(f,r,2);assert.equal(r.options.blThermodynamics,undefined);assert.equal(r.diagnostics.physicalDensityTransfer,undefined);
 assert.ok(r.diagnostics.finiteBaseTransfer);assert.deepEqual(r.system.conditions,f.system.conditions);
});
