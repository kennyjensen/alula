import fs from 'node:fs';
import { createHash } from 'node:crypto';
import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { initializeStreamtubeDensities } from '../src/euler/streamtube-initial-state.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { finiteBaseBodyFixture } from './fixtures/finite-base-body.js';

import { refineStreamtubeBody as refine } from '../src/euler/streamtube-refinement.js';
const baselineBytes=fs.readFileSync(new URL('../docs/rae2822/transonic-refinement-audit/finite-base-draft/production-before/streamtube-refinement.js.txt',import.meta.url));
assert.equal(createHash('sha256').update(baselineBytes).digest('hex'),'8d727a56ff02d6e5556ad9401f978ac8319a05cc1e3cd945b31105c272b415fd');
const baselineText=baselineBytes.toString().replace(/from '([^']+)'/g,(_,p)=>`from '${new URL(p,new URL('../src/euler/streamtube-refinement.js',import.meta.url))}'`);
const { refineStreamtubeBody:currentRefine }=await import('data:text/javascript;base64,'+Buffer.from(baselineText).toString('base64'));

const close=(a,b,t=2e-11)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
const entropy=s=>Math.log(s.enthalpy)/.4-Math.log(s.rho);
const controls={streamwiseFactor:1,normalSubdivisions:[[2],[2]]};
function fixture(mach=1.2,jump=false,mode='hybrid'){
 const input={...intrinsicBodyFixture({bodySegments:4,tubes:1,mach:.5}),streamwiseMode:'hybrid',hybrid:{epsilonP:1e-5},upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}};
 if(mode==='momentum'){input.streamwiseMode='momentum';delete input.hybrid;}
 const source=createStreamtubeBodySystem(input),x=initializeStreamtubeDensities(source,source.initial),before=source.evaluate(x);
 const first=input.bodies[0].trailingIndex+1,last=source.layout.nx-2;
 const q=Math.sqrt(.4*mach*mach*source.conditions.h0/(1+.2*mach*mach));
 for(let i=first;i<=last;i++)for(let g=0;g<2;g++){
  const s=before.sections[i][g][0],speed=jump&&i===last?q*.8:q;
  x[source.layout.densityIndex(i,g,0)]=Math.log(s.rho*s.q/speed);
 }
 assert.ok(source.admissible(x));return{input,source,x,first,last};
}
for(const mach of [.6,1.2])test(`complete Euler hybrid normal refinement preserves the ${mach<1?'sub':'super'}sonic uniform interior`,()=>{
 const f=fixture(mach),before=f.x.slice(),r=refine(f.input,f.source,{...controls,initial:f.x}),v=r.system.evaluate(r.initial),old=f.source.evaluate(f.x);
 assert.deepEqual(r.system.conditions,f.source.conditions);assert.ok(r.diagnostics.quality.valid);assert.deepEqual(f.x,before);
 for(let i=f.first;i<=f.last;i++)for(let g=0;g<2;g++)for(let j=0;j<2;j++){
  assert.equal(r.initial[r.system.layout.densityIndex(i,g,j)],f.x[f.source.layout.densityIndex(i,g,0)]);
  close(v.sections[i][g][j].q,old.sections[i][g][0].q);assert.equal(v.sections[i][g][j].machSquared>1,mach>1);
 }
 const replay=createStreamtubeBodySystem(r.input),x=replay.adoptGeometry(r.initialEuler.x,r.initialEuler.nodes);
 assert.deepEqual(replay.evaluate(x).residual,v.residual);
});
test('complete hybrid Euler normal refinement preserves collocated physical entropy rather than returning a common isentrope',()=>{
 const f=fixture(1.2,true),r=refine(f.input,f.source,{...controls,initial:f.x,streamwiseFactor:1,normalSubdivisions:[[3],[3]]}),v=r.system.evaluate(r.initial),old=f.source.evaluate(f.x);
 assert.ok(r.diagnostics.quality.valid);assert.ok(r.diagnostics.physicalDensityTransfer.collocatedValues>0);
 let span=0;
 for(let i=f.first;i<=f.last;i++)for(let g=0;g<2;g++)for(let j=0;j<3;j++){
  const s=v.sections[i][g][j],p=old.sections[i][g][0];
  assert.equal(r.initial[r.system.layout.densityIndex(i,g,j)],f.x[f.source.layout.densityIndex(i,g,0)]);
  close(entropy(s),entropy(p));span=Math.max(span,Math.abs(entropy(p)-entropy(old.sections[0][g][0])));
 }
 assert.ok(span>.01);assert.ok(v.residual.every(Number.isFinite));
});
test('finite-base hybrid Euler normal refinement retains constant physical base gaps without isentropic inversion',()=>{
 const input={...finiteBaseBodyFixture(),streamwiseMode:'hybrid',hybrid:{epsilonP:1e-5},upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}};
 const source=createStreamtubeBodySystem(input),r=refine(input,source,{streamwiseFactor:1,normalSubdivisions:source.layout.tubes.map(n=>Array(n).fill(2))});
 assert.deepEqual(r.system.displacement,source.displacement);assert.ok(r.diagnostics.physicalDensityTransfer);assert.ok(r.diagnostics.quality.valid);
 const old=source.decode(source.initial),v=r.system.decode(r.initial);
 for(let g=0;g<v.nodes.length;g++)for(let i=0;i<v.nodes[g].length;i++){
  assert.deepEqual(v.nodes[g][i][0],old.nodes[g][i][0]);assert.deepEqual(v.nodes[g][i].at(-1),old.nodes[g][i].at(-1));
 }
 assert.throws(()=>refine(input,source,{streamwiseFactor:2}),/normal subdivisions only/);
});
test('all omitted-mode and geometry-only Euler refinement outputs remain byte-exact in the draft',()=>{
 for(const mode of ['isentropic','momentum','hybrid'])for(const initializeFlow of [false,true]){
  if(mode==='hybrid'&&initializeFlow)continue;
  const f=fixture();if(mode!=='hybrid'){f.input.streamwiseMode=mode;delete f.input.hybrid;delete f.input.upwind;f.source=createStreamtubeBodySystem(f.input);f.x=f.source.initial;}
  const c={...controls,initial:f.x,initializeFlow},a=currentRefine(f.input,f.source,c),b=refine(f.input,f.source,c);
  for(const k of ['input','initial','initialEuler','diagnostics'])assert.deepEqual(b[k],a[k]);
 }
});

test('full streamwise Euler refinement rejects a negative interface pressure even when interpolated gas remains thermal',()=>{
 const f=fixture(1.2,true);
 assert.throws(()=>refine(f.input,f.source,{...controls,initial:f.x,streamwiseFactor:3}),/Nonpositive or nonfinite streamline interface pressure/);
});

test('supersonic momentum Euler with explicit upwinding uses the same physical-density prolongation',()=>{
 const f=fixture(1.2,true,'momentum'),r=refine(f.input,f.source,{...controls,initial:f.x}),v=r.system.evaluate(r.initial),old=f.source.evaluate(f.x);
 assert.deepEqual(r.system.conditions,f.source.conditions);assert.ok(r.diagnostics.physicalDensityTransfer);
 for(let i=f.first;i<=f.last;i++)for(let g=0;g<2;g++)for(let j=0;j<2;j++){
  assert.equal(r.initial[r.system.layout.densityIndex(i,g,j)],f.x[f.source.layout.densityIndex(i,g,0)]);
  close(entropy(v.sections[i][g][j]),entropy(old.sections[i][g][0]));
  assert.equal(v.sections[i][g][j].machSquared>1,old.sections[i][g][0].machSquared>1);
 }
});
