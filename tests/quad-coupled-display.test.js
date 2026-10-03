import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';
import {quadCoupledResultForDisplay} from '../src/ui/quad-coupled-result.js';
import {createQuadMeshProgress} from '../src/ui/quad-mesh-progress.js';
const saved=JSON.parse(fs.readFileSync(new URL('../docs/current-coupled-startup-default-browser.json',import.meta.url)));
const raw=saved.result,input=saved.input;
const near=(a,b,tol=1e-12)=>assert.ok(Math.abs(a-b)<=tol,`${a} != ${b}`);
test('coupled display preserves solved state, maps all bodies and matches physical wall/wake displacement',()=>{
 const before=JSON.stringify(raw),view=quadCoupledResultForDisplay(raw,input);
 assert.equal(JSON.stringify(raw),before);assert.deepEqual(view.boundaryLayer.surfaces.map(s=>[s.element,s.side]),[[0,'upper'],[0,'lower'],[1,'upper'],[1,'lower']]);
 assert.deepEqual(view.boundaryLayer.wakes.map(s=>s.element),[0,1]);
 assert.equal(view.boundaryLayer.stations.length,raw.boundaryLayer.stations.length);
 assert.equal(view.numericalBoundaryLayer,raw.boundaryLayer);assert.equal(view.restart,raw.restart);
 for(const surface of view.boundaryLayer.surfaces)for(const p of surface.stations){
  near(Math.hypot(p.displacement.x-p.x,p.displacement.y-p.y),p.deltaStar,2e-15);
  assert.ok([p.cf,p.h,p.hk,p.cp].every(Number.isFinite));
  const M2=raw.mach**2,gamma=1.4,T=1+.5*(gamma-1)*M2*(1-p.ue**2);
  near(p.cp,2/(gamma*M2)*(T**(gamma/(gamma-1))-1),2e-13);
 }
 for(const w of view.boundaryLayer.wakes)for(const p of w.stations.slice(1)){
  const b=p.body,i=p.i,right=Math.min(raw.flow.nodes[b].length-1,i+1);
  const a=raw.flow.nodes[b][i].at(-1),c=raw.flow.nodes[b+1][i][0];
  const center=k=>({x:.5*(raw.flow.nodes[b][k].at(-1).x+raw.flow.nodes[b+1][k][0].x),y:.5*(raw.flow.nodes[b][k].at(-1).y+raw.flow.nodes[b+1][k][0].y)});
  const lo=center(i-1),hi=center(right),dx=hi.x-lo.x,dy=hi.y-lo.y;
  near(((c.y-a.y)*dx-(c.x-a.x)*dy)/Math.hypot(dx,dy),p.deltaStar,1e-9);
 }
 assert.ok([view.cl,view.cd,view.cm].every(Number.isFinite));assert.ok(view.cd>0);assert.equal(view.physicalAcceptance,false);
 assert.equal(view.coefficientStatus,'research-unvalidated');assert.equal(view.coefficients.wakes.length,2);
 near(view.cd,view.coefficients.wakes.reduce((sum,w)=>sum+w.cd,0));
});
test('geometric scaling doubles displayed lengths and preserves dimensionless BL/pressure profiles',()=>{
 const scaled=structuredClone(raw),scaledInput=structuredClone(input),k=2;
 const scale=p=>({x:k*p.x,y:k*p.y});
 scaled.solverLength*=k;scaled.referenceChord*=k;scaledInput.referenceChord*=k;
 scaledInput.elements.forEach(e=>{e.points=e.points.map(scale);});
 scaled.solverInput.bodies.forEach(b=>{b.points=b.points.map(scale);});
 scaled.boundaryLayer.surfaces.forEach(s=>{s.tripParameter*=k;});
 for(const key of ['nodes','undisplacedNodes'])scaled.flow[key]=scaled.flow[key].map(g=>g.map(row=>row.map(scale)));
 scaled.mesh.vertices=scaled.mesh.vertices.map(scale);
 const a=quadCoupledResultForDisplay(raw,input),b=quadCoupledResultForDisplay(scaled,scaledInput);
 for(const key of ['cl','cd','cm'])near(a[key],b[key]);
 for(let i=0;i<a.boundaryLayer.stations.length;i++){
  const p=a.boundaryLayer.stations[i],q=b.boundaryLayer.stations[i];
  for(const key of ['theta','deltaStar','s','x','y'])near(q[key],k*p[key]);
  for(const key of ['cf','h','ue'])near(q[key],p[key]);
 }
 for(let i=0;i<a.boundaryLayer.surfaces.length;i++){
  const p=a.boundaryLayer.surfaces[i].transitionPoint,q=b.boundaryLayer.surfaces[i].transitionPoint;near(q.x,2*p.x);near(q.y,2*p.y);
 }
});
test('mesh progress measures real movement and resets comparisons between coupled attempts',()=>{
 const p=createQuadMeshProgress(2),line={group:0,tube:0,points:[{x:0,y:0},{x:1,y:0}],speedRatios:[1]};
 const mesh=(dx,q)=>({vertices:[{x:dx,y:0},{x:1+dx,y:0}],cells:[],initialization:{},iteration:{iteration:dx?1:0},flow:{lines:[{...line,speedRatios:[q]}],iteration:dx?1:0,residual:0}});
 p.stage({stage:'coupled',startupAttempt:1});const a=p.mesh(mesh(0,1)),b=p.mesh(mesh(.1,1.2));
 near(a.iteration.maximumNodeMovement,0);near(b.iteration.maximumNodeMovement,.1);near(b.flow.maximumSpeedChangeFromPrevious,.2);
 p.stage({stage:'coupled',startupAttempt:2});const c=p.mesh(mesh(.5,.8));
 near(c.iteration.maximumNodeMovement,0);near(c.flow.maximumSpeedChangeFromPrevious,0);assert.equal(c.iteration.startupAttempt,2);
});

test('restored mesh iteration Ncrit takes precedence over the failed stage label',()=>{
 const p=createQuadMeshProgress(1);
 p.stage({stage:'coupled',actualNcrit:9,targetNcrit:9});
 const mesh={vertices:[{x:0,y:0},{x:1,y:0}],cells:[],initialization:{},
  iteration:{iteration:4,actualNcrit:8,targetNcrit:9},
  flow:{lines:[{group:0,tube:0,points:[{x:0,y:0},{x:1,y:0}],speedRatios:[1]}]}};
 const before=structuredClone(mesh),restored=p.mesh(mesh);
 assert.equal(restored.actualNcrit,8);assert.equal(restored.iteration.actualNcrit,8);
 assert.equal(restored.targetNcrit,9);assert.deepEqual(mesh,before);
});
