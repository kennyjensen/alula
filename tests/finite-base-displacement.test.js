import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createDisplacementOperator, polynomialPotential } from '../src/inviscid/displacement.js';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';
import { streamfunctionAt } from '../src/inviscid/streamfunction.js';
import { transform, naca4 } from '../src/geometry/airfoil.js';
import { finiteBaseBodyFixture } from './fixtures/finite-base-body.js';

const element = () => {
  const b=finiteBaseBodyFixture().bodies[0];
  return {points:b.points,trailingEdge:b.trailingEdge};
};
const max=values=>Math.max(...values.map(Math.abs));
const input={elements:[element()],alpha:2,wakeCount:12,wakeInitialization:'inviscid'};
const outer=createDisplacementOperator(input);

test('finite displacement operator recovers the checked inviscid field without altering source geometry',()=>{
  const reference=solveInviscid(input),zero=new Float64Array(outer.total),field=outer.velocityField(zero);
  assert.deepEqual(outer.bodies[0].solidPoints,input.elements[0].points);
  assert.equal(outer.bodies[0].points.length,41);
  assert.equal(outer.bodies[0].basePanels.length,2);
  assert.ok(max(Array.from(reference.field.gamma,(q,i)=>q-outer.q0[i]))<2e-12);
  for(const p of [{x:-.2,y:.2},{x:1.2,y:.1},{x:.5,y:-.2}]){
    const a=field(p),b=velocityAt(p,reference.field);
    assert.ok(Math.hypot(a.u-b.u,a.v-b.v)<2e-12);
  }
});

test('coupled finite-base strengths satisfy body streamfunction constants and the Kutta relation',()=>{
  const mass=Float64Array.from({length:outer.total},(_,i)=>.0002*Math.sin(i*.7));
  const q=outer.evaluate(mass),sigma=outer.sources(mass),nc=sigma.length;
  const basePanels=outer.bodies.flatMap(b=>b.basePanels??[]).map(p=>({...p,
    sourceStrength:p.sourceCoefficient*(q[p.upperNode]-q[p.lowerNode]),
    vortexStrength:p.vortexCoefficient*(q[p.upperNode]-q[p.lowerNode])}));
  const field={...outer.field,gamma:q.slice(0,outer.surfaceCount),basePanels};
  for(const [e,body] of outer.bodies.entries()){
    const psi=body.points.map(point=>{
      let value=streamfunctionAt(point,field);
      for(const [j,p] of outer.sourcePanels.entries()){
        const dx=body.centroid.x-p.a.x,dy=body.centroid.y-p.a.y;
        const chart=j<outer.panels.length&&p.element===e?Math.PI/2:Math.atan2(-dx*p.ty+dy*p.tx,dx*p.tx+dy*p.ty);
        const basis=polynomialPotential(point,p,chart);
        for(let k=0;k<=p.degree;k++)value+=sigma[3*j+k]*basis[k].imag;
      }
      return value;
    });
    assert.ok(max(psi.map(v=>v-psi[0]))<2e-12,`psi spread ${max(psi.map(v=>v-psi[0]))}`);
    assert.ok(Math.abs(q[body.start]+q[body.end])<2e-12);
  }
  assert.equal(nc,3*outer.sourcePanels.length);
});

test('integrated base and displacement source budget equals wake exit flux plus the unresolved TE matching defect',()=>{
  const mass=Float64Array.from({length:outer.total},(_,i)=>.0002*Math.sin(i*.7));
  const q=outer.evaluate(mass),sigma=outer.sources(mass);
  // Two-point Gaussian integration is exact for each quadratic source.
  const nodes=[.5-.5/Math.sqrt(3),.5+.5/Math.sqrt(3)];
  const sourceFlux=outer.sourcePanels.reduce((sum,p,j)=>sum+.5*p.length*nodes.reduce((s,t)=>
    s+sigma[3*j]+t*sigma[3*j+1]+t*t*sigma[3*j+2],0),0);
  let baseFlux=0,expected=0;
  for(const w of outer.wakes){
    const b=w.body,g=b.baseGeometry.width;
    const base=-.5*g*(q[b.start]-q[b.end]);baseFlux+=base;
    expected+=mass[w.end]+base+mass[b.end]-mass[b.start]-mass[w.start];
    // Native PSWLIN first endpoint uses its own total-mass secant.
    assert.ok(Math.abs(sigma[3*w.sourceStart]-(mass[w.start+1]-mass[w.start])/(w.s[1]-w.s[0]))<2e-14);
    assert.equal(q[w.start],q[b.end]);
  }
  assert.ok(Math.abs(baseFlux+sourceFlux-expected)<2e-14);
  // Independent divergence check on the computed velocity field.
  const velocity=outer.velocityField(mass),radius=1e4,count=64;
  let circleFlux=0;
  for(let i=0;i<count;i++){
    const t=2*Math.PI*(i+.5)/count,c=Math.cos(t),s=Math.sin(t),v=velocity({x:radius*c,y:radius*s});
    circleFlux+=2*Math.PI*radius/count*((v.u-outer.field.u)*c+(v.v-outer.field.v)*s);
  }
  assert.ok(Math.abs(circleFlux-expected)<2e-8,`circle ${circleFlux}, budget ${expected}`);
});

test('sharp displacement geometry and complete influence remain identical to the archived implementation',async()=>{
  const path=new URL('../docs/nlr-panel-bl/before/src__inviscid__displacement.js.txt',import.meta.url);
  const source=fs.readFileSync(path,'utf8').replace(/from '([^']+)'/g,(_,p)=>`from '${new URL(p,new URL('../src/inviscid/displacement.js',import.meta.url)).href}'`);
  const archived=await import('data:text/javascript;base64,'+Buffer.from(source).toString('base64'));
  const input={elements:[{points:naca4('2412',40)}],alpha:2,wakeCount:12,wakeInitialization:'inviscid'};
  const a=createDisplacementOperator(input),b=archived.createDisplacementOperator(input);
  for(const key of ['bodies','wakes','panels','sourcePanels','sourceMatrix','q0','influence','field','diagnostics'])assert.deepEqual(a[key],b[key]);
});

test('finite-base displacement response is covariant under rotation, translation and common length scaling',()=>{
  const scale=2.3,angle=17,move={chord:scale,angle,x:.2,y:-.1};
  const transformed={...input,alpha:input.alpha+angle,elements:input.elements.map(e=>({...e,points:transform(e.points,move)}))};
  const b=createDisplacementOperator(transformed);
  const m=Float64Array.from({length:outer.total},(_,i)=>.0001*Math.sin(i*.3)),q=outer.evaluate(m),qb=b.evaluate(m.map(v=>v*scale));
  assert.ok(max(Array.from(q,(v,i)=>v-qb[i]))<2e-10);
  for(const [k,w] of outer.wakes.entries()){
    const points=transform(w.points,move);
    assert.ok(max(points.flatMap((p,i)=>[p.x-b.wakes[k].points[i].x,p.y-b.wakes[k].points[i].y]))<2e-11);
  }
});
