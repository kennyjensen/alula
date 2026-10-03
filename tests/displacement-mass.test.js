import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4,transform } from '../src/geometry/airfoil.js';
import { createDisplacementOperator } from '../src/inviscid/displacement.js';
import { multielementMesh } from '../src/euler/multielement-mesh.js';
import { displacementMass,meshDisplacementSources } from '../src/viscous/mass.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import {nodalMassSlopes,integrateMassDerivative} from '../src/viscous/mass-interpolation.js';

test('continuous body sources recover a nonuniform quadratic mass profile and preserve each subinterval mass exactly',()=>{
  const s=[0,.001,.02,.11,.3,.55,1],mass=t=>.002*(1+3*t+2*t*t),derivative=t=>.002*(3+4*t);
  const slopes=nodalMassSlopes(s),apply=w=>[...w].reduce((v,[i,c])=>v+c*mass(s[i]),0);
  slopes.forEach((w,i)=>assert.ok(Math.abs(apply(w)-derivative(s[i]))<2e-14));
  for(let i=0;i<s.length-1;i++)for(const [lo,hi]of [[0,.13],[.13,.7],[.7,1],[0,1]]){
    const h=s[i+1]-s[i],weights=integrateMassDerivative(i,h,slopes[i],slopes[i+1],lo,hi);
    assert.ok(Math.abs(apply(weights)-(mass(s[i]+hi*h)-mass(s[i]+lo*h)))<1e-16);
  }
});

test('absolute per-element wake extents keep small-element wakes in the common downstream region',()=>{
  const elements=[{points:naca4('0012',40)},{points:transform(naca4('0012',40),{chord:.3,x:1.05,y:-.2})}],target=3.35;
  const lengths=elements.map(e=>target-e.points[0].x);
  for(const wakeCount of [4,24]){
    const o=createDisplacementOperator({elements,wakeLengths:lengths,wakeCount});
    o.wakes.forEach((w,e)=>{assert.ok(Math.abs(w.s.at(-1)-lengths[e])<1e-12);assert.ok(Math.abs(w.points.at(-1).x-target)<1e-12);});
  }
  assert.throws(()=>createDisplacementOperator({elements,wakeLengths:[2]}),/one positive/);
});

test('compressible mass deficit uses the BL edge density and its complete derivatives',()=>{
  for(const mach of [0,.4])for(const q of [-1.2,.4,1.1]){
    const delta=.002,state=displacementMass(q,delta,{mach}),h=1e-6;
    const kernel=createIntegralKernel({mach}),bl=kernel.station({s:.2,theta:.0008,deltaStar:delta,ue:Math.abs(q),aux:0});
    assert.ok(Math.abs(state.value-bl.rho*q*delta)<1e-15);
    assert.ok(Math.abs(state.velocityDerivative-(displacementMass(q+h,delta,{mach}).value-displacementMass(q-h,delta,{mach}).value)/(2*h))<1e-12);
    assert.ok(Math.abs(state.thicknessDerivative-(displacementMass(q,delta+h,{mach}).value-displacementMass(q,delta-h,{mach}).value)/(2*h))<1e-12);
  }
});

for(const bodyMassInterpolation of ['linear','hermite'])test(`${bodyMassInterpolation} subdivided body and curved-wake sources preserve each element budget and their complete Jacobian`,()=>{
  const elements=[{points:naca4('0012',40)},{points:transform(naca4('0012',40),{chord:.3,x:1.05,y:-.2})}];
  const outer=createDisplacementOperator({elements,alpha:3,wakeCount:12,wakeLength:1,wakeInitialization:'inviscid'});
  const x=new Float64Array(4*outer.total),offset=0,thicknessScale=.001,mach=.3;
  for(let i=0;i<outer.total;i++){x[4*i+2]=2+.4*Math.sin(.2*i);x[4*i+3]=outer.q0[i];}
  const mesh=multielementMesh(elements.map(e=>e.points),{rows:2,padding:2,subdivisions:2,wakePaths:outer.wakes.map(w=>w.points)});
  const {bodySources,wakeSources}=meshDisplacementSources(mesh,outer,{offset,thicknessScale,mach,bodyMassInterpolation});
  const mass=i=>displacementMass(x[4*i+3],x[4*i+2]*thicknessScale,{mach}).value;
  outer.bodies.forEach((b,e)=>{
    let sum=0;for(const [i,s] of bodySources)if(mesh.faces[i].boundary.element===e)sum+=s.evaluate(x).value;
    assert.ok(Math.abs(sum-(mass(b.end)-mass(b.start)))<1e-13);
    sum=0;for(const cut of mesh.cuts.filter(c=>c.type==='wake-cut'&&c.element===e))sum+=wakeSources.get(cut.face).evaluate(x).value;
    const w=outer.wakes[e];assert.ok(Math.abs(sum-(mass(w.end)-mass(w.start)))<1e-13);
  });
  const d=x.map((_,i)=>.1*Math.sin(.33*i)),h=1e-5;
  for(const source of [...bodySources.values(),...wakeSources.values()]){
    const a=source.evaluate(x),fd=(source.evaluate(x.map((v,i)=>v+h*d[i])).value-source.evaluate(x.map((v,i)=>v-h*d[i])).value)/(2*h);
    let derivative=0;for(const [col,v] of a.derivatives)derivative+=v*d[col];
    assert.ok(Math.abs(fd-derivative)<1e-12);
  }
});
