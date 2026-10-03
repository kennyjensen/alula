import test from 'node:test';
import assert from 'node:assert/strict';
import {naca4Standard,transform} from '../src/geometry/airfoil.js';
import {prepareAirfoilElement} from '../src/geometry/airfoil-element.js';
import {builtinAirfoils} from '../src/geometry/builtin-airfoils.js';
import {solveInviscid} from '../src/inviscid/linear-vortex.js';

test('standard NACA sections retain XFOIL polynomial and distinct trailing edges',()=>{
 for(const code of ['0012','2412','4415']){
  const p=naca4Standard(code,160),t=Number(code.slice(2))/100,m=Number(code[0])/100,c=Number(code[1])/10;
  assert.equal(p.length,161);
  assert.equal(p[0].x,1);assert.equal(p.at(-1).x,1);
  assert.ok(Math.abs(p[0].y-p.at(-1).y-.021*t)<1e-14);
  for(let i=0;i<=80;i++){
   const upper=p[80-i],lower=p[80+i],x=upper.x;
   const yt=(.29690*Math.sqrt(x)-.12600*x-.35160*x*x+.28430*x**3-.10150*x**4)*t/.20;
   const yc=m===0?0:x<c?m/c**2*(2*c*x-x*x):m/(1-c)**2*((1-2*c)+2*c*x-x*x);
   assert.ok(Math.abs(upper.y-(yc+yt))<1e-14);assert.ok(Math.abs(lower.y-(yc-yt))<1e-14);
  }
  const e=prepareAirfoilElement({points:p});
  assert.deepEqual(e.trailingEdge,{kind:'finite-base',upperIndex:0,lowerIndex:160});
  assert.deepEqual(e.sourcePoints,p);assert.deepEqual(e.points.at(-1),e.points[0]);
  assert.equal(e.points.length,162);
 }
});

test('all generated built-in assemblies preserve transformed finite bases and solve panel flow',()=>{
 for(const [preset,defs] of Object.entries(builtinAirfoils)){
  const elements=defs.map(d=>prepareAirfoilElement({name:d.name,points:transform(naca4Standard(d.code,80),{chord:d.chord,x:d.x,y:d.y,angle:-d.deflection})}));
  for(let i=0;i<elements.length;i++){
   const e=elements[i],a=e.points[0],b=e.points[e.trailingEdge.lowerIndex];
   assert.equal(e.trailingEdge.kind,'finite-base');
   assert.ok(Math.abs(Math.hypot(a.x-b.x,a.y-b.y)-.021*Number(defs[i].code.slice(2))/100*defs[i].chord)<1e-14);
  }
  const result=solveInviscid({elements,alpha:4,referenceChord:1});
  assert.ok(Number.isFinite(result.cl),preset);
 }
});
