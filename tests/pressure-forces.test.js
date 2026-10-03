import test from 'node:test';
import assert from 'node:assert/strict';
import { pressureForces } from '../src/inviscid/pressure-forces.js';
import { naca4 } from '../src/geometry/airfoil.js';
import { clcalc } from '../src/viscous/xfoil/xfoil.js';

test('pressure quadrature satisfies uniform/linear pressure area and moment identities',()=>{
  const points=[{x:0,y:0},{x:2,y:0},{x:2,y:1},{x:0,y:1},{x:0,y:0}];
  for(const [b,c] of [[0,0],[2,0],[0,3],[2,-3]]){
    const r=pressureForces(points,points.map(p=>1+b*p.x+c*p.y),{momentOrigin:{x:0,y:0}});
    assert.ok(Math.abs(r.cx+2*b)<1e-13);assert.ok(Math.abs(r.cy+2*c)<1e-13);
    assert.ok(Math.abs(r.cm-2*(c-.5*b))<1e-13);
  }
});

test('physical-Cp quadrature reproduces native CLCALC at Mach zero',()=>{
  const points=naca4('2412',160),q=points.map((p,i)=>Math.sin(i*.13)*.2+1),alpha=4;
  const a=pressureForces(points,q.map(v=>1-v*v),{alpha}),b=clcalc(points.length,points.map(p=>p.x),points.map(p=>p.y),q,null,alpha*Math.PI/180,0,1,.25,0);
  assert.ok(Math.abs(a.cl-b.cl)<1e-14);assert.ok(Math.abs(a.cm-b.cm)<1e-14);assert.ok(Math.abs(a.cd-b.cdp)<1e-14);
});
