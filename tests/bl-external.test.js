import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4,transform } from '../src/geometry/airfoil.js';
import { createBoundaryLayerAssembly,createCoupledAssembly } from '../src/viscous/assembly.js';

const input={elements:[{points:naca4('0012',40)},{points:transform(naca4('0012',40),{chord:.3,x:1.05,y:-.15})}],alpha:2,trips:[.05,.05],wakeCount:12};
for(const wakeInitialization of ['straight','inviscid'])test(`external BL block retains all native equations and ${wakeInitialization} wake geometry without unused panel influence`,()=>{
  const a=createCoupledAssembly({...input,wakeInitialization}),b=createBoundaryLayerAssembly({...input,wakeInitialization}),x=a.initial;
  assert.deepEqual(b.initial,a.initial);
  a.updateActive(x);b.updateActive(x);
  const ra=a.residual(x),rb=b.residual(x),ja=a.jacobian(x),jb=b.jacobian(x),n=x.length;
  for(let row=0;row<n;row++){
    if(row%4===3){assert.equal(rb[row],0);continue;}
    assert.equal(ra[row],rb[row]);
    for(let col=0;col<n;col++)assert.equal(ja[row*n+col],jb[row*n+col]);
  }
  assert.equal(b.externalVelocity,true);assert.equal(b.linearSolve,undefined);
  assert.equal(b.outer.influence,undefined);assert.equal(b.outer.evaluate,undefined);
  assert.deepEqual(b.outer.sourceMatrix,a.outer.sourceMatrix);
  assert.deepEqual(b.outer.wakes,a.outer.wakes);
});

test('finite-Mach multielement BL derivatives include physical similarity velocity and moving stagnation distance',()=>{
  const system=createBoundaryLayerAssembly({...input,mach:.3}),x=system.initial.slice();system.updateActive(x);
  const j=system.jacobian(x),d=x.map((v,i)=>Math.sin(.31*i+.2)*(i%4===3?.1:i%4===0?.001:.01)),h=1e-5;
  const plus=system.residual(x.map((v,i)=>v+h*d[i])),minus=system.residual(x.map((v,i)=>v-h*d[i]));
  let largest=0;
  for(let row=0;row<x.length;row++)if(row%4!==3){
    let value=0;for(let col=0;col<x.length;col++)value+=j[row*x.length+col]*d[col];
    const fd=(plus[row]-minus[row])/(2*h),error=Math.abs(value-fd)/Math.max(1,Math.abs(fd));largest=Math.max(largest,error);
  }
  assert.ok(largest<2e-6,String(largest));
});
