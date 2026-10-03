import test from 'node:test';
import assert from 'node:assert/strict';
import { channelMesh } from '../src/euler/mesh.js';
import { affine } from '../src/potential/reconstruction.js';
import { createPotentialSystem,solvePotential } from '../src/potential/system.js';
import { sparseDense } from '../src/numerics/sparse.js';
import { annularVortex } from './oracles/euler.js';

const grid=()=>channelMesh({nx:6,ny:4,map:({x,y,i,j,nx,ny})=>({x:x+.06*Math.sin(Math.PI*i/nx)*Math.sin(Math.PI*j/ny),y:y+.05*Math.sin(2*Math.PI*i/nx)*Math.sin(Math.PI*j/ny)})});
test('full potential reproduces uniform subcritical flow and a quadratic harmonic potential on a skew grid',()=>{
  const mesh=grid(),uniform=createPotentialSystem(mesh,{mach:.5});
  assert.ok(uniform.evaluate(new Float64Array(uniform.n)).diagnostics.residual<1e-13);
  const potential=({x,y})=>.1*(x*x-y*y)+.2*x*y;
  const velocity=({x,y})=>[1+.2*x+.2*y,-.2*y+.2*x];
  const system=createPotentialSystem(mesh,{boundaryPotential:p=>affine(potential(p)),boundaryNormalVelocity:f=>{const [u,v]=velocity(f);return affine(u*f.nx+v*f.ny);}});
  const exact=Float64Array.from(mesh.cells,potential),a=system.evaluate(exact);
  assert.ok(a.diagnostics.residual<1e-13);assert.ok(a.diagnostics.sharedFluxCancellation<1e-14);
  const r=solvePotential(system);assert.equal(r.converged,true);assert.equal(r.history.length,2);
  assert.ok(Math.max(...r.x.map((v,i)=>Math.abs(v-exact[i])))<1e-11);
});

test('full potential density and circulation Jacobian matches whole-residual differences and CSR',()=>{
  const mesh=grid(),count=mesh.cells.length;
  const options={mach:.4,baseVelocity:()=>[affine(1,[[count,.1]]),affine(0,[[count,.2]])],constraints:()=>[affine(-.01,[[count,1],[0,.1]])]};
  const dense=createPotentialSystem(mesh,options),sparse=createPotentialSystem(mesh,{...options,sparse:true});
  const x=Float64Array.from({length:dense.n},(_,i)=>.002*Math.sin(.37*i)),d=x.map((_,i)=>Math.cos(.71*i)),h=1e-6;
  const a=dense.evaluate(x,{jacobian:true}),b=sparseDense(sparse.evaluate(x,{jacobian:true}).jacobian);
  assert.deepEqual(a.jacobian,b);
  const plus=dense.evaluate(x.map((v,i)=>v+h*d[i])).residual,minus=dense.evaluate(x.map((v,i)=>v-h*d[i])).residual;
  for(let row=0;row<dense.n;row++){
    let jv=0;for(let col=0;col<dense.n;col++)jv+=a.jacobian[row*dense.n+col]*d[col];
    assert.ok(Math.abs(jv-(plus[row]-minus[row])/(2*h))<1e-7);
  }
  assert.ok(a.diagnostics.sharedFluxCancellation<1e-14);
});

test('the incompressible initialization shares operators and exactly preserves independent residuals, derivatives and forces',()=>{
  const mesh=grid(),count=mesh.cells.length;
  for(const sparse of [false,true])for(const order of [1,2]){
    const options={mach:.4,sparse,order,
      baseVelocity:()=>[affine(1,[[count,.1]]),affine(0,[[count,.2]])],
      boundaryPotential:p=>affine(.002*p.x*p.y),
      constraints:()=>[affine(-.01,[[count,1],[0,.1]])]};
    const physical=createPotentialSystem(mesh,options),reference=createPotentialSystem(mesh,{...options,mach:0});
    const x=Float64Array.from({length:physical.n},(_,i)=>.002*Math.sin(.37*i));
    const before=physical.evaluate(x,{jacobian:true}),guess=physical.incompressibleInitialization();
    for(const key of ['mesh','reconstruction','faces','cellVelocities','constraints'])assert.equal(guess[key],physical[key],key);
    assert.equal(guess.conditions.mach,0);assert.equal(physical.conditions.mach,.4);
    assert.deepEqual(guess.evaluate(x,{jacobian:true}),reference.evaluate(x,{jacobian:true}));
    assert.deepEqual(guess.preconditioner(x),reference.preconditioner(x));
    assert.deepEqual(physical.evaluate(x,{jacobian:true}),before);
    const fast=x.slice();fast[count]=20;
    assert.equal(guess.admissible(fast),reference.admissible(fast));
    assert.equal(guess.admissible(fast),true);assert.equal(physical.admissible(fast),false);
  }
});

test('compressible potential converges to the independent irrotational vortex',()=>{
  const errors=[];
  for(const [nx,ny] of [[6,2],[12,4],[24,8]]){
    const {mesh,exact}=annularVortex(nx,ny),potential=({x,y})=>1.5*Math.atan2(y,x);
    const system=createPotentialSystem(mesh,{mach:.3,baseVelocity:()=>[affine(),affine()],boundaryPotential:p=>affine(potential(p))});
    const r=solvePotential(system,{initial:mesh.cells.map(potential)});
    assert.equal(r.converged,true,r.reason);assert.ok(r.diagnostics.residual<1e-10);assert.ok(r.diagnostics.relativeMassImbalance<1e-9);
    let error=0,area=0;
    r.states.forEach((s,i)=>{const c=mesh.cells[i],e=exact(c);area+=c.area;error+=c.area*Math.hypot(s.u-e.u,s.v-e.v);});
    errors.push(error/area);
  }
  for(let i=1;i<errors.length;i++)assert.ok(errors[i]<.4*errors[i-1],String(errors));
  assert.ok(errors.at(-1)<.001,String(errors));
});

test('analytic circulation remains divergence-free on skew faces with a nearby vortex and retains its exact Jacobian',()=>{
  const mesh=grid();for(const f of mesh.faces)if(f.boundary)f.boundary.type='farfield';
  const count=mesh.cells.length,center={x:-.08,y:.43};
  const baseVelocity=p=>{
    const x=p.x-center.x,y=p.y-center.y,c=1/(2*Math.PI*(x*x+y*y));
    return[affine(1,[[count,-y*c]]),affine(0,[[count,x*c]])];
  };
  const baseFluxIntegral=f=>{
    const a=mesh.vertices[f.a],b=mesh.vertices[f.b];
    const ra=Math.hypot(a.x-center.x,a.y-center.y),rb=Math.hypot(b.x-center.x,b.y-center.y);
    return affine(f.nx*f.length,[[count,Math.log(ra/rb)/(2*Math.PI)]]);
  };
  for(const fluxQuadrature of [1,2,3]){
    const system=createPotentialSystem(mesh,{baseVelocity,baseFluxIntegral,fluxQuadrature,constraints:()=>[affine(-.3,[[count,1]])]});
    const x=new Float64Array(system.n);x[count]=.3;
    const actual=system.evaluate(x,{jacobian:true});
    assert.ok(actual.diagnostics.residual<2e-14);
    for(let row=0;row<count;row++)assert.ok(Math.abs(actual.jacobian[row*system.n+count])<2e-14);
    const solved=solvePotential(system);
    assert.equal(solved.converged,true);assert.ok(Math.max(...solved.x.map((v,i)=>Math.abs(v-x[i])))<1e-12);
  }
});

test('face pressure and its moment integrate the exact Bernoulli polynomial instead of a midpoint sample',()=>{
  const mesh=grid(),phi=p=>.1*(p.x*p.x-p.y*p.y);
  const velocity=p=>[1+.2*p.x,-.2*p.y];
  const s=createPotentialSystem(mesh,{boundaryPotential:p=>affine(phi(p)),boundaryNormalVelocity:f=>{const [u,v]=velocity(f);return affine(u*f.nx+v*f.ny);}});
  const r=s.evaluate(Float64Array.from(mesh.cells,phi));
  mesh.faces.forEach((f,i)=>{
    const a=mesh.vertices[f.a],b=mesh.vertices[f.b],[u,v]=velocity(f),va=velocity(a),vb=velocity(b),du=vb[0]-va[0],dv=vb[1]-va[1];
    const c0=1-u*u-v*v,c1=-2*(u*du+v*dv),c2=-du*du-dv*dv;
    const pressure=f.length*(c0+c2/12);
    const moment=f.length*((f.y*f.nx-f.x*f.ny)*(c0+c2/12)+((b.y-a.y)*f.nx-(b.x-a.x)*f.ny)*c1/12);
    assert.ok(Math.abs(r.pressureIntegrals[i]-pressure)<2e-14);
    assert.ok(Math.abs(r.pressureMomentIntegrals[i]-moment)<2e-14);
  });
});
