import test from 'node:test';
import assert from 'node:assert/strict';
import { channelMesh } from '../src/euler/mesh.js';
import { linearReconstruction } from '../src/euler/reconstruction.js';
import { eulerResidual } from '../src/euler/residual.js';
import { eulerJacobian } from '../src/euler/jacobian.js';
import { freestream } from '../src/euler/gas.js';
import { encodeStates,decodeStates,solveEuler } from '../src/euler/solve.js';
import { annularVortex,weightedError } from './oracles/euler.js';

const grid=()=>channelMesh({nx:5,ny:3,map:({x,y,i,j,nx,ny})=>({x:x+.04*Math.sin(Math.PI*i/nx)*Math.sin(Math.PI*j/ny),y:y+.06*Math.sin(2*Math.PI*i/nx)*Math.sin(Math.PI*j/ny)})});
test('Euler reconstruction reproduces affine fields at skew-grid faces, including walls',()=>{
  const mesh=grid(),r=linearReconstruction(mesh);
  const exact=({x,y})=>({rho:1+.03*x-.02*y,u:.9+.1*x+.05*y,v:.2-.03*x+.02*y,p:10+.2*x-.1*y});
  const states=mesh.cells.map(exact);
  mesh.faces.forEach((f,i)=>{
    const a=mesh.vertices[f.a],b=mesh.vertices[f.b],expected=exact({x:.5*(a.x+b.x),y:.5*(a.y+b.y)});
    for(const s of [r.faces[i].left,r.faces[i].right])if(s){
      const actual=r.sample(s,states);for(const key of ['rho','u','v','p'])assert.ok(Math.abs(actual[key]-expected[key])<1e-12);
    }
  });
  const reference=freestream(),uniform=eulerResidual(mesh,mesh.cells.map(()=>reference),reference,{reconstruction:r});
  assert.ok(uniform.diagnostics.residual<5e-15);
});

test('reconstructed Euler Jacobian includes the extended stencil and preserves shared flux cancellation',()=>{
  const mesh=grid(),reconstruction=linearReconstruction(mesh),reference=freestream();
  const states=mesh.cells.map((_,i)=>({rho:1+.02*Math.sin(.77*i+.1),u:1+.04*Math.sin(1.13*i+.4),v:.02+.04*Math.cos(.93*i+.2),p:reference.p+.05*Math.cos(.53*i+.2)}));
  const x=encodeStates(states,reference),n=x.length,jac=eulerJacobian(mesh,states,reference,{reconstruction});
  const direction=x.map((_,i)=>Math.sin(.31*i+.2)),h=2e-6;
  const plus=eulerResidual(mesh,decodeStates(x.map((v,i)=>v+h*direction[i]),states.length,reference),reference,{reconstruction});
  const minus=eulerResidual(mesh,decodeStates(x.map((v,i)=>v-h*direction[i]),states.length,reference),reference,{reconstruction});
  assert.ok(plus.diagnostics.sharedFluxCancellation<1e-13);
  for(let row=0;row<n;row++){
    let value=0;for(let col=0;col<n;col++)value+=jac[row*n+col]*direction[col];
    assert.ok(Math.abs(value-(plus.residual[row]-minus.residual[row])/(2*h))<3e-7,`row ${row}`);
  }
});

test('linear Euler reconstruction gives second-order pressure convergence for the exact compressible vortex',()=>{
  const errors=[];
  for(const [nx,ny] of [[6,2],[12,4],[24,8]]){
    const {mesh,exact}=annularVortex(nx,ny),r=solveEuler(mesh,{initial:mesh.cells.map(exact),alpha:120,spatialOrder:2});
    assert.equal(r.converged,true);assert.ok(r.diagnostics.residual<1e-9);
    assert.ok(r.diagnostics.relativeMassImbalance<1e-8);assert.ok(r.diagnostics.wallLeakage<1e-12);
    errors.push({p:weightedError(r,exact,'p'),u:weightedError(r,exact,'u')});
  }
  for(let i=1;i<errors.length;i++)for(const key of ['p','u'])assert.ok(errors[i][key]<.3*errors[i-1][key]);
  assert.ok(errors.at(-1).p<.0014);assert.ok(errors.at(-1).u<.0008);
});
