import test from 'node:test';
import assert from 'node:assert/strict';
import { sparseMatrix,sparseAdd,sparseDense,sparseProduct,solveSparse } from '../src/numerics/sparse.js';
import { solveLinear } from '../src/numerics/linear.js';
import { annularVortex } from './oracles/euler.js';
import { eulerJacobian } from '../src/euler/jacobian.js';
import { linearReconstruction } from '../src/euler/reconstruction.js';
import { temporalJacobian } from '../src/euler/pseudo-transient.js';
import { freestream,conserved,totalConditions } from '../src/euler/gas.js';
import { solveEuler,encodeStates,decodeStates } from '../src/euler/solve.js';

test('sparse GMRES solves a nonsymmetric two-dimensional stencil and checks its true residual',()=>{
  const width=20,n=width*width,columns=Array.from({length:n},(_,i)=>[i,...(i%width?[i-1]:[]),...(i%width<width-1?[i+1]:[]),...(i>=width?[i-width]:[]),...(i<n-width?[i+width]:[])]);
  const a=sparseMatrix(columns);
  for(let i=0;i<n;i++)for(const j of columns[i])sparseAdd(a,i,j,i===j?4.1:j===i-1?-1.2:j===i+1?-.8:-1);
  const exact=Float64Array.from({length:n},(_,i)=>Math.sin(.37*i)),rhs=sparseProduct(a,exact);
  for(const fillLevel of [0,1,2]){
    const r=solveSparse(a,rhs,{fillLevel});assert.ok(r.iterations>1);assert.ok(r.relativeResidual<1e-10);
    assert.ok(Math.max(...r.x.map((v,i)=>Math.abs(v-exact[i])))<2e-9);
  }
  let interrupted;
  try{solveSparse(a,rhs,{maxIterations:1});}catch(error){interrupted=error;}
  assert.equal(interrupted.code,'GMRES_ITERATION_LIMIT');
  const resumed=solveSparse(a,rhs,{initial:interrupted.solution});
  assert.ok(resumed.relativeResidual<1e-10);
  assert.ok(Math.max(...resumed.x.map((v,i)=>Math.abs(v-exact[i])))<2e-9);
  assert.deepEqual(solveSparse(a,new Float64Array(n),{initial:exact}).x,new Float64Array(n));
});

test('sparse Euler assembly is identical to dense assembly and gives the same regularized Newton step',()=>{
  const {mesh,exact}=annularVortex(8,3),states=mesh.cells.map(exact),reference=freestream({alpha:120}),reconstruction=linearReconstruction(mesh);
  const dense=eulerJacobian(mesh,states,reference,{reconstruction}),sparse=eulerJacobian(mesh,states,reference,{reconstruction,sparse:true});
  assert.deepEqual(sparseDense(sparse),dense);
  const n=dense.length**.5,blocks=temporalJacobian(states,reference,30);
  for(let cell=0;cell<states.length;cell++)for(let row=0;row<4;row++)for(let col=0;col<4;col++){
    const i=4*cell+row,j=4*cell+col,value=blocks[cell][4*row+col];dense[i*n+j]+=value;sparseAdd(sparse,i,j,value);
  }
  const rhs=Float64Array.from({length:n},(_,i)=>Math.sin(.23*i)),expected=solveLinear(dense,rhs);
  for(const blockSize of [1,4]){
    const r=solveSparse(sparse,rhs,{blockSize});assert.ok(r.relativeResidual<1e-10);
    for(let i=0;i<n;i++)assert.ok(Math.abs(r.x[i]-expected[i])<2e-8*Math.max(1,Math.abs(expected[i])));
  }
});

test('pseudo-time Jacobian differentiates conserved variables in the encoded Euler coordinates',()=>{
  const reference=freestream(),state={rho:.9,u:1.12,v:.07,p:reference.p-.2},x=encodeStates([state],reference),h=1e-6,cfl=17;
  const j=temporalJacobian([state],reference,cfl)[0],speed=Math.hypot(reference.u,reference.v),mass=reference.rho*speed;
  const scales=[mass,mass*speed,mass*speed,mass*totalConditions(reference).h0],wave=Math.hypot(state.u,state.v)+Math.sqrt(reference.gamma*state.p/state.rho);
  for(let col=0;col<4;col++){
    const plus=x.slice(),minus=x.slice();plus[col]+=h;minus[col]-=h;
    const a=conserved(decodeStates(plus,1,reference)[0]),b=conserved(decodeStates(minus,1,reference)[0]);
    for(let row=0;row<4;row++)assert.ok(Math.abs(j[4*row+col]-(a[row]-b[row])*wave/(2*h*cfl*scales[row]))<2e-10);
  }
});

test('sparse Euler solve converges to the same second-order vortex solution as full Newton',()=>{
  const {mesh,exact}=annularVortex(8,3),input={initial:mesh.cells.map(exact),alpha:120,spatialOrder:2};
  const dense=solveEuler(mesh,input),sparse=solveEuler(mesh,{...input,linearBackend:'sparse'});
  assert.equal(sparse.converged,true);assert.equal(dense.converged,true);
  for(let i=0;i<mesh.cells.length;i++)for(const key of ['rho','u','v','p'])assert.ok(Math.abs(sparse.states[i][key]-dense.states[i][key])<1e-9);
  const direct=solveEuler(mesh,{...input,linearBackend:'klu'});
  assert.equal(direct.converged,true);
  for(let i=0;i<mesh.cells.length;i++)for(const key of ['rho','u','v','p'])assert.ok(Math.abs(direct.states[i][key]-dense.states[i][key])<1e-9);
});
