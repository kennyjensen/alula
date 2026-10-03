import test from 'node:test';
import assert from 'node:assert/strict';
import { solveBlSchur } from '../src/numerics/bl-schur.js';
import { solveLinear, normInf } from '../src/numerics/linear.js';
import { createCoupledAssembly } from '../src/viscous/assembly.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

test('BL Schur step matches full simultaneous factorization with all element interactions',()=>{
  const elements=[{points:naca4('0012',40)},{points:transform(naca4('0012',40),{chord:.3,x:1.05,y:-.1})}];
  const s=createCoupledAssembly({elements,alpha:4,trips:[.05,.1],wakeCount:16}),x=s.initial;
  s.updateActive(x);
  const jac=s.jacobian(x),rhs=s.residual(x).map(v=>-v),full=solveLinear(jac,rhs),fast=s.linearSolve(jac,rhs,x);
  assert.equal(s.linearHistory.at(-1).method,'schur');
  assert.ok(s.linearHistory.at(-1).backwardError<1e-12);
  assert.ok(normInf(fast.map((v,i)=>v-full[i]))<1e-9*Math.max(1,normInf(full)));
});

test('singular local BL block falls back to the invertible coupled matrix',()=>{
  const a=Float64Array.from([0,0,0,1, 0,1,0,0, 0,0,1,0, 1,0,0,1]);
  const r=solveBlSchur(a,[4,2,3,5],[0]);
  assert.equal(r.method,'full');assert.deepEqual([...r.x],[1,2,3,4]);assert.equal(r.backwardError,0);
});
