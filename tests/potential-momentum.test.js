import test from 'node:test';
import assert from 'node:assert/strict';
import {channelMesh} from '../src/euler/mesh.js';
import {affine} from '../src/potential/reconstruction.js';
import {createPotentialSystem,solvePotential} from '../src/potential/system.js';
import {potentialMomentumAudit} from '../scripts/validation/potential-momentum.js';
import {annularVortex} from './oracles/euler.js';

test('independent momentum balance closes for an exact harmonic polynomial, with nonzero wall transpiration',()=>{
  const mesh=channelMesh({nx:5,ny:4}),phi=({x,y})=>.1*(x*x-y*y);
  const velocity=({x,y})=>[1+.2*x,-.2*y];
  const s=createPotentialSystem(mesh,{boundaryPotential:p=>affine(phi(p)),boundaryNormalVelocity:f=>{const [u,v]=velocity(f);return affine(u*f.nx+v*f.ny);}});
  const audit=potentialMomentumAudit(s,Float64Array.from(mesh.cells,phi));
  assert.ok(audit.momentumDefect<1e-13);assert.ok(Math.abs(audit.total.mass)<1e-13);
  assert.ok(audit.facePressureDifference<1e-14);
});

test('compressible vortex momentum error decreases independently of the converged mass equation',()=>{
  const errors=[];
  for(const [nx,ny]of [[6,2],[12,4],[24,8],[48,16],[96,32]]){
    const {mesh}=annularVortex(nx,ny),phi=({x,y})=>1.5*Math.atan2(y,x);
    const s=createPotentialSystem(mesh,{mach:.3,sparse:true,baseVelocity:()=>[affine(),affine()],boundaryPotential:p=>affine(phi(p))});
    const r=solvePotential(s,{initial:mesh.cells.map(phi)});assert.equal(r.converged,true);
    const audit=potentialMomentumAudit(s,r.x);assert.ok(audit.facePressureDifference<1e-13);errors.push(audit.momentumDefect);
  }
  for(let i=1;i<errors.length;i++)assert.ok(errors[i]<.4*errors[i-1],String(errors));
  assert.ok(errors.at(-1)<.003,String(errors));
});
