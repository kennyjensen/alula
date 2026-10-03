import test from 'node:test';
import assert from 'node:assert/strict';
import { channelMesh } from '../src/euler/mesh.js';
import { affine } from '../src/potential/reconstruction.js';
import { partitionPotentialInterfaces,potentialFluxConstraint } from '../src/potential/interfaces.js';
import { createPotentialSystem,solvePotential } from '../src/potential/system.js';

function problem(mach){
  const original=channelMesh({nx:6,ny:4});
  const ids=original.faces.map((f,i)=>f.neighbor!==null&&Math.abs(f.y-.5)<1e-12&&Math.abs(f.ny)>0.9?i:-1).filter(i=>i>=0);
  const mesh=partitionPotentialInterfaces(original,ids),count=mesh.cells.length;
  const phi=p=>.1*Math.abs(p.y-.5),vn=p=>.1+(p.y<.5?-.1:.1),u=.8;
  const rho=v=>(1+.2*mach*mach*(1-u*u-v*v))**2.5;
  const source=rho(.2)*.2-rho(0)*0;
  const system=createPotentialSystem(mesh,{mach,baseVelocity:()=>[affine(u),affine(.1)],
    boundaryPotential:f=>f.boundary.type==='potential-interface'?affine(0,[[count+f.boundary.interface,1]]):affine(phi(f)),
    boundaryNormalVelocity:f=>affine(u*f.nx+vn(f)*f.ny),
    constraints:({fluxColumns})=>mesh.potentialInterfaces.map(p=>{
      const length=mesh.faces[p.faces[0]].length;
      return potentialFluxConstraint(mesh,fluxColumns,p.faces,{source:affine(source*length),scale:length});
    })});
  return{mesh,system,phi,source,exact:Float64Array.from([...mesh.cells,...mesh.potentialInterfaces.map(p=>mesh.faces[p.faces[0]])],phi)};
}

test('a wake has continuous potential, two normal velocities, and an exactly conserved compressible mass jump',()=>{
  for(const mach of [0,.4]){
    const {mesh,system,exact,source}=problem(mach),a=system.evaluate(exact);
    assert.ok(a.diagnostics.residual<1e-12);assert.ok(Math.abs(a.diagnostics.interfaceMass+source*2)<1e-12);
    assert.ok(Math.abs(a.diagnostics.externalMass-source*2)<1e-12);
    const r=solvePotential(system,{tolerance:1e-12});assert.equal(r.converged,true,r.reason);
    assert.ok(Math.max(...r.x.map((v,i)=>Math.abs(v-exact[i])))<1e-10);
    for(const p of mesh.potentialInterfaces){
      const [a,b]=p.faces.map(i=>r.faces[i]);assert.ok(Math.abs(a.u-.8)<1e-11);assert.ok(Math.abs(b.u-.8)<1e-11);
      assert.ok(Math.abs(Math.abs(a.v-b.v)-.2)<1e-11);
    }
  }
});

test('wake interface Jacobian includes both densities and the shared face potential',()=>{
  const {system,exact}=problem(.4),x=exact.map((v,i)=>v+.001*Math.sin(.63*i)),d=x.map((_,i)=>Math.cos(.31*i)),h=1e-6;
  const a=system.evaluate(x,{jacobian:true}),plus=system.evaluate(x.map((v,i)=>v+h*d[i])).residual,minus=system.evaluate(x.map((v,i)=>v-h*d[i])).residual;
  for(let row=0;row<system.n;row++){
    let jv=0;for(let col=0;col<system.n;col++)jv+=a.jacobian[row*system.n+col]*d[col];
    assert.ok(Math.abs(jv-(plus[row]-minus[row])/(2*h))<1e-7,`row ${row}`);
  }
});
