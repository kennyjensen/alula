import test from 'node:test';
import assert from 'node:assert/strict';
import { channelMesh } from '../src/euler/mesh.js';
import { makePanel } from '../src/inviscid/panel.js';
import { polynomialPotential,polynomialSourceBasis } from '../src/inviscid/displacement.js';
import { affine } from '../src/potential/reconstruction.js';
import { partitionPotentialInterfaces,potentialFluxConstraint } from '../src/potential/interfaces.js';
import { createPotentialSystem,solvePotential } from '../src/potential/system.js';

test('a finite source wake ends inside the fluid and converges to its independent analytic field without a downstream cut',()=>{
  const panel=makePanel({x:.5,y:.5},{x:1.5,y:.5}),coeff=[0,.6,-.6];
  const exactPhi=p=>polynomialPotential(p,panel).reduce((s,v,i)=>s+coeff[i]*v.real,0);
  const exactVelocity=p=>polynomialSourceBasis(p,panel).reduce((s,v,i)=>({u:s.u+coeff[i]*v.u,v:s.v+coeff[i]*v.v}),{u:1,v:0});
  const errors=[];
  for(const [nx,ny] of [[8,4],[16,8],[32,16]]){
    const original=channelMesh({nx,ny});for(const f of original.faces)if(f.boundary)f.boundary.type='farfield';
    const ids=original.faces.map((f,i)=>f.neighbor!==null&&Math.abs(f.y-.5)<1e-12&&f.x>.5&&f.x<1.5&&Math.abs(f.ny)>.9?i:-1).filter(i=>i>=0);
    const mesh=partitionPotentialInterfaces(original,ids),count=mesh.cells.length;
    const system=createPotentialSystem(mesh,{sparse:true,boundaryPotential:f=>f.boundary.type==='potential-interface'?affine(0,[[count+f.boundary.interface,1]]):affine(exactPhi(f)),
      constraints:({fluxColumns})=>mesh.potentialInterfaces.map(p=>{
        const f=mesh.faces[p.faces[0]],a=mesh.vertices[f.a].x-.5,b=mesh.vertices[f.b].x-.5,lo=Math.min(a,b),hi=Math.max(a,b);
        const source=.6*((hi*hi-lo*lo)/2-(hi**3-lo**3)/3);
        return potentialFluxConstraint(mesh,fluxColumns,p.faces,{source:affine(source),scale:f.length});
      })});
    const initial=Float64Array.from([...mesh.cells,...mesh.potentialInterfaces.map(p=>mesh.faces[p.faces[0]])],exactPhi);
    const r=solvePotential(system,{initial});assert.equal(r.converged,true,r.reason);
    assert.ok(Math.abs(r.diagnostics.externalMass-.1)<1e-9);assert.ok(Math.abs(r.diagnostics.interfaceMass+.1)<1e-9);
    let error=0,area=0;r.states.forEach((v,i)=>{const c=mesh.cells[i],e=exactVelocity(c);area+=c.area;error+=c.area*Math.hypot(v.u-e.u,v.v-e.v);});errors.push(error/area);
  }
  for(let i=1;i<errors.length;i++)assert.ok(errors[i]<.65*errors[i-1],String(errors));
  assert.ok(errors.at(-1)<.001,String(errors));
});
