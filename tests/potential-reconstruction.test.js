import test from 'node:test';
import assert from 'node:assert/strict';
import { channelMesh } from '../src/euler/mesh.js';
import { affine,evaluateAffine,packAffine,evaluatePackedAffine,potentialReconstruction } from '../src/potential/reconstruction.js';

test('packed affine evaluation preserves cancellation order and is an immutable coefficient snapshot',()=>{
  const form=affine(.25,[[3,1e20],[0,-1e20],[2,.75],[1,1e-20]]),x=Float64Array.from([1,1,1,1]),packed=packAffine(form);
  assert.equal(evaluatePackedAffine(packed,x),evaluateAffine(form,x));
  const value=evaluatePackedAffine(packed,x);form.coefficients.set(3,0);
  assert.equal(evaluatePackedAffine(packed,x),value);
});

test('potential gradients reproduce a quadratic field on a skew grid, including boundary derivative data',()=>{
  const mesh=channelMesh({nx:6,ny:4,map:({x,y,i,j,nx,ny})=>({x:x+.06*Math.sin(Math.PI*i/nx)*Math.sin(Math.PI*j/ny),y:y+.05*Math.sin(2*Math.PI*i/nx)*Math.sin(Math.PI*j/ny)})});
  const exact=({x,y})=>.3+x-.2*y+.4*x*x+.3*x*y-.1*y*y;
  const derivative=({x,y})=>[1+.8*x+.3*y,-.2+.3*x-.2*y];
  const r=potentialReconstruction(mesh,{boundaryPotential:p=>affine(exact(p)),boundaryNormalVelocity:f=>{const [u,v]=derivative(f);return affine(u*f.nx+v*f.ny);}});
  const x=Float64Array.from(mesh.cells,exact);assert.equal(r.linearFallbacks,0);
  for(const face of mesh.faces)for(const cell of [face.owner,...(face.neighbor===null?[]:[face.neighbor])]){
    const expected=derivative(face),g=r.gradient(cell,face);
    for(let k=0;k<2;k++)assert.ok(Math.abs(evaluateAffine(g[k],x)-expected[k])<2e-12);
  }
});

test('potential reconstruction retains affine dependence on circulation and normal boundary velocity',()=>{
  const mesh=channelMesh({nx:5,ny:3}),n=mesh.cells.length;
  const baseVelocity=()=>[affine(.7,[[n,.3]]),affine(.1,[[n,-.2]])];
  const r=potentialReconstruction(mesh,{baseVelocity,boundaryNormalVelocity:f=>affine(0,[[n+1,f.ny]])});
  const x=Float64Array.from({length:n+2},(_,i)=>.03*Math.sin(.17*i)),d=x.map((_,i)=>Math.cos(.37*i)),h=1e-5;
  for(const [cell,c] of mesh.cells.entries())for(const form of r.velocity(cell,c)){
    const plus=evaluateAffine(form,x.map((v,i)=>v+h*d[i])),minus=evaluateAffine(form,x.map((v,i)=>v-h*d[i]));
    let derivative=0;for(const [col,value] of form.coefficients)derivative+=value*d[col];
    assert.ok(Math.abs((plus-minus)/(2*h)-derivative)<1e-10);
  }
});
