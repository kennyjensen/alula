import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {buildMesh} from '../src/euler/mesh.js';
import {potentialReconstruction,affine,evaluateAffine} from '../src/potential/reconstruction.js';

test('an interior triangular stencil retains quadratic accuracy when its first two neighbor rings are rank deficient',()=>{
  const fixture=JSON.parse(readFileSync(new URL('./fixtures/rank-deficient-potential-mesh.json',import.meta.url)));
  const mesh=buildMesh(fixture.vertices,fixture.connectivity,fixture.boundaries);
  const potential=p=>p.x+.3*p.y+.2*(p.x*p.x-p.y*p.y)+.1*p.x*p.y;
  const exact=p=>[1+.4*p.x+.1*p.y,.3-.4*p.y+.1*p.x];
  const reconstruction=potentialReconstruction(mesh,{boundaryPotential:p=>affine(potential(p))});
  const state=Float64Array.from(mesh.cells,potential),cell=fixture.targetCell;
  for(const p of [mesh.cells[cell],...mesh.cells[cell].vertices.map(i=>mesh.vertices[i])]){
    const actual=reconstruction.gradient(cell,p).map(form=>evaluateAffine(form,state));
    actual.forEach((v,k)=>assert.ok(Math.abs(v-exact(p)[k])<2e-11,`quadratic gradient component ${k}`));
  }
});
