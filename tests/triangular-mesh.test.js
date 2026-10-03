import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4,transform,signedArea,pointInside } from '../src/geometry/airfoil.js';
import { triangularAirfoilMesh,deformTriangularWakes } from '../src/euler/triangular-mesh.js';
import { affine } from '../src/potential/reconstruction.js';
import { createPotentialSystem } from '../src/potential/system.js';

test('constrained triangular mesh preserves all body holes, curved wakes, and the exact fluid area',()=>{
  const contours=[naca4('0012',80),transform(naca4('0012',40),{chord:.3,x:1.05,y:-.2})];
  const wakePaths=contours.map(p=>[p[0],{x:p[0].x+.02,y:p[0].y+.001},{x:p[0].x+.2,y:p[0].y+.01},{x:p[0].x+1,y:p[0].y+.02}]);
  const m=triangularAirfoilMesh(contours,{wakePaths,surfaceScale:.8}),b=m.bounds;
  const expected=(b.xmax-b.xmin)*(b.ymax-b.ymin)-contours.reduce((s,c)=>s+Math.abs(signedArea(c)),0);
  assert.ok(Math.abs(m.cells.reduce((s,c)=>s+c.area,0)-expected)<1e-10);
  assert.ok(m.cells.every(c=>!contours.some(p=>pointInside(c,p))));assert.ok(m.diagnostics.metricClosure<1e-14);
  assert.equal(m.faces.filter(f=>f.boundary?.type==='wall').length,240);
  assert.equal(m.cuts.length,12);
  for(const cut of m.cuts)assert.notEqual(m.faces[cut.face].neighbor,null);
  for(const path of wakePaths)for(const p of path)assert.ok(m.vertices.some(v=>Math.hypot(v.x-p.x,v.y-p.y)<1e-12));
  const moved=wakePaths.map(path=>path.map((p,i)=>({...p,y:p.y+(i===0?0:1e-5)}))),next=deformTriangularWakes(m,moved);
  assert.ok(next);assert.deepEqual(next.cells.map(c=>c.vertices),m.cells.map(c=>c.vertices));
  assert.ok(Math.abs(next.cells.reduce((s,c)=>s+c.area,0)-expected)<1e-10);
  for(const cut of next.cuts){
    const a=moved[cut.element][cut.segment],b=moved[cut.element][cut.segment+1],f=next.faces[cut.face];
    for(const id of [f.a,f.b]){
      const p=next.vertices[id];
      assert.ok(Math.abs((p.x-a.x)*(b.y-a.y)-(p.y-a.y)*(b.x-a.x))<1e-13);
    }
  }
  // Uniform flow with its exact normal data exercises the full geometric
  // conservation law independently of the triangulation connectivity.
  const s=createPotentialSystem(m,{mach:.4,boundaryNormalVelocity:f=>affine(f.nx)});
  assert.ok(s.evaluate(new Float64Array(s.n)).diagnostics.residual<1e-12);
});

test('triangular mesh accepts a slat/main/flap layout without fictitious inlet dividing walls',()=>{
  const contours=[transform(naca4('0012',40),{chord:.2,x:-.22,y:.05,angle:-12}),naca4('2412',80),transform(naca4('0012',40),{chord:.3,x:.94,y:-.08,angle:-15})];
  const m=triangularAirfoilMesh(contours,{surfaceScale:1,padding:2});
  assert.equal(new Set(m.faces.filter(f=>f.boundary?.type==='wall').map(f=>f.boundary.element)).size,3);
  assert.equal(m.cuts.length,0);assert.ok(m.cells.length>1000);
  assert.throws(()=>triangularAirfoilMesh(contours,{maxVertices:100}),/limit/);
});
