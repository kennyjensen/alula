import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4,transform } from '../src/geometry/airfoil.js';
import { multielementMesh } from '../src/euler/multielement-mesh.js';

test('curved wakes remain exact connected interfaces on a multielement mesh',()=>{
  const contours=[naca4('0012',40),transform(naca4('0012',40),{chord:.3,x:1.05,y:-.2})];
  const wakePaths=contours.map((p,e)=>[p[0],{x:p[0].x+.1,y:p[0].y+.03},{x:p[0].x+.6,y:p[0].y+.05},{x:p[0].x+1,y:p[0].y+.05+.01*e}]);
  const mesh=multielementMesh(contours,{rows:3,padding:2,wakePaths,outerGrowth:1.7});
  assert.ok(mesh.diagnostics.metricClosure<1e-14);
  for(const [e,path] of wakePaths.entries()){
    for(const p of path)assert.ok(mesh.vertices.some(v=>Math.hypot(v.x-p.x,v.y-p.y)<1e-12));
    for(const cut of mesh.cuts.filter(c=>c.element===e&&c.type==='wake-cut')){
      const f=mesh.faces[cut.face];assert.notEqual(f.neighbor,null);
      let j=1;while(j<path.length-1&&path[j].x<f.x)j++;
      const a=path[j-1],b=path[j],y=a.y+(b.y-a.y)*(f.x-a.x)/(b.x-a.x);
      assert.ok(Math.abs(f.y-y)<1e-12);
    }
  }
  const crossing=wakePaths.map(w=>w.map(p=>({...p})));
  crossing[1][1].y=.5;
  assert.throws(()=>multielementMesh(contours,{wakePaths:crossing,padding:2}),/cross/);
});

test('geometric normal spacing limits the jump between the first and second off-wall cells',()=>{
  const mesh=multielementMesh([naca4('0012',40)],{rows:24,padding:4,normalGrowth:1.5});
  const x=.5,column=mesh.vertices.filter(v=>Math.abs(v.x-x)<1e-12&&v.y>=0).map(v=>v.y).sort((a,b)=>a-b);
  const spacings=column.slice(1).map((y,i)=>y-column[i]);
  assert.equal(spacings.length,24);
  for(let i=1;i<spacings.length;i++)assert.ok(Math.abs(spacings[i]/spacings[i-1]-1.5)<1e-9);
  assert.ok(mesh.diagnostics.metricClosure<1e-14);
});
