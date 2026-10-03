import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4,transform } from '../src/geometry/airfoil.js';
import { solveSubcriticalAssembly } from '../src/potential/result.js';
import { buildSubcriticalMeshPreview } from '../src/potential/mesh-preview.js';

const elements=[{name:'Main',points:naca4('0012',40)},{name:'Flap',trips:[.08,.12],points:transform(naca4('0012',40),{chord:.3,x:1.05,y:-.2})}];
const input={elements,alpha:2,mach:.2,trips:[.05,.05],wakeCount:24,wakeLength:1};
test('subcritical result exports four BLs, two wakes, consistent reference scaling and the finite-density mass budget',()=>{
  const controls={mesh:{padding:4}},a=solveSubcriticalAssembly(input,controls);
  const b=solveSubcriticalAssembly({...input,referenceChord:2,elements:elements.map(e=>({...e,points:e.points.map(p=>({x:2*p.x,y:2*p.y}))}))},controls);
  for(const r of [a,b]){
    assert.equal(r.status,'solved');assert.equal(r.model,'multielement-coupled-subcritical');
    assert.equal(r.solverSettings.edgeVelocitySampling,'station-average');
    assert.equal(r.solverSettings.bodyMassInterpolation,'hermite');
    assert.equal(r.boundaryLayer.surfaces.length,4);assert.equal(r.boundaryLayer.wakes.length,2);
    assert.equal(r.mesh.cells.length,r.diagnostics.cells);
    assert.ok(r.mesh.cells.every(cell=>cell.length===3&&cell.every(i=>Number.isInteger(i)&&i>=0&&i<r.mesh.vertices.length)));
    assert.ok(r.diagnostics.equationResidual<1e-8);assert.ok(r.diagnostics.wakeResidual<1e-6);assert.ok(r.diagnostics.massBudgetError<1e-10);
    assert.ok(r.diagnostics.maxMach>.2&&r.diagnostics.maxMach<1);
    assert.deepEqual(r.solverSettings.elementTrips,[[.05,.05],[.08,.12]]);
    r.diagnostics.transition.forEach((v,i)=>assert.ok(Math.abs(v-[.05,.05,.08,.12][i])<2e-4,'independent material trips on every surface'));
    assert.ok(Math.abs(r.cd-r.boundaryLayer.wakes.reduce((s,w)=>s+w.drag,0))<1e-14);
    for(const w of r.boundaryLayer.wakes)for(const value of Object.values(w.matching))assert.ok(Math.abs(value)<1e-11);
    assert.doesNotThrow(()=>JSON.stringify(structuredClone(r)));
  }
  for(const key of ['cl','cm','cd','cdf'])assert.ok(Math.abs(a[key]-b[key])<1e-12,key);
  assert.equal(b.boundaryLayer.wakes[0].stations.at(-1).deltaStar,2*a.boundaryLayer.wakes[0].stations.at(-1).deltaStar);
  assert.deepEqual(b.mesh.cells,a.mesh.cells);
  a.mesh.vertices.forEach((p,i)=>{assert.equal(b.mesh.vertices[i].x,2*p.x);assert.equal(b.mesh.vertices[i].y,2*p.y);});
});

test('subcritical API rejects invalid conditions and labels finite last-iterate forces as unconverged',()=>{
  for(const change of [{mach:1},{mach:NaN},{wakeCount:200},{wakeLengths:[1]},{trips:[0,1]}])assert.throws(()=>solveSubcriticalAssembly({...input,...change}));
  const preview=buildSubcriticalMeshPreview(input,{mesh:{padding:4}}),events=[],meshes=[];
  const r=solveSubcriticalAssembly({...input,maxIterations:1},{mesh:{padding:4},
    onMesh:mesh=>{events.push('mesh');meshes.push(mesh);},onIteration:()=>events.push('iteration')});
  assert.equal(events[0],'mesh');assert.ok(events.includes('iteration'));
  assert.deepEqual(meshes[0],preview,'preview is the actual initial computational grid');
  assert.deepEqual(meshes.at(-1),r.mesh,'latest mesh event matches the returned computational grid');
  assert.equal(r.status,'unconverged'); assert.equal(r.coefficientStatus,'unconverged');
  for(const key of ['cl','cd','cm','cdf'])assert.ok(Number.isFinite(r[key]),key);
  assert.ok(r.warnings.some(w=>w.includes('provisional'))); assert.ok(r.diagnostics.equationResidual>1e-8);
});
