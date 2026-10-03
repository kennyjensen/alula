import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { solveViscousAssembly } from '../src/viscous/result.js';
import { readFile } from 'node:fs/promises';
import { parseCoordinates } from '../src/geometry/parse.js';

test('multielement result preserves reference scaling, per-element forces and wake matching',()=>{
  const elements=[{name:'Main',points:naca4('0012',80)},{name:'Flap',points:transform(naca4('0012',40),{chord:.3,x:1.05,y:-.1})}];
  const input={elements,alpha:4,trips:[.05,.1],wakeLength:1,wakeLengths:[1,.3]};
  const a=solveViscousAssembly(input);
  const b=solveViscousAssembly({...input,referenceChord:2,wakeLengths:input.wakeLengths.map(v=>2*v),elements:elements.map(e=>({...e,points:e.points.map(p=>({x:2*p.x,y:2*p.y}))}))});
  assert.equal(a.status,'solved');assert.equal(b.status,'solved');
  assert.equal(a.boundaryLayer.surfaces.length,4);assert.equal(a.boundaryLayer.wakes.length,2);
  for(const key of ['cl','cm','cd','cdf'])assert.ok(Math.abs(a[key]-b[key])<1e-10,key);
  assert.ok(Math.abs(a.cl-a.elements.reduce((n,e)=>n+e.cl,0))<1e-12);
  assert.ok(Math.abs(a.cd-a.boundaryLayer.wakes.reduce((n,w)=>n+w.drag,0))<1e-12);
  for(const w of a.boundaryLayer.wakes)for(const value of Object.values(w.matching))assert.ok(Math.abs(value)<1e-10);
  for(const [e,s] of a.boundaryLayer.surfaces.entries())assert.ok(Math.abs(2*s.stations[0].theta-b.boundaryLayer.surfaces[e].stations[0].theta)<1e-10);
  assert.doesNotThrow(()=>structuredClone(a));assert.doesNotThrow(()=>JSON.stringify(a));
});

test('three interacting elements produce six BL surfaces and three separately closed wakes',()=>{
  const elements=[0,-2,2].map((y,i)=>({name:`Element ${i+1}`,points:transform(naca4('0012',80),{y})}));
  const r=solveViscousAssembly({elements,alpha:4,trips:[.1,.1],wakeLength:1});
  assert.equal(r.status,'solved');assert.equal(r.boundaryLayer.surfaces.length,6);assert.equal(r.boundaryLayer.wakes.length,3);
  assert.ok(r.diagnostics.equationResidual<1e-8);assert.ok(r.diagnostics.wakeResidual<1e-6);
  assert.ok(new Set(r.elements.map(e=>e.cl.toFixed(5))).size>1,'all elements must respond to their different neighbors');
});

test('public multielement API refuses finite Mach and labels finite last-iterate forces as unconverged',()=>{
  const elements=[{points:naca4('0012',80)}];
  assert.throws(()=>solveViscousAssembly({elements,mach:.2}),/Mach 0/);
  const r=solveViscousAssembly({elements,alpha:4,maxIterations:1});
  assert.equal(r.status,'unconverged'); assert.equal(r.coefficientStatus,'unconverged');
  for(const key of ['cl','cm','cd','cdf'])assert.ok(Number.isFinite(r[key]),key);
  assert.ok(r.warnings.some(w=>w.includes('provisional'))); assert.ok(r.diagnostics.equationResidual>1e-8);
});

test('original NASA NLR main/flap geometry closes four BL surfaces and both wakes',async()=>{
  const input=parseCoordinates(await readFile(new URL('../public/examples/blade.nlr7301-wind',import.meta.url),'utf8'));
  assert.deepEqual(input.elements.map(e=>e.points.length-1),[148,138]);
  const r=solveViscousAssembly({elements:input.elements,alpha:0,reynolds:1e6,trips:[.05,.05]});
  assert.equal(r.status,'solved');assert.equal(r.boundaryLayer.surfaces.length,4);assert.equal(r.boundaryLayer.wakes.length,2);
  assert.ok(r.diagnostics.equationResidual<1e-8);assert.ok(r.diagnostics.wakeResidual<1e-6);
  // The archive has conflicting flow-condition labels. This verifies solving
  // the actual geometry; its archived CFD forces are not an accuracy oracle.
  assert.ok(r.diagnostics.maxRelativeDisplacement<.1);
  for(const transition of r.diagnostics.transition)assert.ok(Math.abs(transition-.05)<1e-4,'all requested trips must be recovered after the assembly stagnation points move');
});
