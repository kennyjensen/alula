// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import { solveViscousAssembly } from '../src/viscous/result.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

const fixture=JSON.parse(await readFile('tests/fixtures/fortran/coupled.json','utf8'));
const report={generatedAt:new Date().toISOString(),runtime:process.version,status:'passed',
  scope:'Simultaneous multielement incompressible BL/displacement/wake coupling. Native single-element comparisons, surface/wake refinement, and asymptotic decoupling; not experimental validation or compressible Euler validation.',
  nativeProvenance:fixture.provenance,native:[],surfaceRefinement:[],wakeRefinement:[],separationLimit:[],sha256:{}};
for(const [path,hash] of Object.entries(fixture.provenance.sha256))assert.equal(createHash('sha256').update(await readFile(path)).digest('hex'),hash,path);
const summarize=r=>({cl:r.cl,cm:r.cm,cd:r.cd,diagnostics:r.diagnostics});
const solve=input=>{const r=solveViscousAssembly(input);assert.equal(r.status,'solved',JSON.stringify({alpha:input.alpha,trips:input.trips,elements:input.elements.map(e=>({panels:e.points.length-1,te:e.points[0]})),diagnostics:r.diagnostics}));assert.ok(r.diagnostics.equationResidual<1e-8);assert.ok(r.diagnostics.wakeResidual<1e-6);return r;};
for(const f of fixture.cases){
  if(!f.expected.converged){
    const r=solveViscousAssembly({elements:[{points:f.points}],...f.options});
    assert.equal(r.status,'unconverged');assert.equal(r.cl,null);assert.equal(r.cd,null);
    report.native.push({name:f.name,conditions:f.options,converged:false,coefficientsWithheld:true});
    continue;
  }
  const r=solve({elements:[{points:f.points}],...f.options,wakeLength:1});
  const errors={cl:Math.abs(r.cl-f.expected.cl),cm:Math.abs(r.cm-f.expected.cm),cd:Math.abs(r.cd-f.expected.cd),
    cpMax:Math.max(...f.expected.cp.filter(p=>p.index<=f.points.length).map(p=>Math.abs(r.elements[0].cp[p.index-1].cp-p.cp))),
    transitionMax:Math.max(...f.expected.transition.map((x,i)=>Math.abs(x-r.diagnostics.transition[i])))};
  assert.ok(errors.cl<.001&&errors.cm<.0003&&errors.cd<1e-5,`${f.name}: ${JSON.stringify(errors)}`);
  report.native.push({name:f.name,conditions:f.options,panels:f.points.length-1,native:{cl:f.expected.cl,cm:f.expected.cm,cd:f.expected.cd},...summarize(r),errors});
  console.log(`Native comparison passed: ${f.name}`);
}
const closeCase=panels=>({elements:[{name:'Main',points:naca4('0012',panels)},
  {name:'Flap',points:transform(naca4('0012',panels/2),{chord:.3,x:1.05,y:-.1})}],alpha:4,trips:[.05,.1]});
for(const panels of [80,160,240]){
  const r=solve({...closeCase(panels),wakeLength:1});report.surfaceRefinement.push({panels:r.panelCount,...summarize(r)});
  console.log(`Surface refinement passed: ${r.panelCount} panels`);
}
for(const key of ['cl','cm','cd']){
  const [a,b,c]=report.surfaceRefinement;assert.ok(Math.abs(c[key]-b[key])<Math.abs(b[key]-a[key]),`${key} successive refinement changes must decrease`);
}
for(const [wakeLength,wakeCount] of [[1,24],[2,32],[4,40]]){
  const r=solve({...closeCase(160),wakeLength,wakeCount});report.wakeRefinement.push({wakeLength,wakeCount,...summarize(r)});
}
for(const key of ['cl','cd']){
  const [a,b,c]=report.wakeRefinement;assert.ok(Math.abs(c[key]-b[key])<Math.abs(b[key]-a[key]));
}
const base={elements:[{points:naca4('0012',80)}],alpha:4,trips:[.1,.1],wakeLength:1};
const single=solve(base);report.isolatedReference=summarize(single);
for(const distance of [4,8,16]){
  const r=solve({...base,elements:[base.elements[0],{points:transform(base.elements[0].points,{y:distance})}]});
  report.separationLimit.push({distance,...summarize(r),liftError:Math.abs(r.cl-2*single.cl),dragError:Math.abs(r.cd-2*single.cd)});
}
for(let i=1;i<report.separationLimit.length;i++)for(const key of ['liftError','dragError'])assert.ok(report.separationLimit[i][key]<report.separationLimit[i-1][key]);
const reversed=solve({...closeCase(80),elements:closeCase(80).elements.reverse(),wakeLength:1});
report.elementOrder={...summarize(reversed),liftDifference:Math.abs(reversed.cl-report.surfaceRefinement[0].cl),dragDifference:Math.abs(reversed.cd-report.surfaceRefinement[0].cd)};
assert.ok(report.elementOrder.liftDifference<1e-7&&report.elementOrder.dragDifference<1e-9);
for(const path of ['src/viscous/assembly.js','src/viscous/active-newton.js','src/viscous/integral.js','src/viscous/trips.js','src/viscous/multielement.js','src/viscous/result.js','src/inviscid/displacement.js','src/numerics/bl-schur.js','tests/fixtures/fortran/coupled.json','tests/fixtures/fortran/integral.json'])report.sha256[path]=createHash('sha256').update(await readFile(path)).digest('hex');
await mkdir('docs', { recursive: true });
await writeFile('docs/multielement-validation-results.json',JSON.stringify(report,null,2)+'\n');
console.log('Multielement validation passed. See docs/multielement-validation-results.json.');
