import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4,transform } from '../src/geometry/airfoil.js';
import { createCoupledPotential,solveCoupledPotential,solveMultielementPotentialViscous } from '../src/potential/coupled.js';
import { displacementMass } from '../src/viscous/mass.js';
import { sparseDense } from '../src/numerics/sparse.js';

const input={elements:[{points:naca4('0012',40)}],alpha:2,mach:.2,trips:[.05,.05],wakeCount:12,wakeLength:1};
test('the frozen BL incompressible guess shares the potential graph without changing the physical problem',()=>{
  const s=createCoupledPotential(input,{mesh:{rows:2,padding:2},farfield:'multipole',edgeVelocitySampling:'station-average'});
  assert.throws(()=>s.incompressibleInitialization(),/only an outer-flow initialization/);
  s.setFrozen(true);
  const before=s.evaluate(s.initial,{jacobian:true}),guess=s.incompressibleInitialization();
  for(const key of ['mesh','reconstruction','faces','cellVelocities','constraints'])assert.equal(guess[key],s[key],key);
  assert.equal(guess.conditions.mach,0);assert.equal(s.conditions.mach,.2);
  const initialized=guess.evaluate(s.initial,{jacobian:true});
  assert.ok(initialized.residual.every(Number.isFinite));
  assert.deepEqual(s.evaluate(s.initial,{jacobian:true}),before);
  s.setFrozen(false);
  assert.throws(()=>s.incompressibleInitialization(),/only an outer-flow initialization/);
});
test('complete potential/BL Jacobian includes reciprocal displacement, density, Kutta, and wake coupling',()=>{
  const s=createCoupledPotential(input,{mesh:{rows:2,padding:2},sparse:true}),x=s.initial.slice();s.updateActive(x);
  const d=x.map((v,i)=>Math.sin(.31*i+.2)*(i<s.blOffset?.0001:(i-s.blOffset)%4===3?.01:.001)),h=1e-5;
  const j=sparseDense(s.evaluate(x,{jacobian:true}).jacobian),rp=s.evaluate(x.map((v,i)=>v+h*d[i])).residual,rm=s.evaluate(x.map((v,i)=>v-h*d[i])).residual;
  for(let i=0;i<s.n;i++){
    let y=0;for(let k=0;k<s.n;k++)y+=j[i*s.n+k]*d[k];
    const fd=(rp[i]-rm[i])/(2*h);assert.ok(Math.abs(fd-y)/Math.max(1,Math.abs(fd))<1e-6,`row ${i}`);
  }
});

for(const edgeVelocitySampling of ['panel-average','station-average'])for(const bodyMassInterpolation of ['linear','hermite'])test(`the coupled Jacobian includes fitted farfield, ${edgeVelocitySampling} velocity and ${bodyMassInterpolation} displacement mass`,()=>{
  const s=createCoupledPotential(input,{mesh:{type:'triangular',surfaceScale:1.2,padding:2,growth:.5},
    farfield:'multipole',edgeVelocitySampling,bodyMassInterpolation}),x=s.initial.slice();
  x[s.count]=-.04;x[s.farfield.columns[0]]=.004;x[s.farfield.columns[1]]=.03;
  x[s.farfield.columns[2]]=-.01;x[s.farfield.columns[3]]=.04**2;s.updateActive(x);
  const d=x.map((v,i)=>Math.sin(.31*i+.2)*(i<s.blOffset?.0001:(i-s.blOffset)%4===3?.01:.001)),h=1e-5;
  const j=s.evaluate(x,{jacobian:true}).jacobian,rp=s.evaluate(x.map((v,i)=>v+h*d[i])).residual,rm=s.evaluate(x.map((v,i)=>v-h*d[i])).residual;
  for(let i=0;i<s.n;i++){
    let y=0;for(let k=j.rowPtr[i];k<j.rowPtr[i+1];k++)y+=j.values[k]*d[j.colIndex[k]];
    const fd=(rp[i]-rm[i])/(2*h);assert.ok(Math.abs(fd-y)/Math.max(1,Math.abs(fd))<1e-6,`row ${i}`);
  }
});

for(const edgeVelocitySampling of ['panel-average','station-average'])test(`finite-Mach two-element flow with ${edgeVelocitySampling} closes all four BLs with a conservative mass budget`,()=>{
  const elements=[...input.elements,{points:transform(naca4('0012',40),{chord:.3,x:1.05,y:-.2})}];
  const r=solveMultielementPotentialViscous({...input,elements},{mesh:{rows:3,padding:2},wakeRelaxation:.8,machContinuation:true,
    farfield:'multipole',edgeVelocitySampling,edgeVelocityContinuation:true});
  assert.equal(r.converged,true,r.reason);assert.equal(r.stage,'coupled');assert.ok(r.diagnostics.residual<1e-8);
  assert.equal(r.system.conditions.mach,input.mach);
  assert.ok(r.history.some(h=>h.continuationMach===0));assert.ok(r.history.some(h=>h.continuationMach===input.mach));
  assert.ok(r.history.some(h=>h.edgeVelocityFraction===0));assert.ok(r.history.some(h=>h.edgeVelocityFraction===1));
  assert.ok(r.system.farfield);const ff=r.system.farfield;
  const circulation=r.x.subarray(r.system.count,r.system.count+elements.length).reduce((s,v)=>s+v,0);
  assert.ok(Math.abs(r.x[ff.columns[3]]-circulation*circulation)<1e-10);
  assert.equal(r.wakeConverged,true);assert.ok(r.wakeResidual<1e-6);
  const {bl,outer,blOffset,mesh}=r.system,x=r.x.subarray(blOffset);
  assert.equal(bl.surfaces.length,4);assert.equal(outer.wakes.length,2);
  const mass=id=>displacementMass(x[4*id+3],x[4*id+2]*bl.thicknessScale,{mach:input.mach}).value;
  let farfield=0;mesh.faces.forEach((f,i)=>{if(f.boundary?.type==='farfield')farfield+=r.fluxes[i];});
  const downstream=outer.wakes.reduce((sum,w)=>sum+mass(w.end),0);
  assert.ok(Math.abs(farfield-downstream)<1e-10,'farfield flow must equal total surviving wake displacement mass');
  for(const w of outer.wakes)assert.ok(Math.abs(mass(w.start)-(mass(w.body.end)-mass(w.body.start)))<1e-11);
  // Changing the outer mesh must not destroy an already closed BL block
  // before the simultaneous Newton solve has taken its first step.
  const seed=bl.exportSeed(x),remeshed=solveCoupledPotential({...input,elements,seed,wakePaths:outer.wakes.map(w=>w.points)},
    {mesh:{rows:4,padding:2},maxIterations:0,farfield:'multipole',edgeVelocitySampling});
  const preserved=remeshed.x.subarray(remeshed.system.blOffset),block=remeshed.system.bl.residual(preserved);
  assert.ok(Math.max(...block.map(Math.abs))<1e-8);
  assert.ok(Math.max(...preserved.map((v,i)=>Math.abs(v-seed.x[i])))<1e-10);
});
