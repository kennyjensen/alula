import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createCoupledAssembly, solveCoupledAssembly } from '../src/viscous/assembly.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { clcalc } from '../src/viscous/xfoil/xfoil.js';
import { normInf } from '../src/numerics/linear.js';
import { solveMultielementViscous } from '../src/viscous/multielement.js';

test('multielement BL Jacobian includes displacement, wake, stagnation and trip coordinate chains',()=>{
  const elements=[{points:naca4('0012',40)},{points:transform(naca4('0012',40),{y:-2})}];
  const s=createCoupledAssembly({elements,alpha:4,wakeCount:16,trips:[.1,.1]}),x=s.initial,n=x.length,j=s.jacobian(x);
  const d=x.map((v,i)=>Math.sin(.271*i)*Math.max(Math.abs(v),i%4===0?.03:.05)),h=1e-6;
  const rp=s.residual(x.map((v,i)=>v+h*d[i])),rm=s.residual(x.map((v,i)=>v-h*d[i]));
  for(let row=0;row<n;row++){
    let a=0;for(let col=0;col<n;col++)a+=j[row*n+col]*d[col];
    const b=(rp[row]-rm[row])/(2*h);
    assert.ok(Math.abs(a-b)<1e-6*Math.max(1,Math.abs(a),Math.abs(b)),`row ${row}: ${a} vs ${b}`);
  }
});

test('new simultaneous assembly agrees with the independently compiled Fortran forced-transition case',async()=>{
  const fixture=JSON.parse(await readFile(new URL('./fixtures/fortran/coupled.json',import.meta.url)));
  const c=fixture.cases.find(c=>c.name==='forced-transition');
  const r=solveCoupledAssembly({elements:[{points:c.points}],...c.options,wakeLength:1});
  assert.equal(r.converged,true);assert.ok(normInf(r.system.residual(r.x))<1e-8);
  const q=r.x.filter((_,i)=>i%4===3),forces=clcalc(c.points.length,c.points.map(p=>p.x),c.points.map(p=>p.y),q,null,c.options.alpha*Math.PI/180,0,1);
  const w=r.states.at(-1),cd=2*w.theta*w.ue**(.5*(5+w.deltaStar/w.theta));
  // Different source and wake discretizations: these are resolution-error
  // budgets, not bitwise parity tolerances for the underlying BL kernels.
  assert.ok(Math.abs(forces.cl-c.expected.cl)<.001);
  assert.ok(Math.abs(forces.cm-c.expected.cm)<.0002);
  assert.ok(Math.abs(cd-c.expected.cd)<1e-5);
});

test('two elements solve all four BL surfaces and both wakes in one system',()=>{
  const elements=[{points:naca4('0012',80)},{points:transform(naca4('0012',80),{y:-2})}];
  const r=solveCoupledAssembly({elements,alpha:4,trips:[.1,.1],wakeLength:1});
  assert.equal(r.converged,true);assert.equal(r.system.surfaces.length,4);assert.equal(r.system.outer.wakes.length,2);
  assert.ok(normInf(r.system.residual(r.x))<1e-8);
  const {outer}=r.system;
  for(const w of outer.wakes){
    const upper=r.states[w.body.start],lower=r.states[w.body.end],first=r.states[w.start];
    assert.ok(Math.abs(first.theta-upper.theta-lower.theta)<1e-10);
    assert.ok(Math.abs(first.ue*first.deltaStar-upper.ue*upper.deltaStar-lower.ue*lower.deltaStar)<1e-10);
  }
  // A frozen single-body velocity field cannot satisfy the shared equations.
  const frozen=r.x.slice();for(let i=0;i<outer.total;i++)frozen[4*i+3]=r.system.initial[4*i+3];
  assert.ok(normInf(r.system.residual(frozen))>1e-3);
});

test('close main/flap flow crosses a stagnation node and closes both coupled wake shapes',()=>{
  const elements=[{points:naca4('0012',160)},{points:transform(naca4('0012',80),{chord:.3,x:1.05,y:-.1})}];
  let crossed=false;
  const r=solveMultielementViscous({elements,alpha:4,trips:[.05,.1],wakeLength:1},{wakeRelaxation:.8,
    onIteration:h=>{crossed ||= h.activeChange===true;}});
  assert.equal(r.converged,true);assert.equal(crossed,true);assert.equal(r.wakeConverged,true);
  assert.ok(r.wakeResidual<1e-6);assert.ok(normInf(r.system.residual(r.x))<1e-8);
  const mass=Float64Array.from(r.states,(s,i)=>r.x[4*i+3]*s.deltaStar),field=r.system.outer.velocityField(mass);
  for(const wake of r.system.outer.wakes){
    let curvature=0;
    for(const p of wake.segments){const v=field(p);assert.ok(Math.abs(v.u*p.nx+v.v*p.ny)<1e-6);curvature=Math.max(curvature,Math.abs(p.ty));}
    assert.ok(curvature>.005,'wake must respond to the coupled flow');
  }
});

test('symmetric stagnation has zero flux without an artificial edge-velocity floor',()=>{
  const r=solveMultielementViscous({elements:[{points:naca4('0012',80)}],alpha:0,wakeLength:1});
  assert.equal(r.converged,true);assert.equal(r.system.zeroNodes[0],40);assert.ok(Math.abs(r.x[4*40+3])<1e-10);
  for(let i=0;i<40;i++)assert.ok(Math.abs(r.x[4*i+3]+r.x[4*(80-i)+3])<1e-6);
});

test('a converged BL with an unclosed wake cannot report coupled convergence',()=>{
  const r=solveMultielementViscous({elements:[{points:naca4('0012',80)}],alpha:4,trips:[.1,.1],wakeLength:1},{maxWakeIterations:0});
  assert.equal(r.converged,false);assert.equal(r.wakeConverged,false);assert.equal(r.reason,'wake iteration limit');
  assert.ok(r.history.at(-1).residual<1e-8);assert.ok(r.wakeResidual>1e-3);
});

test('an extra edge-flow reversal cannot be hidden by taking absolute velocities',()=>{
  const s=createCoupledAssembly({elements:[{points:naca4('0012',80)}],alpha:4});
  assert.equal(s.admissible(s.initial),true);
  const x=s.initial.slice(),{st}=s.decode(x);x[4*(st[0].right+3)+3]*=-1;
  assert.equal(s.admissible(x),false);
});

test('free transition and a near-uniform wake converge without relaxing the equation tolerance',()=>{
  const elements=[{points:naca4('2412',80)},{points:transform(naca4('0012',80),{chord:.3,x:1.05,y:-.1,angle:-5})}];
  const r=solveMultielementViscous({elements,alpha:4},{maxIterations:60,wakeRelaxation:.8});
  assert.equal(r.converged,true);assert.ok(normInf(r.system.residual(r.x))<1e-8);assert.ok(r.wakeResidual<1e-6);
});
