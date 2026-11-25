// SPDX-License-Identifier: GPL-2.0-or-later
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { createStreamtubeChannel, solveStreamtubeChannel } from '../src/euler/tests/streamtube-channel.js';
import { createUpwindStreamtubeChannel } from '../src/euler/tests/streamtube-upwind-channel.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { channelGas, nozzleChannel } from '../tests/oracles/streamtube.js';
import { transonicNozzleReference } from '../tests/oracles/transonic-nozzle.js';

const binary=fileURLToPath(new URL('./check-channel',import.meta.url));
function native(p, mode, initial, direction, iterations=40) {
  const nx=p.x.length-1, nt=p.massFlows.length;
  const vector=(v,n)=>Array.isArray(v) ? v : Array(n).fill(v ?? 0);
  const input=[
    [nx,nt,mode,iterations,+Boolean(direction),1e-10],
    [p.gamma,p.referenceDensity ?? 1,p.referencePressure,p.pressureCorrectionFactor ?? .1,
      p.upwind?.mucon ?? 1,p.upwind?.mcrit ?? .99,p.hybrid?.epsilonP ?? 0],
    p.x,p.lower,p.upper,p.massFlows,vector(p.stagnationEnthalpy,nt),vector(p.stagnationDensity,nt),
    vector(p.inletSlopes,nt-1),vector(p.outletSlopes,nt-1),vector(p.outletPressure,nt),
    Array.from(initial), ...(direction ? [Array.from(direction)] : []),
  ].map(row=>row.join(' ')).join('\n')+'\n';
  const rows=execFileSync(binary,{input,encoding:'utf8'}).trim().split('\n').map(row=>row.split(' ').map(Number));
  const [status,count,backtracks,residual,maxMach,nonzeros,factorNonzeros,factors,refactors,relativeLinearResidual,ms,reuses]=rows.at(-2);
  return {status,count,backtracks,residual,maxMach,nonzeros,factorNonzeros,factors,refactors,relativeLinearResidual,ms,reuses,
    x:rows.at(-1), ...(direction ? {r:rows[0],jv:rows[1]} : {})};
}

function close(actual, expected, tol, label) {
  assert.equal(actual.length,expected.length);
  for(const [i,x] of actual.entries())
    assert.ok(Math.abs(x-expected[i])<=tol*Math.max(1,Math.abs(expected[i])),
      `${label} ${i}: ${x} vs ${expected[i]}`);
}

function verify(p,mode,system,initial,label) {
  const direction=initial.map((_,i)=>Math.sin(.4+.9*i));
  const a=native(p,mode,initial,direction);
  close(a.r,system.residual(initial),2e-12,label+' residual');
  close(a.jv,sparseProduct(system.jacobian(initial,{sparse:true}),direction),2e-10,label+' Jacobian product');
  const h=1e-6;
  const plus=system.residual(initial.map((v,i)=>v+h*direction[i]));
  const minus=system.residual(initial.map((v,i)=>v-h*direction[i]));
  close(a.jv,plus.map((v,i)=>(v-minus[i])/(2*h)),2e-7,label+' central difference');
  const t=performance.now();
  const b=solveStreamtubeChannel(system,{initial,tolerance:1e-10});
  const js=performance.now()-t;
  assert.equal(b.converged,true,b.reason);
  assert.equal(a.status,0,`${label}: residual ${a.residual}, updates ${a.count}, factors ${a.factors}`);
  assert.ok(a.residual<=1e-10 && a.relativeLinearResidual<1e-10);
  assert.ok(a.factors>=1 && (a.count<2 || a.refactors+a.reuses>=1));
  assert.ok(a.nonzeros<35*system.n);
  close(a.x,b.x,3e-8,label+' converged state');
  console.log(`${label}: ${system.n} unknowns, ${a.count} updates, ${a.refactors} refactors, ${a.reuses} factor reuses, max Mach ${a.maxMach.toFixed(3)}, C ${a.ms.toFixed(1)} ms / JS ${js.toFixed(1)} ms`);
  return a;
}

for (const mode of [0,1,2]) {
  const p={...nozzleChannel(12,3),...(mode===1 ? {streamwiseMode:'isentropic'} : {})};
  if (mode===2) Object.assign(p,{...channelGas(.55),outletPressure:channelGas(.55).referencePressure,
    upwind:{mucon:1,mcrit:.75,boundary:{kind:'unfiltered-first-two'}}});
  const system=mode<2 ? createStreamtubeChannel(p) : createUpwindStreamtubeChannel(p);
  const initial=system.initial.map((_,i)=>.002*Math.sin(i+.6));
  verify(p,mode,system,initial,['momentum','isentropic','reservoir'][mode]);
  const r=native(p,mode,initial,null,0);
  assert.equal(r.status,1);
  assert.deepEqual(r.x,Array.from(initial));
  const folded=initial.slice();
  folded[system.positionIndex(3,1)]=2;
  assert.equal(native(p,mode,folded,null).status,2);
}

for (const [nx,nt,shear] of [[32,3,0],[64,12,1]]) {
  const exact=transonicNozzleReference({shockMach:1.5,areaRatio:2,throatArea:.02});
  const x=Array.from({length:nx+1},(_,i)=>exact.domain[0]+(exact.domain[1]-exact.domain[0])*i/nx);
  const p={x,lower:x.map(z=>shear*z-.5*exact.area(z)),upper:x.map(z=>shear*z+.5*exact.area(z)),
    inletSlopes:Array(nt-1).fill(shear),massFlows:Array(nt).fill(exact.chokedMass/nt/Math.hypot(1,shear)),
    stagnationEnthalpy:exact.stagnationEnthalpy,stagnationDensity:exact.stagnationDensity,
    gamma:1.4,referenceDensity:1,referencePressure:exact.totalPressure,outletPressure:exact.backPressure,
    streamwiseMode:'hybrid',hybrid:{epsilonP:1e-5},
    upwind:{mucon:1,mcrit:.99,boundary:{kind:'unfiltered-first-two'}}};
  const system=createUpwindStreamtubeChannel(p), initial=system.initial;
  assert.ok(system.evaluate(initial).sections.flat().every(s=>s.machSquared<1));
  const result=verify(p,3,system,initial,`transonic nozzle ${nx}×${nt}`);
  assert.ok(result.maxMach>Math.sqrt(1.5));
  const v=system.evaluate(Float64Array.from(result.x));
  const mach=v.sections.map(row=>Math.sqrt(row[0].machSquared));
  let shock=1;
  for(let i=2;i<nx;i++) if(mach[i-1]-mach[i]>mach[shock-1]-mach[shock]) shock=i;
  assert.ok(Math.abs(x[shock]-exact.shockX)<1.5*(x[1]-x[0]));
  assert.ok(Math.abs(v.massFlows.reduce((a,b)=>a+b,0)*Math.hypot(1,shear)/exact.chokedMass-1)<6/nx**2);
}
