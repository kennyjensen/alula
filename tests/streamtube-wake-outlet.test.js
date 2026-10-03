// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { multipoleBasis } from '../src/potential/farfield.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { directChannelConservation } from './oracles/streamtube.js';
const max = a => Math.max(...Array.from(a, Math.abs));
// Independent physical condition: cross(unit segment, endpoint-average
// multipole velocity). No call to the residual/Jacobian outlet helper.
function bankAngles(system, flow) {
  const { center, alpha, mach, gamma } = system.euler.conditions;
  const c = flow.strengths, coeff = [-c.circulation,c.source,c.doubletX,c.doubletY,c.circulation**2];
  const velocity = p => multipoleBasis(p,{center,alpha,mach,gamma}).velocity.reduce((v,b,k)=>
    [v[0]+coeff[k]*b[0],v[1]+coeff[k]*b[1]], [Math.cos(alpha*Math.PI/180),Math.sin(alpha*Math.PI/180)]);
  return system.bl.wakes.flatMap(w => ['lower','upper'].map(side => {
    const g=side==='lower'?w.body:w.body+1, j=side==='lower'?system.euler.layout.tubes[g]:0;
    const a=flow.nodes[g].at(-2)[j],b=flow.nodes[g].at(-1)[j],va=velocity(a),vb=velocity(b);
    return ((b.x-a.x)*(va[1]+vb[1])-(b.y-a.y)*(va[0]+vb[0]))/(2*Math.hypot(b.x-a.x,b.y-a.y));
  }));
}
for (const wakeGeometry of ['centerline','independent-banks']) test(`${wakeGeometry}: both outlet banks close in a simultaneous two-element Euler/BL solve`, t => {
  const input={...intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2}),streamwiseMode:'isentropic',wakeGeometry,wakeOutlet:'banks'};
  const options={reynolds:1e5,edgeMatching:'section-velocity-distance'};
  const s=createCoupledStreamtubeBody(input,options),x=s.initial;
  const terminal=s.bl.wakes.map(w=>s.ne+4*w.ids.at(-1)+3);
  const matrix=s.jacobian(x),d=x.map((v,i)=>i<s.ne?.001*Math.sin(i*.7):Math.max(.01,Math.abs(v))*.01*Math.cos(i*.3));
  const jv=sparseProduct(matrix,d);
  for(const h of [2e-6,6e-7]) {
    const p=s.residual(x.map((v,i)=>v+h*d[i])),m=s.residual(x.map((v,i)=>v-h*d[i]));
    let error=0; jv.forEach((v,i)=>{const fd=(p[i]-m[i])/(2*h);error=Math.max(error,Math.abs(fd-v)/Math.max(1,Math.abs(v),Math.abs(fd)));});
    assert.ok(error<5e-6,JSON.stringify({h,error}));
  }
  // Every state column in both new terminal rows, including all independent
  // wake positions and multipole coefficients, is checked at two step sizes.
  const dense=s.jacobian(x,{sparse:false});
  for(let col=0;col<s.n;col++) for(const h of [2e-6,6e-7]) {
    const p=x.slice(),m=x.slice();p[col]+=h;m[col]-=h;
    const rp=s.residual(p),rm=s.residual(m);
    for(const row of terminal){const fd=(rp[row]-rm[row])/(2*h),v=dense[row*s.n+col];
      assert.ok(Math.abs(fd-v)/Math.max(1,Math.abs(fd),Math.abs(v))<3e-7,JSON.stringify({row,col,h,fd,v}));}
  }
  const before=s.evaluate(x),rebased=s.rebase(x);assert.ok(max(s.residual(rebased).map((v,i)=>v-before.residual[i]))<1e-11);
  const r=solveCoupledStreamtubeBody(s,{initial:rebased,maxIterations:10,tolerance:1e-10});
  assert.equal(r.converged,true,r.reason);assert.ok(directStreamtubeVolumeGeometry(r.flow.nodes).valid);
  assert.equal(r.boundaryLayer.surfaces.length,4);assert.equal(r.boundaryLayer.wakes.length,2);
  assert.ok(max(bankAngles(s,r.flow))<2e-10);
  for(let g=0;g<r.flow.nodes.length;g++) {
    const c=directChannelConservation({nodes:r.flow.nodes[g],sections:r.flow.sections.map(row=>row[g]),cells:r.flow.cells.map(row=>row[g])});
    for(const key of ['total','maxLocal','internalCancellation']) for(const k of [0,3])assert.ok(Math.abs(c[key][k])<2e-9);
  }
  const restored=createCoupledStreamtubeBody(input,{...options,initialEuler:r.flow,initialBL:r.x.slice(s.ne)});
  assert.ok(max(restored.residual(restored.initial))<1e-10);
  t.diagnostic(JSON.stringify({unknowns:s.n,iterations:r.history.length-1,families:r.families,angles:bankAngles(s,r.flow)}));
});
