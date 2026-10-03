import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';
import { evaluateLeadingTransitionInterval } from '../src/viscous/leading-transition-interval.js';

const fixture=JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/hk-floor.json',import.meta.url)));
const close=(a,b)=>assert.ok(Math.abs(a-b)<=3e-9*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);

test('RAE wake-floor contact has distinct one-sided derivatives within the physical domain', () => {
 const c=JSON.parse(fs.readFileSync(new URL('./fixtures/rae8x11-wake-floor-interval.json',import.meta.url)));
 assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${c.provenance.checkpoint}`,import.meta.url))).digest('hex'),c.provenance.sha256);
 const kernel=createIntegralKernel(c.parameters),base=kernel.interval(c.input);
 for(const side of ['upstream','downstream']){
  const state=c.input[side],h=state.theta*1e-7;
  assert.ok(Math.abs(kernel.station(state,'wake').rawHk-1.00005)<1e-14);
  const plus={...state,deltaStar:state.deltaStar+h},minus={...state,deltaStar:state.deltaStar-h};
  assert.ok(kernel.station(minus,'wake').rawHk>1,'Both perturbations must stay physically admissible.');
  const p=kernel.interval({...c.input,[side]:plus}),m=kernel.interval({...c.input,[side]:minus});
  // The energy row exposes the floor without the shear closure's extreme
  // curvature near Hk=1. The contact Jacobian intentionally uses the
  // increasing-Hk branch; a centered difference straddles two branches.
  const analytic=base[side][2][2],forward=(p.residual[2]-base.residual[2])/h;
  const backward=(base.residual[2]-m.residual[2])/h;
  assert.ok(Math.abs(forward-analytic)<1e-5*Math.abs(analytic));
  assert.ok(Math.abs(backward)<.001*Math.abs(analytic));
  assert.ok(Math.abs(.5*(forward+backward)-analytic)>.4*Math.abs(analytic));
 }
});

test('clipped Hk blocks preserve the original Fortran residuals and native partials',()=>{
 for(const [file,hash] of Object.entries(fixture.provenance.sha256))assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${file}`,import.meta.url))).digest('hex'),hash,file);
 for(const c of fixture.cases){const native=createIntegralKernel({...c.parameters,exactJacobian:false}).interval(c.input);
  for(const key of ['residual','upstream','downstream'])native[key].flat().forEach((v,i)=>close(v,c.expected[key].flat()[i]));
  const exact=createIntegralKernel({...c.parameters,exactJacobian:true}).interval(c.input), continued=createIntegralKernel({...c.parameters,exactJacobian:true,hkFloorLinearization:'native'}).interval(c.input);
  assert.deepEqual(continued.residual,exact.residual);assert.deepEqual(continued.properties,exact.properties);
  assert.equal(continued.hkFloorLinearizationUsed,true);assert.notDeepEqual(continued.downstream,exact.downstream);
  // Momentum and energy rows have no shear-lag ReTheta approximation.
  // These surface cases also avoid the separate extreme-upwind clamp.
  for(const key of ['upstream','downstream'])for(let r=1;r<3;r++)continued[key][r].forEach((v,i)=>close(v,c.expected[key][r][i]));
 }
});

test('mixed transition adds only native Hk-floor continuation to resolved-root partials',()=>{
 const c=fixture.cases[0],parameters={...c.parameters,exactJacobian:true,transitionTolerance:1e-12};
 const exactKernel=createIntegralKernel(parameters),nativeKernel=createIntegralKernel({...parameters,hkFloorLinearization:'native'});
 const exact=evaluateTransitionInterval(exactKernel,c.input,{jacobian:true}),continued=evaluateTransitionInterval(nativeKernel,c.input,{jacobian:true});
 assert.deepEqual(continued.residual,exact.residual);assert.deepEqual(continued.transition,exact.transition);assert.deepEqual(continued.partials.location,exact.partials.location);
 const correctedInput={...c.input,tripS:continued.transition.forced?c.input.tripS:Number.MAX_VALUE};
 const nativeAnalytic=nativeKernel.interval(correctedInput),exactAnalytic=exactKernel.interval(correctedInput);
 for(const side of ['upstream','downstream'])for(let r=0;r<3;r++)for(let k=0;k<5;k++)assert.equal(continued.partials[side][r][k],exact.partials[side][r][k]+(nativeAnalytic[side][r][k]-exactAnalytic[side][r][k]));
 for(let r=0;r<3;r++)assert.equal(continued.partials.trip[r],exact.partials.trip[r]+(nativeAnalytic.tripDerivative[r]-exactAnalytic.tripDerivative[r]));
 assert.notDeepEqual(continued.partials.downstream,exact.partials.downstream);
});

test('unclipped blocks remain byte-identical and omitted policy keeps exact behavior',()=>{
 const c=fixture.cases[2],input={...c.input,upstream:{...c.input.upstream,deltaStar:3*c.input.upstream.theta},downstream:{...c.input.downstream,deltaStar:3*c.input.downstream.theta}};
 const p={...c.parameters,exactJacobian:true},a=createIntegralKernel(p),b=createIntegralKernel({...p,hkFloorLinearization:'exact'}),n=createIntegralKernel({...p,hkFloorLinearization:'native'});
 assert.deepEqual(a.parameters,b.parameters);assert.deepEqual(a.interval(input),b.interval(input));assert.deepEqual(a.interval(input),n.interval(input));
 assert.equal(n.parameters.hkFloorLinearization,'native');assert.equal(a.parameters.hkFloorLinearization,undefined);
 assert.throws(()=>createIntegralKernel({hkFloorLinearization:'floor-off'}),/Invalid Hk-floor/);
});

test('unclipped mixed intervals retain every resolved-root partial unchanged',()=>{
 const original=JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/integral.json',import.meta.url)));
 let checked=0;
 for(const {input:c} of original.cases.filter(c=>c.input.regime==='transition')){
  const p={...c.parameters,exactJacobian:true,transitionTolerance:1e-12},input={upstream:c.states[0],downstream:c.states[1],regime:'transition',tripS:c.tripS};
  const exact=evaluateTransitionInterval(createIntegralKernel(p),input,{jacobian:true}),native=evaluateTransitionInterval(createIntegralKernel({...p,hkFloorLinearization:'native'}),input,{jacobian:true});
  if(native.hkFloorLinearizationUsed)continue;
  checked++;assert.deepEqual(native,exact);
 }
 assert.ok(checked>=2,'Check natural and forced unclipped mixed intervals.');
});

test('native floor policy keeps the independent exact shear-lag correction active',()=>{
 const c=fixture.cases[1],p={...c.parameters,exactJacobian:true,hkFloorLinearization:'native'};
 const continued=createIntegralKernel(p).interval(c.input),fullyNative=createIntegralKernel({...p,exactJacobian:false}).interval(c.input);
 assert.deepEqual(continued.residual,fullyNative.residual);
 assert.notDeepEqual(continued.downstream[0],fullyNative.downstream[0]);
 // Internal baseline re-evaluation must not leave the context in exact-floor mode.
 const kernel=createIntegralKernel(p);kernel.interval(c.input,{hkFloorCorrection:true});assert.deepEqual(kernel.interval(c.input),continued);
});

test('forced interior and terminal transitions preserve roots and add only the floor extension',()=>{
 const c=fixture.cases[0],p={...c.parameters,ncrit:100,exactJacobian:true,transitionTolerance:1e-12};
 for(const f of [.3,1]){
  const input={...c.input,tripS:c.input.upstream.s+f*(c.input.downstream.s-c.input.upstream.s)};
  const exactKernel=createIntegralKernel(p),nativeKernel=createIntegralKernel({...p,hkFloorLinearization:'native'});
  const exact=evaluateTransitionInterval(exactKernel,input,{jacobian:true}),native=evaluateTransitionInterval(nativeKernel,input,{jacobian:true});
  assert.equal(native.transition.forced,true);assert.deepEqual(native.residual,exact.residual);assert.deepEqual(native.transition,exact.transition);assert.deepEqual(native.partials.location,exact.partials.location);
  const correction=nativeKernel.interval(input,{hkFloorCorrection:true}).hkFloorJacobianCorrection;assert.ok(correction);
  for(const side of ['upstream','downstream'])for(let r=0;r<3;r++)for(let k=0;k<5;k++)assert.equal(native.partials[side][r][k],exact.partials[side][r][k]+correction[side][r][k]);
  for(let r=0;r<3;r++)assert.equal(native.partials.trip[r],exact.partials.trip[r]+correction.trip[r]);
 }
});

test('a clipped station or correction evaluation cannot leak its branch into the next interval',()=>{
 const c=fixture.cases[2],p={...c.parameters,exactJacobian:true},kernel=createIntegralKernel({...p,hkFloorLinearization:'native'});
 const input={...c.input,upstream:{...c.input.upstream,deltaStar:3*c.input.upstream.theta},downstream:{...c.input.downstream,deltaStar:3*c.input.downstream.theta}},expected=createIntegralKernel(p).interval(input);
 kernel.station(c.input.downstream);assert.deepEqual(kernel.interval(input),expected);
 kernel.interval(c.input,{hkFloorCorrection:true});assert.deepEqual(kernel.interval(input),expected);
 assert.throws(()=>kernel.interval({...c.input,downstream:{...c.input.downstream,theta:-1}},{hkFloorCorrection:true}),/Inadmissible/);assert.deepEqual(kernel.interval(input),expected);
});

test('virtual leading-transition composite retains its exact partials under either floor policy',()=>{
 const c=fixture.cases[2],p={...c.parameters,ncrit:100,exactJacobian:true,transitionTolerance:1e-12};
 const input={downstream:{...c.input.downstream,aux:.03},tripS:.05};
 const exact=evaluateLeadingTransitionInterval(createIntegralKernel(p),input,{jacobian:true});
 const native=evaluateLeadingTransitionInterval(createIntegralKernel({...p,hkFloorLinearization:'native'}),input,{jacobian:true});
 assert.ok(native.properties.rawHk<1.05);assert.equal(native.transition.forced,true);
 for(const key of ['residual','transition','partials','upstreamState'])assert.deepEqual(native[key],exact[key]);
});
