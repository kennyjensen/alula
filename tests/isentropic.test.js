import test from 'node:test';
import assert from 'node:assert/strict';
import { isentropicState,isentropicMassFlux } from '../src/potential/isentropic.js';
import { totalConditions,freestream } from '../src/euler/gas.js';

test('subcritical gas closure preserves entropy and total enthalpy and has a regular Mach-zero limit',()=>{
  for(const mach of [0,1e-8,.1,.4])for(const [u,v] of [[.2,.1],[1,0],[1.2,.3]]){
    const s=isentropicState(u,v,{mach});
    if(mach>1e-6){
      const reference=freestream({mach}),flow={rho:s.rho,u,v,p:reference.p+.5*s.cp};
      const a=totalConditions(flow),b=totalConditions(reference);
      assert.ok(Math.abs(a.entropy-b.entropy)<1e-14);assert.ok(Math.abs(a.h0/b.h0-1)<1e-14);
    }else{assert.ok(Math.abs(s.cp-(1-u*u-v*v))<1e-14);assert.ok(Math.abs(s.rho-1)<1e-14);}
  }
});

test('isentropic mass-flux derivatives match differences and retain ellipticity up to the sonic limit',()=>{
  for(const mach of [0,.2,.6]){
    const u=1.1,v=.3,h=1e-6,nx=.6,ny=.8,conditions={mach},s=isentropicMassFlux(u,v,nx,ny,conditions);
    for(let k=0;k<2;k++){
      const a=isentropicMassFlux(u+(k===0?h:0),v+(k===1?h:0),nx,ny,conditions);
      const b=isentropicMassFlux(u-(k===0?h:0),v-(k===1?h:0),nx,ny,conditions);
      assert.ok(Math.abs((a.flux-b.flux)/(2*h)-s.derivative[k])<2e-10);
      assert.ok(Math.abs((a.cp-b.cp)/(2*h)+2*s.rho*(k===0?u:v))<3e-10);
    }
    const x=isentropicMassFlux(u,v,1,0,conditions),y=isentropicMassFlux(u,v,0,1,conditions);
    const determinant=x.derivative[0]*y.derivative[1]-x.derivative[1]*y.derivative[0];
    assert.ok(determinant>0);assert.ok(Math.abs(determinant-s.rho**2*(1-s.machSquared))<1e-14);
  }
  assert.throws(()=>isentropicState(2,0,{mach:.6}),/sonic/);
});
