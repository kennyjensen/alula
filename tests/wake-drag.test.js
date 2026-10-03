import test from 'node:test';
import assert from 'node:assert/strict';
import { wakeMomentumDrag } from '../src/viscous/wake-drag.js';

test('compressible wake drag agrees with independent integration of the momentum equation',()=>{
  for(const mach of [0,.4,.7])for(const ue of [.6,.9,1.05]){
    const theta=.004,h=1.7,gamma=1.4,start=Math.log(ue),hInfinity=1+(gamma-1)*mach**2;
    // Integrate d(theta)/d(log U), retaining the local Mach number, with RK4.
    // This does not use the density or power expression under test.
    const slope=(t,value)=>{
      const u=Math.exp(t),temperature=1+.5*(gamma-1)*mach**2*(1-u*u);
      const localMachSquared=mach**2*u*u/temperature,H=hInfinity+(h-hInfinity)*t/start;
      return -(H+2-localMachSquared)*value;
    };
    let y=theta;const step=-start/1024;
    for(let i=0;i<1024;i++){
      const t=start+i*step,k1=slope(t,y),k2=slope(t+step/2,y+step*k1/2),k3=slope(t+step/2,y+step*k2/2),k4=slope(t+step,y+step*k3);
      y+=step*(k1+2*k2+2*k3+k4)/6;
    }
    const result=wakeMomentumDrag({theta,deltaStar:h*theta,ue},{mach,gamma});
    assert.ok(Math.abs(result.thetaInfinity-y)/y<1e-11);
    assert.equal(result.cd,2*result.thetaInfinity);
    if(mach===0)assert.equal(result.cd,2*theta*ue**((5+h)/2));
  }
});

test('far-wake compressible shape factor follows mass-deficit and momentum-deficit definitions',()=>{
  for(const mach of [0,.2,.7]){
    const gamma=1.4,epsilon=1e-6;let displacement=0,momentum=0;
    for(let j=0;j<1000;j++){
      const y=(j+.5)/100,u=1-epsilon*Math.exp(-y*y);
      // Across an adiabatic wake at constant pressure, rho/rho_e = T_e/T.
      const rho=1/(1+.5*(gamma-1)*mach**2*(1-u*u));
      displacement+=1-rho*u;momentum+=rho*u*(1-u);
    }
    const result=wakeMomentumDrag({theta:.001,deltaStar:.0012,ue:1},{mach,gamma});
    assert.ok(Math.abs(displacement/momentum-result.hInfinity)<1e-6);
    assert.equal(result.cd,.002);
  }
});
