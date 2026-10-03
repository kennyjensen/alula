// SPDX-License-Identifier: GPL-2.0-or-later
// Subcritical, uniform-entropy and uniform-total-enthalpy gas relation. This
// is the thermodynamic closure for the MSES ISMOM=2/4 subcritical limit, not a
// pressure correction applied after an incompressible flow calculation.
export function isentropicState(u,v,{mach=0,gamma=1.4}={}){
  if(![u,v,mach,gamma].every(Number.isFinite)||mach<0||mach>=1||gamma<=1)throw new Error('Invalid isentropic state.');
  const speedSquared=u*u+v*v,m2=mach*mach,gm1=gamma-1;
  const increment=.5*gm1*m2*(1-speedSquared),temperatureRatio=1+increment;
  if(!(temperatureRatio>0))throw new Error('Velocity exceeds the available stagnation enthalpy.');
  const logTemperature=Math.log1p(increment),rho=Math.exp(logTemperature/gm1);
  const localMachSquared=m2*speedSquared/temperatureRatio;
  if(localMachSquared>=1)throw new Error('Local sonic flow is outside the subcritical potential formulation.');
  const cp=m2===0?1-speedSquared:2*Math.expm1(gamma/gm1*logTemperature)/(gamma*m2);
  const rhoSpeedSquared=-.5*m2*rho/temperatureRatio;
  return{rho,cp,temperatureRatio,machSquared:localMachSquared,rhoSpeedSquared,cpSpeedSquared:-rho};
}

export function isentropicMassFlux(u,v,nx,ny,conditions){
  if(!Number.isFinite(nx)||!Number.isFinite(ny)||Math.abs(Math.hypot(nx,ny)-1)>1e-10)throw new Error('Mass flux needs a unit normal.');
  const state=isentropicState(u,v,conditions),normal=u*nx+v*ny;
  return{...state,flux:state.rho*normal,
    derivative:[state.rho*nx+2*state.rhoSpeedSquared*u*normal,state.rho*ny+2*state.rhoSpeedSquared*v*normal]};
}
