// SPDX-License-Identifier: GPL-2.0-or-later
import { isentropicState } from '../potential/isentropic.js';

// Squire–Young extrapolation of a compressible, adiabatic thin wake.
// With zero wall shear: d log(theta) = -(H + 2 - M_e^2) d log(Ue).
// Isentropic edge flow gives d log(rho_e) = -M_e^2 d log(Ue).
// Integrate H linearly in log(Ue), from the last station to freestream.
// A weak velocity deficit at constant pressure/total enthalpy has
// H_infinity = 1 + (gamma - 1) M_infinity^2, not 1 at finite Mach.
// The extrapolation is approximate; demonstrate wake-length independence.
// Drela, TASOPT 2.16, equations 338–343 (momentum-area convention).
export function wakeMomentumDrag({theta,deltaStar,ue},{mach=0,gamma=1.4}={}){
  if(![theta,deltaStar,ue].every(Number.isFinite)||theta<=0||deltaStar<=theta||ue<=0)throw new Error('Invalid downstream wake state.');
  const {rho}=isentropicState(ue,0,{mach,gamma});
  const h=deltaStar/theta,hInfinity=1+(gamma-1)*mach*mach;
  const thetaInfinity=rho*theta*ue**(2+.5*(h+hInfinity));
  return{cd:2*thetaInfinity,thetaInfinity,hInfinity,rho,edgeVelocity:ue};
}
