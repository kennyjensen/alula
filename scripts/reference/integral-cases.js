// SPDX-License-Identifier: GPL-2.0-or-later
const laminar = [{ s: .04, aux: .3, theta: .00016, deltaStar: .00042, ue: 1.12 },
  { s: .055, aux: .6, theta: .0002, deltaStar: .00056, ue: 1.08 }];
const turbulent = [{ s: .45, aux: .04, theta: .0011, deltaStar: .0018, ue: 1.14 },
  { s: .5, aux: .043, theta: .0013, deltaStar: .0022, ue: 1.10 }];
const wake = [{ s: 1.1, aux: .035, theta: .003, deltaStar: .005, ue: .94 },
  { s: 1.3, aux: .03, theta: .0029, deltaStar: .0046, ue: .96 }];
const transition = [{ s: .23, aux: 6.5, theta: .00045, deltaStar: .0013, ue: 1.12 },
  { s: .27, aux: .04, theta: .00056, deltaStar: .0012, ue: 1.1 }];
export const integralCases = [];
for (const velocityConvention of ['physical','xfoil']) for (const mach of [0,.4]) {
  for (const [regime, states] of [['laminar',laminar], ['turbulent',turbulent], ['wake',wake]]) {
    integralCases.push({ name: `${regime}-${velocityConvention}-${mach}`, regime, states, parameters: { reynolds: 1e6, mach, velocityConvention } });
  }
  const similarity = { s: .001, aux: 0, theta: .00008, deltaStar: .00018, ue: .12 };
  integralCases.push({ name: `similarity-${velocityConvention}-${mach}`, regime: 'similarity', states: [similarity,similarity], parameters: { reynolds: 1e6, mach, velocityConvention } });
  integralCases.push({ name: `forced-transition-${velocityConvention}-${mach}`, regime: 'transition', states: transition, tripS: .25, parameters: { reynolds: 1e6, mach, velocityConvention } });
}
integralCases.push({ name: 'free-transition-physical', regime: 'transition', states: [
  { s: .2, aux: 8.8, theta: .0007, deltaStar: .0024, ue: 1.2 },
  { s: .24, aux: .035, theta: .0009, deltaStar: .002, ue: 1.15 }], parameters: { reynolds: 1e6, mach: .3, velocityConvention: 'physical' } });
integralCases.push({ name: 'trailing-edge', regime: 'te', states: [wake[0],wake[1]], matched: { aux: .032, theta: .0028, deltaStar: .0048 }, parameters: { reynolds: 1e6, mach: .3, velocityConvention: 'physical' } });

// Finite-base prerequisites: native wake equations use total displacement
// thickness in the mass deficit and subtract the prescribed gap for shape
// closures. These local blocks do not implement finite-base outer geometry.
for(const mach of [0,.3]){
  const parameters={reynolds:1e6,mach,velocityConvention:'physical'};
  integralCases.push({name:`finite-gap-wake-${mach}`,regime:'wake',parameters,states:[
    {s:1.02,aux:.028,theta:.0012,deltaStar:.0029,wakeGap:.0006,ue:1.01},
    {s:1.035,aux:.0285,theta:.00124,deltaStar:.0027,wakeGap:.0002,ue:1.004}]});
  integralCases.push({name:`finite-base-trailing-edge-${mach}`,regime:'te',parameters,
    states:[wake[0],{s:1.1,aux:.034,theta:.0029,deltaStar:.0053,wakeGap:.0004,ue:1.02}],
    gap:.0004,surfaceStates:[
      {aux:.025,theta:.0012,deltaStar:.0022},
      {aux:.04,theta:.0016,deltaStar:.0026}],
    // Independently specified totals supplied directly to original TESYS.
    matched:{aux:.03357142857142857,theta:.0028,deltaStar:.0052}});
}
