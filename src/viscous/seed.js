// SPDX-License-Identifier: GPL-2.0-or-later
import { isentropicState } from '../potential/isentropic.js';

// Change only an INITIAL guess between Mach numbers. Preserve the native
// kinematic shape factor Hk, since an almost mixed incompressible wake has
// H near one and is inadmissible under a finite-Mach closure without this
// conversion. Every subsequent residual uses the requested physical state.
export function remapInitialMach(x,fromMach,toMach){
  const result=Float64Array.from(x);
  for(let i=0;i<x.length;i+=4){
    const oldM2=isentropicState(x[i+3],0,{mach:fromMach}).machSquared;
    const newM2=isentropicState(x[i+3],0,{mach:toMach}).machSquared;
    const hk=(x[i+2]/x[i+1]-.29*oldM2)/(1+.113*oldM2);
    result[i+2]=x[i+1]*(hk*(1+.113*newM2)+.29*newM2);
  }
  return result;
}

export function remapSeedMach(seed,mach){
  const key=JSON.parse(seed.key),x=remapInitialMach(seed.x,key.mach,mach);key.mach=mach;
  return{...seed,key:JSON.stringify(key),x};
}
