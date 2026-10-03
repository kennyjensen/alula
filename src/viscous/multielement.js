// SPDX-License-Identifier: GPL-2.0-or-later
import { solveCoupledAssembly } from './assembly.js';

// Close the wake centerline using its mean velocity, not either side of its
// source-sheet jump. All BL/mass-displacement unknowns are solved together for
// each geometry. Success requires both the BL and wake-shape residuals.
export function solveMultielementViscous(input,{onIteration,maxIterations=40,
  tolerance=1e-8,wakeTolerance=1e-6,maxWakeIterations=40,wakeRelaxation=.5}={}){
  if(!Number.isInteger(maxWakeIterations)||maxWakeIterations<0||!Number.isFinite(wakeTolerance)||wakeTolerance<=0
    ||!Number.isFinite(wakeRelaxation)||wakeRelaxation<=0||wakeRelaxation>1)throw new Error('Invalid wake convergence controls.');
  let result,wakePaths=input.wakePaths,initialState,activeSet,referenceLengths,seed=input.seed;
  const wakeHistory=[],history=[];
  for(let iteration=0;iteration<=maxWakeIterations;iteration++){
    result=solveCoupledAssembly({initialization:input.elements.length>1?'auto':'native',wakeInitialization:'inviscid',...input,wakePaths,seed},{maxIterations,tolerance,initialState,activeSet,
      onIteration:h=>onIteration?.({...h,wakeIteration:iteration})});
    history.push(...result.history.map(h=>({...h,wakeIteration:iteration})));
    result={...result,history};
    if(!result.converged)return{...result,wakeHistory,wakeConverged:false,wakeResidual:Infinity};
    const {outer}=result.system;
    if(!referenceLengths)referenceLengths=outer.wakes.map(w=>w.segments.map(p=>p.length));
    const mass=Float64Array.from(result.states,(s,i)=>result.x[4*i+3]*s.deltaStar);
    const velocity=outer.velocityField(mass);
    let wakeResidual=0;
    const directions=outer.wakes.map(w=>w.segments.map(p=>{
      const v=velocity(p),normal=v.u*p.nx+v.v*p.ny,speed=Math.hypot(v.u,v.v);
      if(!(v.u*p.tx+v.v*p.ty>0))throw new Error('Wake flow reverses; the thin-wake model is inadmissible.');
      wakeResidual=Math.max(wakeResidual,Math.abs(normal));
      return{tx:v.u/speed,ty:v.v/speed};
    }));
    wakeHistory.push({iteration,residual:wakeResidual,blResidual:result.history.at(-1).residual,
      iterations:result.history.filter(h=>h.wakeIteration===iteration&&h.iteration>0).length});
    onIteration?.({wakeIteration:iteration,wakeResidual,kind:'wake'});
    if(wakeResidual<=wakeTolerance)return{...result,wakeHistory,wakeConverged:true,wakeResidual};
    if(iteration===maxWakeIterations)return{...result,converged:false,reason:'wake iteration limit',wakeHistory,wakeConverged:false,wakeResidual};
    wakePaths=outer.wakes.map((w,e)=>{
      const points=[{...w.points[0]}];
      for(let j=0;j<w.segments.length;j++){
        const p=w.segments[j],d=directions[e][j],tx=(1-wakeRelaxation)*p.tx+wakeRelaxation*d.tx,ty=(1-wakeRelaxation)*p.ty+wakeRelaxation*d.ty;
        const scale=referenceLengths[e][j]/Math.hypot(tx,ty),a=points.at(-1);
        points.push({x:a.x+scale*tx,y:a.y+scale*ty});
      }
      return points;
    });
    seed=result.system.exportSeed(result.x);
    // Finite-base seed import restores the prescribed gap on the new wake
    // geometry while retaining fluid displacement. Do not overwrite it with
    // old total thicknesses after that coordinate transfer.
    initialState=seed.wakeGaps?undefined:result.x;activeSet=result.system.snapshot();
  }
  throw new Error('Unreachable wake iteration state.');
}
