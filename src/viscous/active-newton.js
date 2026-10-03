// SPDX-License-Identifier: GPL-2.0-or-later
import { normInf, solveLinear } from '../numerics/linear.js';

// Newton with discrete transition updates. A trial owns a copy of both the
// state and the active set; rejected trials must restore the previous set.
export function solveActiveNewton({initial,residual,jacobian,admissible,updateActive,
  snapshot,restore,linearSolve=solveLinear,trialState=(x,d,step)=>x.map((v,i)=>v+step*d[i]),
  limitedVariables=x=>({offset:0,length:x.length}),maxIterations=40,tolerance=1e-8,eventStepFraction=1,onIteration}){
  if(!Number.isInteger(maxIterations)||maxIterations<0||!Number.isFinite(tolerance)||tolerance<=0)throw new Error('Invalid Newton convergence controls.');
  if(!Number.isFinite(eventStepFraction)||eventStepFraction<=0||eventStepFraction>1)throw new Error('Invalid active-event step fraction.');
  let x=Float64Array.from(initial),history=[];
  if(!admissible(x))throw new Error('Invalid Newton initial state.');
  updateActive(x);
  let r=residual(x);
  const report=(iteration,step,activeChange=false)=>{const h={iteration,residual:normInf(r),step,activeChange};history.push(h);onIteration?.(h);};
  report(0,0);
  for(let iteration=0;iteration<maxIterations;iteration++){
    if(normInf(r)<=tolerance)return{converged:true,x,history,reason:'residual'};
    let delta;
    try{delta=linearSolve(jacobian(x),r.map(v=>-v),x);}
    catch(error){return{converged:false,x,history,reason:error.message};}
    const saved=snapshot();let accepted=false;
    // Limit physical thickness reductions and velocity changes before the
    // residual line search. Surface velocities may cross zero at stagnation.
    let maximum=1;
    const limits=limitedVariables(x);
    for(let i=limits.offset;i<limits.offset+limits.length;i++){
      const k=(i-limits.offset)%4;
      if((k===1||k===2)&&delta[i]<0)maximum=Math.min(maximum,-.5*x[i]/delta[i]);
      if(k===3&&Math.abs(delta[i])>.2)maximum=Math.min(maximum,.2/Math.abs(delta[i]));
    }
    for(let step=maximum;step>=maximum*2**-20;step*=.5){
      restore(saved);
      const candidate=trialState(x,delta,step);
      if(!admissible(candidate))continue;
      let next,changed;
      try{
        changed=updateActive(candidate);if(!admissible(candidate))continue;
        // Limit an event's travel into a different equation branch. Accepting
        // a full step merely because the active set changed can cycle between
        // two transition intervals. Ordinary descent steps retain full length.
        if(changed&&step>maximum*eventStepFraction*(1+1e-12))continue;
        next=residual(candidate);if(!next.every(Number.isFinite))continue;
      }
      catch{continue;}
      const decrease=normInf(next)<=(1-1e-4*step)*normInf(r);
      // Changing the stagnation/transition interval changes the equations and
      // initializes different variables. Their residual norm is not the old
      // merit function. Accept the admissible, physically limited event step,
      // then resume ordinary descent on the new active set.
      if(decrease||changed){
        x=candidate;r=next;accepted=true;report(iteration+1,step,changed);break;
      }
    }
    if(!accepted){restore(saved);return{converged:false,x,history,reason:'line search failed'};}
    const recent=history.slice(-3);
    if(normInf(r)>tolerance&&recent.length===3&&recent.every(h=>h.iteration>0&&!h.activeChange&&h.step<1e-6)
      &&recent.at(-1).residual>.9999*recent[0].residual)return{converged:false,x,history,reason:'stagnated Newton step'};
  }
  return{converged:normInf(r)<=tolerance,x,history,reason:normInf(r)<=tolerance?'residual':'iteration limit'};
}
