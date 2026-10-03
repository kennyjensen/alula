// SPDX-License-Identifier: GPL-2.0-or-later
import { normInf, solveLinear } from '../numerics/linear.js';
import { conserved, totalConditions } from './gas.js';
import { sparseAdd } from '../numerics/sparse.js';

// Pseudo-time regularizes Newton corrections to the steady conservative
// residual. It never changes the residual used to certify convergence.
// Unknowns are [ln(rho/rho_inf), u, v, (p-p_inf)/(rho_inf U_inf^2)].
export function temporalJacobian(states,reference,cfl){
  const speed=Math.hypot(reference.u,reference.v),q2=reference.rho*speed**2;
  const scales=[reference.rho*speed,q2,q2,reference.rho*speed*totalConditions(reference).h0];
  return states.map(s=>{
    const wave=Math.hypot(s.u,s.v)+Math.sqrt(reference.gamma*s.p/s.rho);
    const derivative=[s.rho,0,0,0,
      s.rho*s.u,s.rho,0,0,
      s.rho*s.v,0,s.rho,0,
      .5*s.rho*(s.u*s.u+s.v*s.v),s.rho*s.u,s.rho*s.v,q2/(reference.gamma-1)];
    return Float64Array.from(derivative,(v,i)=>v*wave/(cfl*scales[Math.floor(i/4)]));
  });
}

export function solvePseudoTransient({initial,residual,jacobian,decode,reference,
  admissible,tolerance=1e-9,maxIterations=200,onIteration,initialCfl=1,linearSolve=solveLinear}){
  if(!(initialCfl>0)||!Number.isFinite(initialCfl))throw new Error('Invalid pseudo-time CFL.');
  if(!(tolerance>0)||!Number.isFinite(tolerance)||!Number.isInteger(maxIterations)||maxIterations<0)throw new Error('Invalid pseudo-time convergence controls.');
  if(!initial.length||!Array.from(initial).every(Number.isFinite)||!admissible(initial))throw new Error('Invalid pseudo-time initial state.');
  let x=Float64Array.from(initial),r=residual(x),cfl=initialCfl;
  const history=[],norm2=v=>{let scale=0,sum=0;for(const x of v){const a=Math.abs(x);if(a>scale){sum=1+sum*(scale/a)**2;scale=a;}else if(a>0)sum+=(a/scale)**2;}return scale*Math.sqrt(sum);};
  const report=(iteration,step)=>{const h={iteration,residual:normInf(r),step,cfl};history.push(h);onIteration?.(h);};
  report(0,0);
  for(let iteration=0;iteration<maxIterations;iteration++){
    if(normInf(r)<=tolerance)return{converged:true,x,history,reason:'residual'};
    const base=jacobian(x),n=x.length,states=decode(x),oldNorm=norm2(r);
    const oldConserved=states.map(s=>conserved(s,reference.gamma));
    const speed=Math.hypot(reference.u,reference.v),mass=reference.rho*speed;
    const scales=[mass,mass*speed,mass*speed,mass*totalConditions(reference).h0];
    let accepted=false;
    for(let attempt=0;attempt<6&&!accepted;attempt++){
      const sparse=base.rowPtr!==undefined,matrix=sparse?{...base,values:base.values.slice()}:base.slice(),blocks=temporalJacobian(states,reference,cfl);
      for(let cell=0;cell<states.length;cell++)for(let row=0;row<4;row++)for(let col=0;col<4;col++)
        if(sparse)sparseAdd(matrix,4*cell+row,4*cell+col,blocks[cell][4*row+col]);
        else matrix[(4*cell+row)*n+4*cell+col]+=blocks[cell][4*row+col];
      let delta;
      try{delta=linearSolve(matrix,r.map(v=>-v));}catch{cfl/=10;continue;}
      for(let step=1;step>=2**-14;step*=.5){
        const candidate=x.map((v,i)=>v+step*delta[i]);
        if(!candidate.every(Number.isFinite)||!admissible(candidate))continue;
        let next;try{next=residual(candidate);}catch{continue;}
        const trialStates=decode(candidate),timeResidual=next.slice();
        for(let cell=0;cell<states.length;cell++){
          const state=states[cell],wave=Math.hypot(state.u,state.v)+Math.sqrt(reference.gamma*state.p/state.rho);
          const q=conserved(trialStates[cell],reference.gamma);
          for(let k=0;k<4;k++)timeResidual[4*cell+k]+=wave*(q[k]-oldConserved[cell][k])/(cfl*scales[k]);
        }
        if(next.every(Number.isFinite)&&norm2(timeResidual)<(1-1e-4*step)*oldNorm){
          x=candidate;r=next;accepted=true;report(iteration+1,step);
          cfl=Math.min(1e12,cfl*(step===1?1.5:Math.max(.5,step)));
          break;
        }
      }
      if(!accepted)cfl/=10;
    }
    if(!accepted)return{converged:false,x,history,reason:'pseudo-time line search failed'};
  }
  return{converged:normInf(r)<=tolerance,x,history,reason:normInf(r)<=tolerance?'residual':'iteration limit'};
}
