// SPDX-License-Identifier: GPL-2.0-or-later
import { solveSparse } from '../numerics/sparse.js';
import { solveSparseDirect } from '../numerics/klu.js';

// Difficult stagnation/transition blocks may need a larger Krylov space.
// Reuse the first attempt's correction and strengthen preconditioning; both
// attempts certify the original Jacobian at the same requested tolerance.
export function solvePotentialLinear(a,b,{preconditionerMatrix,linearTolerance=1e-7,fillLevel=1,linearBackend='klu'}={}){
  if(linearBackend==='klu')return solveSparseDirect(a,b,{tolerance:linearTolerance}).x;
  if(linearBackend!=='gmres')throw new Error('Unknown potential linear solver.');
  const options={preconditionerMatrix,tolerance:linearTolerance,fillLevel,restart:150};
  try{return solveSparse(a,b,options).x;}
  catch(error){
    if(error.code!=='GMRES_ITERATION_LIMIT')throw error;
    return solveSparse(a,b,{...options,initial:error.solution,fillLevel:Math.min(3,fillLevel+1),restart:300,maxIterations:2000}).x;
  }
}
