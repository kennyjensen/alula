// SPDX-License-Identifier: GPL-2.0-or-later
import { freestream, totalConditions } from './gas.js';
import { eulerResidual } from './residual.js';
import { solveNewton } from '../numerics/newton.js';
import { eulerJacobian } from './jacobian.js';
import { linearReconstruction } from './reconstruction.js';
import { solvePseudoTransient } from './pseudo-transient.js';
import { solveSparse,sparseIndex,sparseAdd } from '../numerics/sparse.js';
import { solveSparseDirect } from '../numerics/klu.js';

// Density/pressure remain positive under damping. Dynamic-pressure scaling of
// the pressure unknown avoids a large ambient pressure dominating Newton steps.
export function encodeStates(states, reference) {
  const q2 = reference.rho * (reference.u ** 2 + reference.v ** 2);
  return Float64Array.from(states.flatMap(s => [Math.log(s.rho / reference.rho), s.u, s.v, (s.p - reference.p) / q2]));
}
export function decodeStates(x, count, reference) {
  const q2 = reference.rho * (reference.u ** 2 + reference.v ** 2);
  return Array.from({ length: count }, (_, i) => ({ rho: reference.rho * Math.exp(x[4 * i]),
    u: x[4 * i + 1], v: x[4 * i + 2], p: reference.p + q2 * x[4 * i + 3] }));
}

export function solveEuler(mesh, { mach = 0.3, alpha = 0, gamma = 1.4, initial,
  tolerance = 1e-9, maxIterations = 30, onIteration,spatialOrder=1,pseudoTransient=false,linearBackend='dense' } = {}) {
  if(![1,2].includes(spatialOrder))throw new Error('Euler spatial order must be one or two.');
  const reconstruction=spatialOrder===2?linearReconstruction(mesh):undefined;
  const reference = freestream({ mach, alpha, gamma });
  const count = mesh.cells.length;
  if(!['dense','sparse','klu'].includes(linearBackend))throw new Error('Unknown Euler linear backend.');
  const sparse=linearBackend!=='dense';
  if (count > 600&&!sparse) throw new Error('Dense Euler reference is limited to 600 cells.');
  const guess = initial ?? mesh.cells.map(() => reference);
  if (guess.length !== count) throw new Error('Invalid Euler initial-state size.');
  const admissible = x => decodeStates(x, count, reference).every(s => Number.isFinite(s.rho) && s.rho > 0 && s.p > 0);
  let baseJacobian,preconditionerJacobian;
  const jacobian=x=>{
    const states=decodeStates(x,count,reference);
    baseJacobian=eulerJacobian(mesh,states,reference,{reconstruction,sparse});
    if(linearBackend==='sparse'&&reconstruction)preconditionerJacobian=eulerJacobian(mesh,states,reference,{sparse:true});
    return baseJacobian;
  };
  const linearSolve=sparse?(matrix,rhs)=>{
    if(linearBackend==='klu')return solveSparseDirect(matrix,rhs).x;
    let preconditionerMatrix=matrix;
    if(reconstruction){
      // A nearest-neighbor first-order operator supplies a robust approximate
      // inverse for the full reconstructed Jacobian. Include the identical
      // pseudo-time diagonal when it is present in the matrix being solved.
      preconditionerMatrix={...preconditionerJacobian,values:preconditionerJacobian.values.slice()};
      for(let cell=0;cell<count;cell++)for(let row=0;row<4;row++)for(let col=0;col<4;col++){
        const i=4*cell+row,j=4*cell+col,k=sparseIndex(matrix,i,j);
        sparseAdd(preconditionerMatrix,i,j,matrix.values[k]-baseJacobian.values[k]);
      }
    }
    return solveSparse(matrix,rhs,{blockSize:4,preconditionerMatrix}).x;
  }:undefined;
  const solve = (pseudoTransient?solvePseudoTransient:solveNewton)({ initial: encodeStates(guess, reference), admissible, tolerance, maxIterations,
    reference,decode:x=>decodeStates(x,count,reference),
    linearSolve,
    residual: x => eulerResidual(mesh, decodeStates(x, count, reference), reference,{reconstruction}).residual,
    jacobian,onIteration });
  const states = decodeStates(solve.x, count, reference);
  return {...summarizeEuler(mesh, states, reference, solve,{reconstruction}),spatialOrder,linearBackend,pseudoTransient};
}

export function summarizeEuler(mesh, states, reference, solve,{reconstruction}={}) {
  const gamma = reference.gamma;
  const evaluated = eulerResidual(mesh, states, reference,{reconstruction});
  const upstream = totalConditions(reference);
  const entropy = states.map(s => totalConditions(s, gamma).entropy - upstream.entropy);
  const enthalpy = states.map(s => totalConditions(s, gamma).h0 / upstream.h0 - 1);
  return { model: 'euler-fv-reference', converged: solve.converged, reason: solve.reason, history: solve.history,
    accuracyStatus: 'verification-only', cl: null, cd: null, cm: null,
    reference, states, mesh, ...evaluated, diagnostics: { ...evaluated.diagnostics,
      minDensity: Math.min(...states.map(s => s.rho)), minPressure: Math.min(...states.map(s => s.p)),
      maxMach: Math.max(...states.map(s => Math.hypot(s.u, s.v) / Math.sqrt(gamma * s.p / s.rho))),
      maxEntropyError: Math.max(...entropy.map(Math.abs)), maxTotalEnthalpyError: Math.max(...enthalpy.map(Math.abs)) } };
}
