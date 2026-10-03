// SPDX-License-Identifier: GPL-2.0-or-later
import { integratedFaceFlux } from './residual.js';
import { totalConditions } from './gas.js';
import { sparseMatrix,sparseAdd } from '../numerics/sparse.js';

// Assemble derivatives face by face, including every reconstruction neighbor.
// This retains every cross-cell block while avoiding a
// whole-domain evaluation for each perturbed unknown. Geometry derivatives
// in the moving-grid reference still use the independent full residual.
export function eulerJacobian(mesh, states, reference,{reconstruction,sparse=false}={}) {
  const n = 4 * states.length;
  const dependencies=mesh.faces.map((face,index)=>{
    const stencil=reconstruction?.faces[index];
    return reconstruction?[...new Set([stencil.left,...(stencil.right?[stencil.right]:[])].flatMap(s=>[s.cell,...s.neighbors.map(a=>a.cell)]))]
      :face.neighbor===null?[face.owner]:[face.owner,face.neighbor];
  });
  let matrix;
  if(sparse){
    const columns=Array.from({length:n},()=>new Set());
    mesh.faces.forEach((face,index)=>{
      for(const rowCell of [face.owner,...(face.neighbor===null?[]:[face.neighbor])])for(let row=0;row<4;row++)
        for(const colCell of dependencies[index])for(let col=0;col<4;col++)columns[4*rowCell+row].add(4*colCell+col);
    });
    matrix=sparseMatrix(columns);
  }else matrix=new Float64Array(n*n);
  const add=sparse?(row,col,value)=>sparseAdd(matrix,row,col,value):(row,col,value)=>{matrix[row*n+col]+=value;};
  const speed = Math.hypot(reference.u, reference.v); const mass = reference.rho * speed;
  const q2 = mass * speed; const scales = [mass, q2, q2, mass * totalConditions(reference).h0];
  for (const [index,face] of mesh.faces.entries()) {
    const stencil=reconstruction?.faces[index];
    const left = reconstruction?reconstruction.sample(stencil.left,states):states[face.owner];
    const right = reconstruction?reconstruction.sample(stencil.right,states):face.neighbor === null ? null : states[face.neighbor];
    for (const cell of dependencies[index]) {
      const state = states[cell];
      const values = [Math.log(state.rho / reference.rho), state.u, state.v, (state.p - reference.p) / q2];
      for (let variable = 0; variable < 4; variable++) {
        const h = Math.cbrt(Number.EPSILON) * Math.max(1, Math.abs(values[variable]));
        const shifted = sign => {
          const s = { ...state };
          if (variable === 0) s.rho *= Math.exp(sign * h);
          else if (variable === 1) s.u += sign * h;
          else if (variable === 2) s.v += sign * h;
          else s.p += sign * h * q2;
          return s;
        };
        const plus = shifted(1); const minus = shifted(-1);
        const flux=shift=>integratedFaceFlux(reconstruction?reconstruction.sample(stencil.left,states,cell,shift):cell===face.owner?shift:left,
          reconstruction?reconstruction.sample(stencil.right,states,cell,shift):cell===face.neighbor?shift:right,face,reference);
        let fp,fm;
        try{fp=flux(plus);}catch{}
        try{fm=flux(minus);}catch{}
        if(!fp&&!fm)throw new Error('No admissible reconstructed Euler derivative.');
        const width=fp&&fm?2*h:h;
        fp??=integratedFaceFlux(left,right,face,reference);
        fm??=integratedFaceFlux(left,right,face,reference);
        for (let equation = 0; equation < 4; equation++) {
          const d = (fp[equation] - fm[equation]) / (width * scales[equation]);
          add(4*face.owner+equation,4*cell+variable,d/mesh.cells[face.owner].perimeter);
          if(face.neighbor!==null)add(4*face.neighbor+equation,4*cell+variable,-d/mesh.cells[face.neighbor].perimeter);
        }
      }
    }
  }
  return matrix;
}
