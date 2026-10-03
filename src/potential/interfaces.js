// SPDX-License-Identifier: GPL-2.0-or-later
import { affine,evaluateAffine } from './reconstruction.js';

// A thin source wake has continuous potential and different normal
// derivatives on its two sides. Split only its reconstruction topology:
// a common face potential and a mass-jump equation reconnect both sides.
export function partitionPotentialInterfaces(mesh,indices){
  const faces=mesh.faces.map(f=>({...f})),cells=mesh.cells.map(c=>({...c,faces:c.faces.map(f=>({...f}))})),interfaces=[];
  if(new Set(indices).size!==indices.length)throw new Error('Repeated potential interface.');
  for(const index of indices){
    const original=mesh.faces[index];
    if(!original||original.neighbor===null)throw new Error('Potential interface must be an internal fluid face.');
    const id=interfaces.length,rightIndex=faces.length,tag={type:'potential-interface',interface:id};
    faces[index]={...original,neighbor:null,boundary:{...tag,side:0}};
    faces.push({...original,a:original.b,b:original.a,owner:original.neighbor,neighbor:null,nx:-original.nx,ny:-original.ny,boundary:{...tag,side:1}});
    const ref=cells[original.neighbor].faces.find(r=>r.face===index);ref.face=rightIndex;ref.sign=1;
    interfaces.push({faces:[index,rightIndex],original:index});
  }
  return{...mesh,cells,faces,potentialInterfaces:interfaces};
}

// Sum outward one-sided mass fluxes plus a specified sheet source. The
// source may itself depend on external unknowns through an affine form.
// A nonlinear BL mass source will use the same value/derivatives contract.
export function potentialFluxConstraint(mesh,fluxColumns,indices,{source=affine(),scale=1}={}){
  if(!(scale>0))throw new Error('Invalid mass constraint scale.');
  const columns=new Set(source.evaluate?source.columns:source.coefficients.keys());
  for(const i of indices)for(const col of fluxColumns[i])columns.add(col);
  return{columns,evaluate(x,{fluxes,fluxDerivatives,jacobian}){
    const s=source.evaluate?source.evaluate(x):{value:evaluateAffine(source,x),derivatives:source.coefficients};
    let value=s.value;const derivatives=new Map();
    const add=(col,v)=>derivatives.set(col,(derivatives.get(col)??0)+v/scale);
    if(jacobian)for(const [col,v] of s.derivatives)add(col,v);
    for(const i of indices){
      value+=fluxes[i];
      if(jacobian)for(const [col,v] of fluxDerivatives[i])add(col,v);
    }
    return{value:value/scale,derivatives};
  }};
}
