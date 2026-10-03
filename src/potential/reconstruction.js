// SPDX-License-Identifier: GPL-2.0-or-later
import { factorLinear } from '../numerics/linear.js';

// Sparse affine forms let geometry-dependent potential reconstruction be
// assembled once while retaining exact circulation and boundary derivatives.
export const affine=(constant=0,coefficients=[])=>({constant,coefficients:new Map(coefficients)});
export function addAffine(target,source,weight=1){
  target.constant+=weight*source.constant;
  for(const [col,value] of source.coefficients){const next=(target.coefficients.get(col)??0)+weight*value;if(next===0)target.coefficients.delete(col);else target.coefficients.set(col,next);}
  return target;
}
export const evaluateAffine=(form,x)=>{let value=form.constant;for(const [col,weight] of form.coefficients)value+=weight*x[col];return value;};
// Immutable evaluation layout. Preserve Map insertion order so the floating
// point sum is identical, while avoiding an iterator pair for every term.
export const packAffine=form=>({constant:form.constant,columns:Int32Array.from(form.coefficients.keys()),weights:Float64Array.from(form.coefficients.values())});
export const evaluatePackedAffine=(form,x)=>{let value=form.constant;for(let i=0;i<form.columns.length;i++)value+=form.weights[i]*x[form.columns[i]];return value;};
export const sumAffine=(terms)=>{const result=affine();for(const [form,weight] of terms)addAffine(result,form,weight);return result;};

// Quadratic least squares for scalar potential, including boundary equations
// at their actual face locations. Symmetric Hessians preserve irrotationality
// within each reconstructed polynomial. Extend a rank-deficient stencil
// through connected fluid neighbors before falling back to linear order.
// No diagonal regularization changes the polynomial equations.
export function potentialReconstruction(mesh,{baseVelocity=()=>[affine(),affine()],
  boundaryPotential=()=>affine(),boundaryNormalVelocity=()=>affine(),order=2}={}){
  if(![1,2].includes(order))throw new Error('Potential reconstruction order must be one or two.');
  const neighbors=mesh.cells.map(()=>new Set()),boundaryFaces=mesh.cells.map(()=>[]);
  mesh.faces.forEach((f,index)=>{
    if(f.neighbor===null)boundaryFaces[f.owner].push(index);
    else{neighbors[f.owner].add(f.neighbor);neighbors[f.neighbor].add(f.owner);}
  });
  let linearFallbacks=0,expandedStencils=0;
  const jets=mesh.cells.map((cell,index)=>{
    const length=Math.sqrt(cell.area),stencil=new Set(neighbors[index]);
    if(order===2)for(const j of neighbors[index])for(const k of neighbors[j])if(k!==index)stencil.add(k);
    const rows=[];
    const valueRow=(point,target)=>{
      const x=(point.x-cell.x)/length,y=(point.y-cell.y)/length;
      rows.push({a:[x,y,.5*x*x,x*y,.5*y*y],target,weight:1/(x*x+y*y)});
    };
    for(const j of stencil)valueRow(mesh.cells[j],affine(0,[[j,1],[index,-1]]));
    const boundaries=new Set(boundaryFaces[index]);
    if(order===2)for(const j of neighbors[index])for(const f of boundaryFaces[j])boundaries.add(f);
    for(const f of boundaries){
      const face=mesh.faces[f],x=(face.x-cell.x)/length,y=(face.y-cell.y)/length;
      if(face.boundary.type==='wall'){
        const [u,v]=baseVelocity(face),target=sumAffine([[boundaryNormalVelocity(face),length],[u,-length*face.nx],[v,-length*face.ny]]);
        rows.push({a:[face.nx,face.ny,face.nx*x,face.nx*y+face.ny*x,face.ny*y],target,weight:1});
      }else valueRow(face,sumAffine([[boundaryPotential(face),1],[affine(0,[[index,1]]),-1]]));
    }
    const fit=size=>{
      const matrix=new Float64Array(size*size),rhs=Array.from({length:size},()=>affine());
      for(const row of rows)for(let i=0;i<size;i++){
        addAffine(rhs[i],row.target,row.weight*row.a[i]);
        for(let j=0;j<size;j++)matrix[i*size+j]+=row.weight*row.a[i]*row.a[j];
      }
      const solve=factorLinear(matrix,size),columns=new Set(rhs.flatMap(r=>[...r.coefficients.keys()]));
      const constants=solve(rhs.map(r=>r.constant)),result=Array.from(constants,v=>affine(v));
      for(const col of columns){const values=solve(rhs.map(r=>r.coefficients.get(col)??0));for(let i=0;i<size;i++)if(values[i]!==0)result[i].coefficients.set(col,values[i]);}
      return result;
    };
    let result,rings=order===2?2:1;
    if(order===2){
      for(;;){
        try{result=fit(5);break;}
        catch{
          if(rings===4){linearFallbacks++;result=fit(2);break;}
          const added=new Set();
          for(const j of stencil)for(const k of neighbors[j])if(k!==index&&!stencil.has(k))added.add(k);
          if(!added.size){linearFallbacks++;result=fit(2);break;}
          if(rings===2)expandedStencils++;rings++;
          for(const j of added){stencil.add(j);valueRow(mesh.cells[j],affine(0,[[j,1],[index,-1]]));}
        }
      }
    }else result=fit(2);
    return{length,result,rings};
  });
  const gradient=(cellIndex,point=mesh.cells[cellIndex])=>{
    const cell=mesh.cells[cellIndex],{length,result:r}=jets[cellIndex];
    const x=(point.x-cell.x)/length,y=(point.y-cell.y)/length;
    return[sumAffine([[r[0],1/length],...(r.length===5?[[r[2],x/length],[r[3],y/length]]:[])]),
      sumAffine([[r[1],1/length],...(r.length===5?[[r[3],x/length],[r[4],y/length]]:[])])];
  };
  const velocity=(cell,point)=>gradient(cell,point).map((g,k)=>addAffine(g,baseVelocity(point??mesh.cells[cell])[k]));
  return{gradient,velocity,jets,linearFallbacks,expandedStencils};
}
