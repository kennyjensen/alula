// SPDX-License-Identifier: GPL-2.0-or-later
import { affine,addAffine,sumAffine,evaluateAffine,packAffine,evaluatePackedAffine,potentialReconstruction } from './reconstruction.js';
import { isentropicState,isentropicMassFlux } from './isentropic.js';
import { sparseMatrix,sparseAdd } from '../numerics/sparse.js';
import { solvePotentialLinear } from './linear.js';
import { normInf,solveLinear } from '../numerics/linear.js';
import { solveNewton } from '../numerics/newton.js';

// One scalar potential per cell. Extra affine unknowns (e.g. circulation)
// must each have a closing constraint. A shared face has a single mass flux,
// including the density derivative; the two cells receive opposite fluxes.
export function createPotentialSystem(mesh,{mach=0,gamma=1.4,order=2,
  baseVelocity=()=>[affine(1),affine()],baseFluxIntegral,boundaryPotential=()=>affine(),fluxQuadrature=3,
  boundaryNormalVelocity=()=>affine(),constraints=()=>[],sparse=false}={}){
  const rules={1:[[0,1]],2:[[-1/Math.sqrt(3),.5],[1/Math.sqrt(3),.5]],
    3:[[-Math.sqrt(.6),5/18],[0,4/9],[Math.sqrt(.6),5/18]]};
  if(!rules[fluxQuadrature])throw new Error('Use one, two or three face quadrature points.');
  const conditions={mach,gamma};isentropicState(1,0,conditions);
  const reconstruction=potentialReconstruction(mesh,{baseVelocity,boundaryPotential,boundaryNormalVelocity,order});
  const gradients=mesh.cells.map((c,i)=>reconstruction.gradient(i,c));
  const faceCorrection=mesh.faces.map(f=>{
    const left=reconstruction.gradient(f.owner,f),center=mesh.cells[f.owner];
    let delta,dx,dy;
    if(f.neighbor!==null){
      const other=mesh.cells[f.neighbor];
      dx=other.x-center.x;dy=other.y-center.y;
      // Trapezoidal line integral of the two center gradients is exact for a
      // quadratic potential. The defect enforces direct neighbor coupling.
      delta=affine(0,[[f.neighbor,1],[f.owner,-1]]);
      for(let k=0;k<2;k++)for(const cell of [f.owner,f.neighbor])addAffine(delta,gradients[cell][k],-.5*(k===0?dx:dy));
    }else if(f.boundary.type!=='wall'){
      dx=f.x-center.x;dy=f.y-center.y;
      delta=sumAffine([[boundaryPotential(f),1],[affine(0,[[f.owner,1]]),-1]]);
      for(let k=0;k<2;k++)for(const g of [gradients[f.owner][k],left[k]])addAffine(delta,g,-.5*(k===0?dx:dy));
    }else return null;
    const dn=dx*f.nx+dy*f.ny;
    if(!(dn>0))throw new Error('Potential flux requires positive cell-center normal spacing.');
    return sumAffine([[delta,1/dn]]);
  });
  const faceVelocity=(f,i,point)=>{
    let gradient=reconstruction.gradient(f.owner,point);
    if(f.neighbor!==null){const right=reconstruction.gradient(f.neighbor,point);gradient=gradient.map((g,k)=>sumAffine([[g,.5],[right[k],.5]]));}
    const base=baseVelocity(point),velocity=gradient.map((g,k)=>addAffine(g,base[k]));
    if(f.boundary?.type==='wall'){
      // A wall-normal unknown is constant over its face. Reconstruct the
      // varying tangential speed for the density at every integration point.
      const defect=sumAffine([[boundaryNormalVelocity(f),1],[velocity[0],-f.nx],[velocity[1],-f.ny]]);
      return velocity.map((v,k)=>addAffine(v,defect,k===0?f.nx:f.ny));
    }
    // The center-to-center (or center-to-face) defect is a constant normal
    // correction. A shared wake potential is prescribed at its midpoint,
    // not separately at each quadrature location.
    return velocity.map((v,k)=>addAffine(v,faceCorrection[i],k===0?f.nx:f.ny));
  };
  const faces=mesh.faces.map((f,i)=>faceVelocity(f,i,f));
  const quadrature=mesh.faces.map((f,i)=>rules[fluxQuadrature].map(([s,weight])=>{
    const a=mesh.vertices[f.a],b=mesh.vertices[f.b],point={x:f.x+s*(b.x-a.x)/2,y:f.y+s*(b.y-a.y)/2};
    return{weight,point,velocity:s===0?faces[i]:faceVelocity(f,i,point)};
  }));
  const meanVelocities=quadrature.map(qs=>[0,1].map(k=>sumAffine(qs.map(q=>[q.velocity[k],q.weight]))));
  // Replace only the quadrature error in the incompressible analytic carrier.
  // The density-dependent part still uses the actual local velocity. This
  // gives exact carrier mass conservation at Mach zero without changing the
  // nonlinear isentropic closure or its derivatives.
  const carrierCorrections=mesh.faces.map((f,i)=>{
    if(!baseFluxIntegral||f.boundary?.type==='wall')return affine();
    const correction=sumAffine([[baseFluxIntegral(f),1]]);
    for(const q of quadrature[i]){const base=baseVelocity(q.point);addAffine(correction,base[0],-q.weight*f.length*f.nx);addAffine(correction,base[1],-q.weight*f.length*f.ny);}
    return correction;
  });
  const fluxColumns=quadrature.map((qs,i)=>new Set([...carrierCorrections[i].coefficients.keys(),...qs.flatMap(q=>q.velocity.flatMap(f=>[...f.coefficients.keys()]))]));
  const extra=constraints({faces,reconstruction,fluxColumns,meanVelocities}),count=mesh.cells.length,n=count+extra.length;
  const cellVelocities=mesh.cells.map((c,i)=>reconstruction.velocity(i,c));
  const forms=[...quadrature.flatMap(qs=>qs.flatMap(q=>q.velocity)),...carrierCorrections,...extra.filter(f=>!f.evaluate)];
  if(forms.some(f=>[...f.coefficients.keys()].some(i=>!Number.isInteger(i)||i<0||i>=n)))throw new Error('Unclosed potential unknown.');
  const packed=new WeakMap(),pack=form=>{if(!packed.has(form))packed.set(form,packAffine(form));return packed.get(form);};
  const packedFaces=faces.map(pair=>pair.map(pack)),packedCells=cellVelocities.map(pair=>pair.map(pack));
  const admissibleForms=[...packedFaces,...packedCells];
  for(let i=0;i<quadrature.length;i++)for(const q of quadrature[i]){
    q.midpoint=q.velocity===faces[i];q.packed=q.midpoint?packedFaces[i]:q.velocity.map(pack);
    if(!q.midpoint)admissibleForms.push(q.packed);
    // All derivative columns and mean-velocity constraints are already
    // assembled. Non-midpoint Maps can now be collected; evaluation uses
    // their immutable packed coefficients in exactly the same order.
    delete q.velocity;
  }
  const packedCarriers=carrierCorrections.map(pack);
  const pattern=Array.from({length:n},()=>new Set());
  mesh.faces.forEach((f,i)=>{for(const row of [f.owner,...(f.neighbor===null?[]:[f.neighbor])])for(const col of fluxColumns[i])pattern[row].add(col);});
  extra.forEach((form,i)=>{for(const col of form.evaluate?form.columns:form.coefficients.keys()){
    if(!Number.isInteger(col)||col<0||col>=n)throw new Error('Unclosed potential constraint unknown.');
    pattern[count+i].add(col);
  }});
  const matrix=sparse?sparseMatrix(pattern):null;
  const evaluate=(x,{jacobian=false}={},gas=conditions)=>{
    if(x.length!==n||!Array.from(x).every(Number.isFinite))throw new Error('Invalid potential state vector.');
    const residual=new Float64Array(n),j=jacobian?(sparse?{...matrix,values:new Float64Array(matrix.values.length)}:new Float64Array(n*n)):null;
    const states=[],fluxes=new Float64Array(faces.length),fluxDerivatives=[];
    const pressureIntegrals=new Float64Array(faces.length),pressureMomentIntegrals=new Float64Array(faces.length);
    let boundaryMass=0,boundaryThroughput=0,wallLeakage=0,interfaceMass=0,maxMachSquared=0;
    mesh.faces.forEach((f,i)=>{
      const u=evaluatePackedAffine(packedFaces[i][0],x),v=evaluatePackedAffine(packedFaces[i][1],x),state=isentropicMassFlux(u,v,f.nx,f.ny,gas);
      states.push({u,v,...state});maxMachSquared=Math.max(maxMachSquared,state.machSquared);
      let flux=evaluatePackedAffine(packedCarriers[i],x);
      const derivatives=jacobian?new Map(carrierCorrections[i].coefficients):null;
      for(const q of quadrature[i]){
        const uq=q.midpoint?u:evaluatePackedAffine(q.packed[0],x),vq=q.midpoint?v:evaluatePackedAffine(q.packed[1],x);
        const sq=q.midpoint?state:isentropicMassFlux(uq,vq,f.nx,f.ny,gas),weight=q.weight*f.length;
        maxMachSquared=Math.max(maxMachSquared,sq.machSquared);flux+=weight*sq.flux;
        pressureIntegrals[i]+=weight*sq.cp;
        pressureMomentIntegrals[i]+=weight*sq.cp*(q.point.y*f.nx-q.point.x*f.ny);
        if(j)for(let k=0;k<2;k++){
          const form=q.packed[k];
          for(let t=0;t<form.columns.length;t++){const col=form.columns[t];derivatives.set(col,(derivatives.get(col)??0)+weight*sq.derivative[k]*form.weights[t]);}
        }
      }
      fluxes[i]=flux;fluxDerivatives.push(derivatives);
      for(const [row,sign] of [[f.owner,1],...(f.neighbor===null?[]:[[f.neighbor,-1]])]){
        const scale=sign/mesh.cells[row].perimeter;
        residual[row]+=scale*flux;
        if(j)for(const [col,value] of derivatives){
          const entry=scale*value;
          if(sparse)sparseAdd(j,row,col,entry);else j[row*n+col]+=entry;
        }
      }
      if(f.neighbor===null){boundaryMass+=flux;boundaryThroughput+=Math.abs(flux);if(f.boundary.type==='wall')wallLeakage=Math.max(wallLeakage,Math.abs(flux/f.length));if(f.boundary.type==='potential-interface')interfaceMass+=flux;}
    });
    extra.forEach((form,i)=>{
      const value=form.evaluate?form.evaluate(x,{states,fluxes,fluxDerivatives,jacobian}):{value:evaluateAffine(form,x),derivatives:form.coefficients};
      residual[count+i]=value.value;
      if(j)for(const [col,v] of value.derivatives){if(sparse)sparseAdd(j,count+i,col,v);else j[(count+i)*n+col]=v;}
    });
    const cellStates=packedCells.map(forms=>{
      const u=evaluatePackedAffine(forms[0],x),v=evaluatePackedAffine(forms[1],x),state=isentropicState(u,v,gas);
      maxMachSquared=Math.max(maxMachSquared,state.machSquared);return{u,v,...state};
    });
    let volumeMass=0;for(let i=0;i<count;i++)volumeMass+=residual[i]*mesh.cells[i].perimeter;
    return{residual,jacobian:j,faces:states,cellStates,fluxes,pressureIntegrals,pressureMomentIntegrals,diagnostics:{residual:normInf(residual),massResidual:normInf(residual.subarray(0,count)),
      constraintResidual:normInf(residual.subarray(count)),boundaryMass,relativeMassImbalance:Math.abs(boundaryMass)/Math.max(boundaryThroughput,1e-30),
      sharedFluxCancellation:Math.abs(volumeMass-boundaryMass),interfaceMass,externalMass:boundaryMass-interfaceMass,wallLeakage,maxMach:Math.sqrt(maxMachSquared),linearFallbacks:reconstruction.linearFallbacks,expandedStencils:reconstruction.expandedStencils,fluxQuadrature}};
  };
  // Primitive admissibility must not evaluate discrete BL/transition
  // equations: their active set is updated after a trial passes this check.
  const admissible=(x,gas)=>{
    if(x.length!==n||!Array.from(x).every(Number.isFinite))return false;
    try{for(const forms of admissibleForms)isentropicState(evaluatePackedAffine(forms[0],x),evaluatePackedAffine(forms[1],x),gas);return true;}catch{return false;}
  };
  const view=gas=>{
    let lowOrder;
    const preconditioner=x=>{
      if(order===1)return evaluate(x,{jacobian:true},gas).jacobian;
      lowOrder??=createPotentialSystem(mesh,{...gas,order:1,baseVelocity,baseFluxIntegral,boundaryPotential,boundaryNormalVelocity,constraints,sparse,fluxQuadrature});
      return lowOrder.evaluate(x,{jacobian:true}).jacobian;
    };
    return{mesh,n,count,conditions:gas,reconstruction,faces,cellVelocities,constraints:extra,sparse,preconditioner,
      evaluate:(x,options)=>evaluate(x,options,gas),admissible:x=>admissible(x,gas),
      // Keep the prescribed boundary/constraint callbacks and share every
      // geometric operator. This is an initialization of the same problem,
      // not a new physical Mach-zero case (whose boundary data may differ).
      incompressibleInitialization:()=>view({mach:0,gamma})};
  };
  return view(conditions);
}

export function solvePotential(system,{initial=new Float64Array(system.n),tolerance=1e-10,maxIterations=40,onIteration,linearSolve,
  linearTolerance=1e-7,fillLevel=1,linearBackend='klu'}={}){
  let preconditionerMatrix;
  const result=solveNewton({initial,tolerance,maxIterations,onIteration,
    residual:x=>system.evaluate(x).residual,jacobian:x=>{
      if(system.sparse&&!linearSolve&&linearBackend==='gmres')preconditionerMatrix=system.preconditioner(x);
      return system.evaluate(x,{jacobian:true}).jacobian;
    },
    admissible:x=>{try{system.evaluate(x);return true;}catch{return false;}},
    // Inexact Newton corrections are refined by subsequent nonlinear steps;
    // the unmodified full residual still decides final convergence.
    linearSolve:linearSolve??(system.sparse?(a,b)=>solvePotentialLinear(a,b,{preconditionerMatrix,linearTolerance,fillLevel,linearBackend}):solveLinear)});
  const evaluated=system.evaluate(result.x);
  return{...result,...evaluated,mesh:system.mesh,states:evaluated.cellStates,system};
}
