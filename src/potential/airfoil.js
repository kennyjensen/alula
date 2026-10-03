// SPDX-License-Identifier: GPL-2.0-or-later
import { affine,sumAffine } from './reconstruction.js';
import { createPotentialSystem } from './system.js';
import { createMultipoleFarfield } from './farfield.js';

// Circulation is carried by analytic vortices strictly inside the solid
// bodies; the cell potential remains single-valued throughout the fluid.
// Each additional circulation is closed by its body's trailing-edge Kutta
// condition. No arbitrary wake cut becomes a wall or a zero-crossflow line.
export function airfoilCirculation(mesh,{alpha=0}={}){
  if(!mesh.contours?.length)throw new Error('Airfoil potential needs a body mesh with contours.');
  if(!Number.isFinite(alpha))throw new Error('Invalid incidence.');
  const count=mesh.cells.length,angle=alpha*Math.PI/180;
  const centers=mesh.contours.map(points=>{
    const xmin=Math.min(...points.map(p=>p.x)),xmax=Math.max(...points.map(p=>p.x)),x=.5*(xmin+xmax),ys=[];
    for(let i=1;i<points.length;i++){
      const a=points[i-1],b=points[i];
      if((a.x<x&&b.x>=x)||(b.x<x&&a.x>=x))ys.push(a.y+(b.y-a.y)*(x-a.x)/(b.x-a.x));
    }
    if(ys.length!==2)throw new Error('Circulation center needs an x-monotone body.');
    return{x,y:.5*(ys[0]+ys[1])};
  });
  const baseVelocity=point=>{
    const forms=[affine(Math.cos(angle)),affine(Math.sin(angle))];
    centers.forEach((c,i)=>{
      const dx=point.x-c.x,dy=point.y-c.y,r2=dx*dx+dy*dy;
      if(!(r2>0))throw new Error('Circulation singularity in the fluid.');
      forms[0].coefficients.set(count+i,-dy/(2*Math.PI*r2));
      forms[1].coefficients.set(count+i,dx/(2*Math.PI*r2));
    });
    return forms;
  };
  // The analytic circulation carrier is divergence-free in every fluid
  // cell. Its exact straight-face integral preserves that identity even
  // when the local panel length is comparable to distance from the vortex.
  const baseFluxIntegral=face=>{
    const a=mesh.vertices[face.a],b=mesh.vertices[face.b];
    const flux=affine((Math.cos(angle)*face.nx+Math.sin(angle)*face.ny)*face.length);
    centers.forEach((c,i)=>flux.coefficients.set(count+i,
      -Math.log(Math.hypot(b.x-c.x,b.y-c.y)/Math.hypot(a.x-c.x,a.y-c.y))/(2*Math.PI)));
    return flux;
  };
  const trailingFaces=mesh.contours.map((points,element)=>{
    const te=points[0],chord=Math.hypot(te.x-centers[element].x,te.y-centers[element].y)*2;
    const indices=[];
    mesh.faces.forEach((f,i)=>{if(f.boundary?.element===element&&[mesh.vertices[f.a],mesh.vertices[f.b]].some(v=>Math.hypot(v.x-te.x,v.y-te.y)<chord*1e-10))indices.push(i);});
    if(indices.length!==2)throw new Error('Kutta condition requires exactly two sharp trailing-edge wall faces.');
    return indices;
  });
  return{centers,trailingFaces,baseVelocity,baseFluxIntegral,alpha};
}

export function createAirfoilPotential(mesh,{alpha=0,farfield='fixed',...options}={}){
  if(!['fixed','multipole'].includes(farfield))throw new Error('Unknown farfield condition.');
  const circulation=airfoilCirculation(mesh,{alpha}),{baseVelocity,baseFluxIntegral,centers,trailingFaces}=circulation;
  const matching=farfield==='multipole'?createMultipoleFarfield(mesh,{...options,alpha,circulation,offset:mesh.cells.length+centers.length}):null;
  const system=createPotentialSystem(mesh,{...options,baseVelocity,baseFluxIntegral,
    boundaryPotential:matching?matching.potential:options.boundaryPotential,
    constraints:context=>[...trailingFaces.map((indices,e)=>sumAffine(indices.flatMap(i=>{
    // Use the same two limiting TE station velocities as the coupled BL
    // assembly. Refining outer faces must not move the Kutta sampling point.
    const f=mesh.faces[i],velocity=context.reconstruction.velocity(f.owner,mesh.contours[e][0]);
    return[[velocity[0],-f.ny],[velocity[1],f.nx]];
  }))),...(matching?.constraints(context)??[])]});
  return{...system,centers,trailingFaces,alpha,farfield:matching};
}

export function potentialForces(result,{referenceChord=1,momentOrigin={x:.25,y:0},alpha=result.system.alpha??0}={}){
  const elements=result.system.mesh.contours.map(()=>({cx:0,cy:0,cm:0}));
  result.mesh.faces.forEach((f,i)=>{
    if(f.boundary?.element===undefined)return;
    const force=elements[f.boundary.element],pressure=result.pressureIntegrals[i];
    const fx=pressure*f.nx/referenceChord,fy=pressure*f.ny/referenceChord;
    force.cx+=fx;force.cy+=fy;
    force.cm+=(result.pressureMomentIntegrals[i]+pressure*(momentOrigin.x*f.ny-momentOrigin.y*f.nx))/(referenceChord*referenceChord);
  });
  const angle=alpha*Math.PI/180;
  for(const e of elements){e.cl=e.cy*Math.cos(angle)-e.cx*Math.sin(angle);e.cd=e.cx*Math.cos(angle)+e.cy*Math.sin(angle);}
  return{elements,cl:elements.reduce((s,e)=>s+e.cl,0),cd:elements.reduce((s,e)=>s+e.cd,0),cm:elements.reduce((s,e)=>s+e.cm,0)};
}
