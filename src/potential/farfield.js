// SPDX-License-Identifier: GPL-2.0-or-later
import {affine,sumAffine,addAffine} from './reconstruction.js';

// MSES manual section 1.2.1: Prandtl–Glauert vortex, source, two doublets,
// and the second-order circulation correction. Gradients are with respect
// to physical coordinates, including rotation and the PG chain rule.
export function multipoleBasis(point,{center={x:0,y:0},alpha=0,mach=0,gamma=1.4}={}){
  if(![point.x,point.y,center.x,center.y,alpha,mach,gamma].every(Number.isFinite)||mach<0||mach>=1||gamma<=1)throw new Error('Invalid farfield expansion.');
  const angle=alpha*Math.PI/180,c=Math.cos(angle),s=Math.sin(angle),beta=Math.sqrt(1-mach*mach);
  const x=((point.x-center.x)*c+(point.y-center.y)*s)/beta,y=-(point.x-center.x)*s+(point.y-center.y)*c,r2=x*x+y*y;
  if(!(r2>0))throw new Error('Farfield expansion evaluated at its singularity.');
  const r4=r2*r2,r6=r4*r2,logr=.5*Math.log(r2),theta=Math.atan2(y,x),twoPi=2*Math.PI;
  const a=.25*((3-gamma)/beta+(gamma+1)/beta**3),b=(gamma+1)/16*(1/beta-1/beta**3),scale=(mach/twoPi)**2;
  const cubic=x*x*x-3*x*y*y;
  const potential=[theta/twoPi,logr/twoPi,x/(twoPi*r2),y/(twoPi*r2),scale*(a*logr*x/r2+b*cubic/r4)];
  const gradients=[[-y/(twoPi*r2),x/(twoPi*r2)],
    [x/(twoPi*r2),y/(twoPi*r2)],[(y*y-x*x)/(twoPi*r4),-2*x*y/(twoPi*r4)],
    [-2*x*y/(twoPi*r4),(x*x-y*y)/(twoPi*r4)],
    [scale*(a*(x*x+logr*(y*y-x*x))/r4+b*((3*x*x-3*y*y)/r4-4*x*cubic/r6)),
      scale*(a*x*y*(1-2*logr)/r4+b*(-6*x*y/r4-4*y*cubic/r6))]];
  return{potential,velocity:gradients.map(([gx,gy])=>[gx*c/beta-gy*s,gx*s/beta+gy*c]),theta};
}

// The three free farfield strengths are unknowns in the same Newton system.
// Dirichlet data uses the expansion; projected normal-velocity matching
// determines its source and doublets. This is the fixed-grid counterpart of
// MSES's least-squares directional matching on its moving outer streamlines.
// It does not prescribe strengths from a previous or an inviscid solution.
export function createMultipoleFarfield(mesh,{circulation,offset,center,alpha=0,mach=0,gamma=1.4}={}){
  const gammaColumns=circulation.centers.map((_,i)=>mesh.cells.length+i),columns=Array.from({length:4},(_,i)=>offset+i);
  if(!center){
    const bodies=mesh.contours.map(points=>{
      const te=points[0],le=points.reduce((a,p)=>Math.hypot(p.x-te.x,p.y-te.y)>Math.hypot(a.x-te.x,a.y-te.y)?p:a,te);
      return{te,le,chord:Math.hypot(le.x-te.x,le.y-te.y)};
    });
    const body=bodies.reduce((a,b)=>b.chord>a.chord?b:a);
    center={x:.75*body.le.x+.25*body.te.x,y:.75*body.le.y+.25*body.te.y};
  }
  const settings={center,alpha,mach,gamma},angle=alpha*Math.PI/180,c=Math.cos(angle),s=Math.sin(angle);
  const potential=point=>{
    const basis=multipoleBasis(point,settings),value=affine();
    circulation.centers.forEach((p,e)=>{
      const dx=point.x-p.x,dy=point.y-p.y,physicalTheta=Math.atan2(-dx*s+dy*c,dx*c+dy*s);
      // The difference has no circulation and is single-valued outside all
      // bodies. Wrap the difference, rather than subtracting two branch cuts.
      value.coefficients.set(gammaColumns[e],Math.atan2(Math.sin(basis.theta-physicalTheta),Math.cos(basis.theta-physicalTheta))/(2*Math.PI));
    });
    columns.forEach((col,k)=>value.coefficients.set(col,basis.potential[k+1]));
    return value;
  };
  const velocity=point=>{
    const basis=multipoleBasis(point,settings),v=[affine(c),affine(s)];
    for(let k=0;k<2;k++){
      gammaColumns.forEach(col=>v[k].coefficients.set(col,basis.velocity[0][k]));
      columns.forEach((col,j)=>v[k].coefficients.set(col,basis.velocity[j+1][k]));
    }
    return v;
  };
  const farFaces=mesh.faces.map((f,i)=>f.boundary?.type==='farfield'?i:-1).filter(i=>i>=0);
  if(farFaces.length<3)throw new Error('Multipole matching needs an enclosing farfield.');
  const constraints=({faces})=>{
    const rows=Array.from({length:3},()=>affine()),scales=[0,0,0];
    for(const i of farFaces){
      const f=mesh.faces[i],target=velocity(f),basis=multipoleBasis(f,settings);
      const defect=sumAffine([[faces[i][0],f.nx],[faces[i][1],f.ny],[target[0],-f.nx],[target[1],-f.ny]]);
      for(let k=0;k<3;k++){
        const weight=f.length*(basis.velocity[k+1][0]*f.nx+basis.velocity[k+1][1]*f.ny);
        addAffine(rows[k],defect,weight);scales[k]+=Math.abs(weight);
      }
    }
    const matching=rows.map((row,k)=>{
      if(!(scales[k]>0))throw new Error('Degenerate farfield matching mode.');
      return sumAffine([[row,1/scales[k]]]);
    });
    matching.push({columns:new Set([columns[3],...gammaColumns]),evaluate(x){
      const total=gammaColumns.reduce((sum,col)=>sum+x[col],0);
      return{value:x[columns[3]]-total*total,derivatives:new Map([[columns[3],1],...gammaColumns.map(col=>[col,-2*total])])};
    }});
    return matching;
  };
  return{size:4,columns,center,potential,velocity,constraints};
}
