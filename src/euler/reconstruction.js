// SPDX-License-Identifier: GPL-2.0-or-later
// Unlimited linear reconstruction for smooth subcritical Euler verification.
// Gradients use weighted least squares in physical coordinates. Wall pressure
// is extrapolated from the fluid stencil; imposing dp/dn=0 at a curved wall
// would discard the normal pressure gradient required by streamline curvature.
export function linearReconstruction(mesh){
  const neighbors=mesh.cells.map(()=>new Set());
  for(const f of mesh.faces)if(f.neighbor!==null){neighbors[f.owner].add(f.neighbor);neighbors[f.neighbor].add(f.owner);}
  const gradients=mesh.cells.map((center,i)=>{
    const entries=[...neighbors[i]].map(j=>{const dx=mesh.cells[j].x-center.x,dy=mesh.cells[j].y-center.y;return{cell:j,dx,dy,w:1/(dx*dx+dy*dy)};});
    let xx=0,xy=0,yy=0;
    for(const a of entries){xx+=a.w*a.dx*a.dx;xy+=a.w*a.dx*a.dy;yy+=a.w*a.dy*a.dy;}
    const determinant=xx*yy-xy*xy,trace=xx+yy;
    if(determinant>1e-12*trace*trace)return entries.map(a=>({cell:a.cell,gx:a.w*(yy*a.dx-xy*a.dy)/determinant,gy:a.w*(xx*a.dy-xy*a.dx)/determinant}));
    // A one-cell-high channel resolves only one gradient direction. Use the
    // rank-one pseudoinverse rather than inventing a transverse gradient.
    const length=xx>=yy?Math.hypot(xx,xy):Math.hypot(xy,yy);
    if(length===0)return[];
    const ex=(xx>=yy?xx:xy)/length,ey=(xx>=yy?xy:yy)/length;
    return entries.map(a=>{const value=a.w*(a.dx*ex+a.dy*ey)/trace;return{cell:a.cell,gx:value*ex,gy:value*ey};});
  });
  const stencil=(cell,face)=>{
    if(cell===null)return null;
    const center=mesh.cells[cell],a=mesh.vertices[face.a],b=mesh.vertices[face.b];
    const dx=.5*(a.x+b.x)-center.x,dy=.5*(a.y+b.y)-center.y;
    return{cell,neighbors:gradients[cell].map(g=>({cell:g.cell,weight:g.gx*dx+g.gy*dy}))};
  };
  const faces=mesh.faces.map(f=>({left:stencil(f.owner,f),right:stencil(f.neighbor,f)}));
  const sample=(s,states,overrideCell=-1,overrideState)=>{
    if(!s)return null;
    const get=i=>i===overrideCell?overrideState:states[i],base=get(s.cell),result={...base};
    for(const key of ['rho','u','v','p'])for(const n of s.neighbors)result[key]+=n.weight*(get(n.cell)[key]-base[key]);
    return result;
  };
  return{faces,sample,gradients};
}
