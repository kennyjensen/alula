// SPDX-License-Identifier: GPL-2.0-or-later
// Exact integration of linearly interpolated nodal Cp on a CCW polygon.
// Coordinates and thicknesses are normalized by the common reference chord.
export function pressureForces(points,cp,{alpha=0,momentOrigin={x:.25,y:0}}={}){
  if(points.length!==cp.length||points.length<4||!cp.every(Number.isFinite))throw new Error('Pressure integration needs one finite Cp per contour vertex.');
  let cx=0,cy=0,cm=0;
  for(let i=1;i<points.length;i++){
    const a=points[i-1],b=points[i],dx=b.x-a.x,dy=b.y-a.y,average=.5*(cp[i-1]+cp[i]),difference=cp[i]-cp[i-1];
    cx-=average*dy;cy+=average*dx;
    cm-=average*(dx*(.5*(a.x+b.x)-momentOrigin.x)+dy*(.5*(a.y+b.y)-momentOrigin.y))+difference*(dx*dx+dy*dy)/12;
  }
  const angle=alpha*Math.PI/180;
  return{cx,cy,cl:cy*Math.cos(angle)-cx*Math.sin(angle),cd:cx*Math.cos(angle)+cy*Math.sin(angle),cm};
}
