// SPDX-License-Identifier: GPL-2.0-or-later
import {sumAffine} from './reconstruction.js';

// Tangential velocity over any fraction of an original BL panel. Integrate
// each intersected outer-mesh face separately, including when the requested
// half-panel boundary cuts through a face. This keeps the sampling interval
// independent of the outer triangulation and its subdivision parity.
export function integratedPanelVelocity(panel,segments,velocityAt,{lo=0,hi=1}={}){
  if(!(0<=lo&&lo<hi&&hi<=1))throw new Error('Invalid panel sampling interval.');
  const terms=[],project=p=>((p.x-panel.a.x)*panel.tx+(p.y-panel.a.y)*panel.ty)/panel.length;
  const rule=[[-Math.sqrt(.6),5/18],[0,4/9],[Math.sqrt(.6),5/18]];let covered=0;
  for(const face of segments){
    const a=project(face.a),b=project(face.b),left=Math.max(lo,Math.min(a,b)),right=Math.min(hi,Math.max(a,b));
    if(right<=left)continue;
    const fraction=(right-left)/(hi-lo);covered+=fraction;
    for(const [z,w]of rule){
      const t=.5*(left+right+z*(right-left)),point={x:panel.a.x+t*panel.length*panel.tx,y:panel.a.y+t*panel.length*panel.ty};
      const [u,v]=velocityAt(face,point);terms.push([u,fraction*w*panel.tx],[v,fraction*w*panel.ty]);
    }
  }
  if(Math.abs(covered-1)>1e-8)throw new Error('Incomplete BL panel velocity coverage.');
  return sumAffine(terms);
}

// Recover nodal signed speed from interval-averaged tangential velocities.
// A linear field is reproduced exactly on nonuniform intervals, including
// one-sided TE extrapolation. Smooth-field error is second order. Averaging
// over a geometric BL panel avoids sampling a faceted-wall corner directly.
export function interpolatePanelAverages(samples,positions){
  if(samples.length<2||samples.some((p,i)=>!Number.isFinite(p.s)||(i&&p.s<=samples[i-1].s)))throw new Error('Panel velocity samples must have distinct ordered arc coordinates.');
  let left=0;
  return positions.map(s=>{
    if(!Number.isFinite(s))throw new Error('Invalid BL arc coordinate.');
    while(left<samples.length-2&&samples[left+1].s<s)left++;
    while(left>0&&samples[left].s>s)left--;
    const a=samples[left],b=samples[left+1],weight=(s-a.s)/(b.s-a.s);
    return sumAffine([[a.velocity,1-weight],[b.velocity,weight]]);
  });
}
