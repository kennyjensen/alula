// SPDX-License-Identifier: GPL-2.0-or-later
import { isentropicState } from '../potential/isentropic.js';
import {nodalMassSlopes,integrateMassDerivative,addMassWeights} from './mass-interpolation.js';

// Density and speed use freestream units; delta* uses the reference chord.
// Surface speed is signed along its CCW contour, while wake speed is positive
// downstream. Density belongs inside the differentiated mass deficit.
export function displacementMass(q,deltaStar,conditions){
  if(!Number.isFinite(deltaStar)||deltaStar<0)throw new Error('Invalid displacement thickness.');
  const gas=isentropicState(q,0,conditions);
  return{value:gas.rho*q*deltaStar,velocityDerivative:deltaStar*(gas.rho+2*q*q*gas.rhoSpeedSquared),thicknessDerivative:gas.rho*q};
}

export function displacementSource(weights,{offset=0,thicknessScale=1,mach=0,gamma=1.4}={}){
  const coefficients=[...weights].filter(([,v])=>v!==0),columns=new Set(coefficients.flatMap(([id])=>[offset+4*id+2,offset+4*id+3]));
  return{columns,weights:new Map(coefficients),evaluate(x){
    let value=0;const derivatives=new Map();
    for(const [id,weight] of coefficients){
      const i=offset+4*id,state=displacementMass(x[i+3],x[i+2]*thicknessScale,{mach,gamma});
      value+=weight*state.value;
      derivatives.set(i+2,weight*state.thicknessDerivative*thicknessScale);
      derivatives.set(i+3,weight*state.velocityDerivative);
    }
    return{value,derivatives};
  }};
}

// Integrate the selected body source and quadratic wake source
// representation over each mesh face. Subdivision cannot change the mass
// budget, and each side of a sharp TE retains its own signed nodal state.
export function meshDisplacementSources(mesh,outer,{bodyMassInterpolation='linear',...options}={}){
  if(!['linear','hermite'].includes(bodyMassInterpolation))throw new Error('Unknown body mass interpolation.');
  const bodySources=new Map(),wakeSources=new Map(),bodyPanels=new Map(),tolerance=1e-9;
  const slopes=bodyMassInterpolation==='hermite'?outer.bodies.map(b=>nodalMassSlopes(b.s,b.start)):null;
  const bodyPanel=f=>{
    const body=outer.bodies[f.boundary.element];
    for(let k=body.first;k<=body.last;k++){
      const p=outer.panels[k],dx=f.x-p.a.x,dy=f.y-p.a.y,t=(dx*p.tx+dy*p.ty)/p.length;
      if(t>=-tolerance&&t<=1+tolerance&&Math.abs(-dx*p.ty+dy*p.tx)<tolerance*p.length)return p;
    }
    throw new Error('Mesh wall is not a subdivision of the BL contour.');
  };
  mesh.faces.forEach((f,i)=>{
    if(f.boundary?.element===undefined)return;
    const p=bodyPanel(f),weight=f.length/p.length;
    bodyPanels.set(i,p);
    let weights=new Map([[p.node,-weight],[p.node+1,weight]]);
    if(slopes){
      const body=outer.bodies[f.boundary.element],j=p.node-body.start,d=slopes[f.boundary.element];
      const projection=v=>((v.x-p.a.x)*p.tx+(v.y-p.a.y)*p.ty)/p.length;
      const a=projection(mesh.vertices[f.a]),b=projection(mesh.vertices[f.b]);
      weights=integrateMassDerivative(p.node,p.length,d[j],d[j+1],Math.min(a,b),Math.max(a,b));
    }
    bodySources.set(i,displacementSource(weights,options));
  });
  for(const cut of mesh.cuts.filter(c=>c.type==='wake-cut')){
    const w=outer.wakes[cut.element],f=mesh.faces[cut.face],weights=new Map();
    if(f.x>w.points.at(-1).x){wakeSources.set(cut.face,displacementSource(weights,options));continue;}
    let k=0;while(k<w.segments.length-1&&w.points[k+1].x<f.x)k++;
    const p=w.segments[k],projection=v=>((v.x-p.a.x)*p.tx+(v.y-p.a.y)*p.ty)/p.length;
    const a=projection(mesh.vertices[f.a]),b=projection(mesh.vertices[f.b]),lo=Math.min(a,b),hi=Math.max(a,b);
    if(lo<-tolerance||hi>1+tolerance||Math.abs((f.x-p.a.x)*p.ty-(f.y-p.a.y)*p.tx)>tolerance*p.length)throw new Error('Mesh wake face is not aligned with the BL wake segment.');
    for(let degree=0;degree<3;degree++){
      const integral=p.length*(hi**(degree+1)-lo**(degree+1))/(degree+1),row=3*(w.sourceStart+k)+degree;
      for(let col=0;col<outer.total;col++){
        const value=integral*outer.sourceMatrix[row*outer.total+col];
        if(value!==0)weights.set(col,(weights.get(col)??0)+value);
      }
    }
    if(slopes&&k===0){
      // Match the wake's initial source to the two new surface TE limits.
      // Changing only the left Hermite derivative adds (1-4t+3t²) Δd_TE;
      // its whole-segment integral is zero, so TE/downstream mass is exact.
      const d=slopes[cut.element],delta=addMassWeights(new Map(d[0]),d.at(-1));
      for(let col=0;col<outer.total;col++){
        const old=outer.sourceMatrix[3*w.body.first*outer.total+col]+outer.sourceMatrix[3*w.body.last*outer.total+col];
        if(old!==0)delta.set(col,(delta.get(col)??0)-old);
      }
      const primitive=t=>t-2*t*t+t**3;
      addMassWeights(weights,delta,p.length*(primitive(hi)-primitive(lo)));
    }
    wakeSources.set(cut.face,displacementSource(weights,options));
  }
  return{bodySources,wakeSources,bodyPanels};
}
