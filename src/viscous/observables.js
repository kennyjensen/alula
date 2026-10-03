// SPDX-License-Identifier: GPL-2.0-or-later
import { surfaceTransition } from './trips.js';
import { isentropicState } from '../potential/isentropic.js';
import { pressureForces } from '../inviscid/pressure-forces.js';
import { wakeMomentumDrag } from './wake-drag.js';

// Shared dimensional output for a converged BL assembly. This function does
// not decide convergence or certify forces; the public adapter applies those
// gates using its complete outer-flow and boundary-layer residual.
export function viscousObservables({system,states,st,x},{elements,alpha,mach,referenceChord,momentReference,inviscidReference=mach===0}){
  const {outer,kernel}=system;
  const mass=s=>isentropicState(s.ue,0,{mach}).rho*s.ue*s.deltaStar;
  const surfaces=[],wakes=[],outputElements=[];
  let cl=0,cm=0,cd=0,cdf=0,pressureIntegralDrag=0,maxRelativeSurfaceDisplacement=0,maxRelativeWakeDisplacement=0;
  const scalePoint=p=>({x:p.x*referenceChord,y:p.y*referenceChord});
  const station=(id,point,regime)=>{
    const s=states[id],zero=system.zeroNodes.includes(id),properties=zero?{cf:0,rho:isentropicState(0,0,{mach}).rho,hk:s.deltaStar/s.theta}:kernel.station(s,regime);
    return{index:id,...scalePoint(point),s:s.s*referenceChord,ue:s.ue,
      theta:s.theta*referenceChord,deltaStar:s.deltaStar*referenceChord,h:(s.deltaStar-(s.wakeGap??0))/s.theta,hk:properties.hk,
      ...(s.wakeGap===undefined?{}:{wakeGap:s.wakeGap*referenceChord,viscousDeltaStar:(s.deltaStar-s.wakeGap)*referenceChord,totalH:s.deltaStar/s.theta}),
      cf:properties.cf*properties.rho*s.ue**2,amplification:regime==='laminar'?s.aux:null,ctau:regime==='laminar'?null:s.aux,regime};
  };
  outer.bodies.forEach((body,e)=>{
    const finite=body.trailingEdge?.kind==='finite-base';
    const te=finite?body.baseGeometry?.center:body.points[0];
    if(!te||![te.x,te.y].every(Number.isFinite)||finite&&(!Number.isFinite(body.baseGeometry.width)||body.baseGeometry.width<=0))
      throw new Error('Finite-base observables require the physical trailing-edge center and positive projected gap.');
    const chord=Math.max(...body.points.map(p=>Math.hypot(p.x-te.x,p.y-te.y)));
    const le=body.points.reduce((a,p)=>Math.hypot(p.x-te.x,p.y-te.y)>Math.hypot(a.x-te.x,a.y-te.y)?p:a);
    const chx=te.x-le.x,chy=te.y-le.y,chord2=chx*chx+chy*chy;
    const localX=p=>((p.x-le.x)*chx+(p.y-le.y)*chy)/chord2;
    const cf=new Float64Array(body.points.length);
    for(const surf of system.activeSurfaces(st).filter(s=>s.body===body)){
      const stations=surf.ids.map((id,k)=>{
        const j=id-body.start,regime=k<surf.transition?'laminar':'turbulent';
        const result=station(id,body.points[j],regime);cf[j]=result.cf;
        maxRelativeSurfaceDisplacement=Math.max(maxRelativeSurfaceDisplacement,states[id].deltaStar/chord);
        const a=outer.panels[Math.max(body.first,body.first+j-1)],b=outer.panels[Math.min(body.last,body.first+j)];
        const tx=a.tx+b.tx,ty=a.ty+b.ty,length=Math.hypot(tx,ty);
        result.displacement={x:result.x+result.deltaStar*ty/length,y:result.y-result.deltaStar*tx/length};
        return result;
      });
      const transition=surfaceTransition(system,surf,states,st[e].s),{point}=transition;
      surfaces.push({element:e,name:elements[e].name??`Element ${e+1}`,side:surf.side===1?'upper':'lower',
        transition:localX(point),transitionPoint:scalePoint(point),forced:transition.forced,stations});
    }
    const q=Float64Array.from({length:body.points.length},(_,j)=>x[4*(body.start+j)+3]);
    const cpValues=Array.from(q,v=>isentropicState(v,0,{mach}).cp);
    const forceOptions={alpha,momentOrigin:{x:momentReference.x/referenceChord,y:momentReference.y/referenceChord}};
    let forces=pressureForces(body.points,cpValues,forceOptions),baseInfo=null,solidCp=null;
    if(finite){
      const solid=body.solidPoints,{upperIndex,lowerIndex}=body.trailingEdge,n=solid?.length-1;
      if(!Array.isArray(solid)||n<3||![upperIndex,lowerIndex].every(i=>Number.isInteger(i)&&i>=0&&i<n)
        ||solid[0].x!==solid[n].x||solid[0].y!==solid[n].y||body.points.length!==(lowerIndex-upperIndex+n)%n+1)
        throw new Error('Finite-base observables require the complete retained solid and matching surface indices.');
      const sourceNodes=new Map();
      body.points.forEach((p,j)=>{
        const index=(upperIndex+j)%n;
        if(p.x!==solid[index].x||p.y!==solid[index].y)throw new Error('Finite-base surface and solid coordinates differ.');
        sourceNodes.set(index,j);
      });
      const base=Array.from({length:(upperIndex-lowerIndex+n)%n+1},(_,j)=>solid[(lowerIndex+j)%n]);
      const baseCp=.5*(cpValues[0]+cpValues.at(-1));
      // A repeated lower corner represents the one-sided surface/base
      // pressure values. Its zero-length edge contributes no force; every
      // subsequent base edge is an original retained solid segment.
      const full=pressureForces([...body.points,...base],[...cpValues,...base.map(()=>baseCp)],forceOptions);
      baseInfo={normalGap:body.baseGeometry.width*referenceChord,points:base.map(scalePoint),pressure:baseCp,
        pressureModel:'Mean TE endpoint Cp on every retained base segment; no base skin friction',
        pressureForces:Object.fromEntries(['cx','cy','cl','cd','cm'].map(k=>[k,full[k]-forces[k]]))};
      forces=full;
      const inviscidBase=inviscidReference ? .5*((1-outer.q0[body.start]**2)+(1-outer.q0[body.end]**2)) : null;
      solidCp=solid.map((p,i)=>{
        const j=sourceNodes.get(i%n),baseNode=j===undefined;
        return{...scalePoint(p),cp:baseNode?baseCp:cpValues[j],
          ...(inviscidReference?{cpInviscid:baseNode?inviscidBase:1-outer.q0[body.start+j]**2}:{}),
          qt:baseNode?null:q[j],base:baseNode};
      });
    }
    cl+=forces.cl;cm+=forces.cm;pressureIntegralDrag+=forces.cd;
    for(let j=0;j<body.points.length-1;j++){
      const p=outer.panels[body.first+j],direction=Math.sign(q[j]+q[j+1]);
      cdf+=.5*(cf[j]+cf[j+1])*p.length*direction*(p.tx*Math.cos(alpha*Math.PI/180)+p.ty*Math.sin(alpha*Math.PI/180));
    }
    const cp=body.points.map((p,j)=>({...scalePoint(p),cp:cpValues[j],...(inviscidReference?{cpInviscid:1-outer.q0[body.start+j]**2}:{}),qt:q[j]}));
    outputElements.push({name:elements[e].name??`Element ${e+1}`,points:(finite?body.solidPoints:body.points).map(scalePoint),cp:solidCp??cp,
      cl:forces.cl,cm:forces.cm,...(finite?{finiteBase:baseInfo}:{})});
    const w=outer.wakes[e],stations=w.points.map((p,j)=>station(w.start+j,p,'wake')),end=states[w.end];
    for(let id=w.start;id<=w.end;id++)maxRelativeWakeDisplacement=Math.max(maxRelativeWakeDisplacement,states[id].deltaStar/chord);
    if((end.wakeGap??0)>0)throw Object.assign(new Error('Wake exit lies inside the finite-base dead-air closure; extend the wake before estimating drag.'),
      {code:'WAKE_EXIT_GAP_OPEN',element:e,wakeGap:end.wakeGap*referenceChord});
    const drag=wakeMomentumDrag(end,{mach}).cd;cd+=drag;
    const upper=states[body.start],lower=states[body.end],first=states[w.start];
    const gap=first.wakeGap??0,totalMassDifference=mass(first)-mass(upper)-mass(lower),gapMass=isentropicState(first.ue,0,{mach}).rho*first.ue*gap;
    wakes.push({element:e,name:outputElements[e].name,stations,drag,
      matching:{momentum:first.theta-upper.theta-lower.theta,displacement:first.deltaStar-upper.deltaStar-lower.deltaStar-gap,
        massFlux:totalMassDifference-gapMass,...(finite?{geometricGap:gap,geometricMassFlux:gapMass,totalMassFluxDifference:totalMassDifference,
          units:'Lengths and mass-deficit lengths normalized by reference chord; edge density and speed normalized by freestream'}:{})}});
  });
  return{surfaces,wakes,outputElements,cl,cm,cd,cdf,pressureIntegralDrag,maxRelativeSurfaceDisplacement,maxRelativeWakeDisplacement};
}
