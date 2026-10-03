// SPDX-License-Identifier: GPL-2.0-or-later
// Prescribed shared-panel-field MRCHUE startup. This changes the initial
// guess only; the simultaneous BL/source/wake equations remain unchanged.
import {createPanelContext} from './context.js';
import {createIntegralKernel} from './integral.js';
import {geometricTripArc,tripDistance} from './trips.js';
import {mrchue,blpini} from './xfoil/xbl.js';
import {ensureCtx,blprv,blkin,blvar,blmid,trchek,tesys,blsys,hkin,syncComToVars} from './xfoil/xblsys.js';

// Also independently testable against the existing prescribed-profile
// original-Fortran fixture. Geometry is needed only for nonterminal trips.
export function initializePrescribedPanelProfiles(profiles,parameters={}, {geometry,wake=[],normalGap=0,trips=[1,1]}={}){
 const {reynolds,mach,gamma,ncrit,velocityConvention}=createIntegralKernel(parameters).parameters;
 if(velocityConvention!=='physical')throw new Error('Shared MRCHUE requires physical edge velocities.');
 if(!Array.isArray(profiles)||profiles.length!==2||profiles.some(p=>!Array.isArray(p)||p.length<2||p.some((a,i)=>
  ![a.s,a.ue].every(Number.isFinite)||a.s<=0||a.ue<=0||(i>0&&a.s<=p[i-1].s))))throw new Error('Supply two ordered positive surface profiles.');
 if(!Array.isArray(trips)||trips.length!==2||trips.some(t=>!Number.isFinite(t)||t<=0||t>1)||!Number.isFinite(normalGap)||normalGap<0)
  throw new Error('Invalid prescribed MRCHUE trip or base width.');
 if(trips.some(t=>t<1)&&!geometry)throw new Error('Material trips require complete panel geometry.');
 if(wake.some((p,i)=>![p.s,p.ue,p.wakeGap].every(Number.isFinite)||p.ue<=0||p.wakeGap<0
  ||(i===0?p.s!==profiles[1].at(-1).s:p.s<=wake[i-1].s)))throw new Error('Invalid prescribed merged-wake stations.');
 const gm1=gamma-1,totalRatio=1+.5*gm1*mach**2,hstinv=gm1*mach**2/totalRatio,hstinvMs=gm1/totalRatio**2;
 const herat=1-.5*hstinv,heratMs=-.5*hstinvMs,hvrat=.35,factor=herat**1.5*(1+hvrat)/(herat+hvrat),messages=[];
 const ctx={...geometry,QINFBL:1,TKBL:0,TKBL_MS:0,RSTBL:totalRatio**(1/gm1),RSTBL_MS:.5*totalRatio**(1/gm1-1),
  HSTINV:hstinv,HSTINV_MS:hstinvMs,REYBL:reynolds*factor,REYBL_RE:factor,
  REYBL_MS:reynolds*factor*(1.5/herat-1/(herat+hvrat))*heratMs,
  GAMBL:gamma,GM1BL:gm1,HVRAT:hvrat,AMCRIT:ncrit,BULE:1,IDAMPV:0,ANTE:normalGap,
  blprv,blkin,blvar,blmid,trchek,tesys,blsys,hkin,syncComToVars,log:m=>messages.push(m)};
 const size=Math.max(profiles[0].length,profiles[1].length+wake.length)+2;
 for(const name of ['XSSI','UEDG','CTAU','THET','DSTR','MASS','TAU','DIS','CTQ','DELT','TSTR'])ctx[name]=Array.from({length:size},()=>new Float64Array(3));
 for(const name of ['NBL','IBLTE','ITRAN','XSSITR','TFORCE'])ctx[name]=new Float64Array(3);
 ctx.ACRIT=[0,ncrit,ncrit];ctx.XSTRIP=[0,...trips];ctx.WGAP=Float64Array.from([0,...wake.map(p=>p.wakeGap)]);
 profiles.forEach((profile,side)=>{const is=side+1;ctx.NBL[is]=ctx.IBLTE[is]=profile.length+1;
  profile.forEach((p,i)=>{ctx.XSSI[i+2][is]=p.s;ctx.UEDG[i+2][is]=p.ue;});});
 ctx.NBL[2]+=wake.length;
 wake.forEach((p,j)=>{const i=ctx.IBLTE[2]+1+j;ctx.XSSI[i][2]=p.s;ctx.UEDG[i][2]=p.ue;});
 ensureCtx(ctx);blpini(ctx);mrchue(ctx);
 const read=(is,i,p)=>({s:p.s,ue:ctx.UEDG[i][is],aux:ctx.CTAU[i][is],theta:ctx.THET[i][is],deltaStar:ctx.DSTR[i][is],
  ...(p.wakeGap===undefined?{}:{wakeGap:p.wakeGap})});
 const surfaces=profiles.map((p,side)=>({states:p.map((a,i)=>read(side+1,i+2,a)),transition:ctx.ITRAN[side+1]-2,
  s:ctx.XSSITR[side+1],forced:Boolean(ctx.TFORCE[side+1]),targetHK:p.map((_,i)=>ctx.HTARG[side+1][i+2])}));
 const wakeStates=wake.map((p,j)=>read(2,ctx.IBLTE[2]+1+j,p));
 if(surfaces.some((s,k)=>!Number.isInteger(s.transition)||s.transition<1||s.transition>=profiles[k].length)
  ||[...surfaces.flatMap(s=>s.states),...wakeStates].some(p=>!Object.values(p).every(Number.isFinite)||p.ue<=0||p.theta<=0||p.deltaStar-(p.wakeGap??0)<=p.theta))
  throw Object.assign(new Error('Shared prescribed MRCHUE returned an unusable initial guess.'),{diagnostics:{messages,surfaces,wakeStates}});
 return{surfaces,wakeStates,messages,localConvergenceWarnings:messages.filter(m=>m.includes('Convergence failed')),
  method:'shared-field MRCHUE',flowSolved:false};
}

export function initializeSharedPanelBoundaryLayers(outer,{reynolds=1e6,mach=0,ncrit=9,elementTrips,wakeGaps}={}){
 const {bodies,wakes,total,q0}=outer;
 if(!q0||q0.length!==total||!q0.every(Number.isFinite))throw new Error('Shared MRCHUE requires the complete assembly inviscid edge field.');
 const scale=1/Math.sqrt(reynolds),initial=new Float64Array(4*total),turbulent=new Uint8Array(total),zeroNodes=new Int32Array(bodies.length).fill(-1),surfaces=[],diagnostics=[];
 for(let e=0;e<bodies.length;e++){
  const b=bodies[e],w=wakes[e],trips=elementTrips[e],zeros=[];for(let id=b.start;id<=b.end;id++)if(q0[id]===0)zeros.push(id);
  let left=b.start;while(left<b.end&&q0[left+1]<0)left++;
  if(zeros.length>1)throw new Error('Shared edge field has multiple exact stagnation nodes.');
  let right=left+1,sst,distances;
  if(zeros.length){const id=zeros[0];left=id-1;right=id+1;zeroNodes[e]=id;sst=b.s[id-b.start];}
  else{const ua=-q0[left],ub=q0[right],a=b.s[left-b.start],z=b.s[right-b.start];
   sst=(ub*a+ua*z)/(ua+ub);distances=[ua*(z-a)/(ua+ub),ub*(z-a)/(ua+ub)];}
  if(left<=b.start||right>=b.end||!Number.isFinite(sst))throw new Error('Unresolved shared-field leading-edge stagnation point.');
  for(let id=b.start;id<=b.end;id++)if(id!==zeroNodes[e]&&!(q0[id]*(id<=left?-1:1)>0))throw new Error('Shared edge field has an additional surface reversal.');
  const panel=createPanelContext(b.points,0,{geometryOnly:true}),geometry={N:b.points.length,SST:sst,SLE:panel.SLE,XLE:panel.XLE,YLE:panel.YLE,XTE:panel.XTE,YTE:panel.YTE};
  for(const key of ['X','Y','S','W1','W2','W3','W4'])geometry[key]=new Float64Array(b.points.length+1);
  b.points.forEach((p,i)=>{geometry.X[i+1]=p.x;geometry.Y[i+1]=p.y;geometry.S[i+1]=b.s[i];});geometry.XSTRIP=[0,...trips];
  const maps=[Array.from({length:left-b.start+1},(_,i)=>left-i),Array.from({length:b.end-right+1},(_,i)=>right+i)];
  const descriptors=maps.map((ids,side)=>({body:b,side:side+1,ids,tripArc:geometricTripArc(geometry,side+1),forced:trips[side]<1}));
  const profiles=maps.map((ids,side)=>ids.map((id,k)=>({id,s:k===0&&distances?distances[side]:Math.abs(b.s[id-b.start]-sst),ue:Math.abs(q0[id])})));
  descriptors.forEach((d,k)=>{const trip=tripDistance(d,sst);if(!(trip>0))throw new Error('Material trip coincides with shared stagnation.');
   if(trip<=profiles[k][0].s){const a=profiles[k][0],s=.5*trip;profiles[k].unshift({id:null,s,ue:a.ue*s/a.s});}});
  const wake=w.s.map((arc,j)=>({s:profiles[1].at(-1).s+arc,ue:q0[w.start+j],wakeGap:wakeGaps?.[w.start+j]??0}));
  const native=initializePrescribedPanelProfiles(profiles,{reynolds,mach,ncrit},{geometry,wake,normalGap:b.baseGeometry?.width??0,trips});
  descriptors.forEach((d,k)=>{const r=native.surfaces[k],transitionId=profiles[k][r.transition].id;
   if(transitionId===null)throw new Error('Shared MRCHUE transition falls at an unresolved virtual station.');
   d.transitionId=transitionId;d.transition=d.ids.indexOf(transitionId);surfaces.push(d);
   profiles[k].forEach((p,j)=>{if(p.id===null)return;const a=r.states[j],id=p.id;initial.set([a.aux,a.theta/scale,a.deltaStar/scale,(k===0?-1:1)*a.ue],4*id);turbulent[id]=j>=r.transition?1:0;});});
  native.wakeStates.forEach((a,j)=>{const id=w.start+j;initial.set([a.aux,a.theta/scale,a.deltaStar/scale,a.ue],4*id);turbulent[id]=1;});
  if(zeroNodes[e]>=0){const id=zeroNodes[e];for(let k=1;k<=2;k++)initial[4*id+k]=.5*(initial[4*(id-1)+k]+initial[4*(id+1)+k]);}
  diagnostics.push({element:e,stagnation:{left,right,s:sst,zero:zeroNodes[e]},surfaces:descriptors.map((d,k)=>({side:d.side,tripArc:d.tripArc,tripS:tripDistance(d,sst),transition:d.transition,transitionId:d.transitionId,
   virtualStation:profiles[k][0].id===null,requestedProfile:profiles[k],states:native.surfaces[k].states,targetHK:native.surfaces[k].targetHK})),wakeStates:native.wakeStates,
   messages:native.messages,localConvergenceWarnings:native.localConvergenceWarnings});
 }
 return{initial,turbulent,zeroNodes,surfaces,diagnostics:{method:'shared-field MRCHUE',flowSolved:false,isolatedSolves:0,elements:diagnostics}};
}
