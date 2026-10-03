// SPDX-License-Identifier: GPL-2.0-or-later
// Incompressible Euler limit: a shared vortex/source system for every body.
// Surface inputs are signed Ue*delta*, ordered along the CCW contour; wake
// inputs are positive downstream. No cross-element influence is discarded.
import { prepareContour, validateAssembly, pointInside, segmentsTouch } from '../geometry/airfoil.js';
import { makePanel, sourceVelocity } from './panel.js';
import { vortexBasis, velocityAt } from './linear-vortex.js';
import { factorLinear, normInf, linearResidual } from '../numerics/linear.js';
import { createContourTopology } from '../geometry/contour-topology.js';
import { createContourCurve } from '../geometry/contour-curve.js';
import { panelTrailingEdgeGeometry } from '../viscous/context.js';
import { vortexStreamfunctionBasis } from './streamfunction.js';
import { basePanelStreamfunctionBasis } from './finite-base-influence.js';
import { assemblyFiniteBaseSourceChart } from './finite-base-source-chart.js';
import { prepareAirfoilElement } from '../geometry/airfoil-element.js';

const dot = (a,b) => a.u*b.tx+a.v*b.ty;

// Complex potential of a polynomial source. A common angle chart on each
// target body prevents a source's streamfunction branch cut crossing another
// element. Real parts are branch independent and give vortex streamfunctions.
export function polynomialPotential(point,p,chart=Math.PI/2){
  const x=((point.x-p.a.x)*p.tx+(point.y-p.a.y)*p.ty)/p.length;
  const y=(-(point.x-p.a.x)*p.ty+(point.y-p.a.y)*p.tx)/p.length;
  const r0=Math.hypot(x,y),r1=Math.hypot(x-1,y),phi=Math.atan2(p.ty,p.tx);
  const wrap=(angle,center)=>angle+2*Math.PI*Math.round((center-angle)/(2*Math.PI));
  let t0=wrap(Math.atan2(y,x),chart),t1=wrap(Math.atan2(y,x-1),t0);
  if(r0>2){
    const zr=x/(x*x+y*y),zi=-y/(x*x+y*y),result=[];
    for(let k=0;k<3;k++){
      let real=Math.log(r0)/(k+1),imag=t0/(k+1),pr=1,pi=0;
      for(let n=1;n<=60;n++){
        const next=pr*zr-pi*zi;pi=pr*zi+pi*zr;pr=next;
        const weight=1/(n*(k+n+1));real-=pr*weight;imag-=pi*weight;
        if(Math.hypot(pr,pi)*weight<Number.EPSILON*1e-2)break;
      }
      result.push({real:p.length*(real+Math.log(p.length)/(k+1))/(2*Math.PI),imag:p.length*(imag+phi/(k+1))/(2*Math.PI)});
    }
    return result;
  }
  if(r0<1e-14)t1=wrap(Math.PI,chart);
  const result=[];
  let ir=Math.log(r0/r1),ii=t0-t1;
  for(let k=0;k<3;k++){
    let fr,fi;
    if(r0<1e-14){fr=-1/(k+1)**2;fi=t1/(k+1);}
    else if(r1<1e-14){let harmonic=0;for(let j=1;j<=k+1;j++)harmonic+=1/j;fr=-harmonic/(k+1);fi=t0/(k+1);}
    else{const next=x*ir-y*ii-1/(k+1);ii=x*ii+y*ir;ir=next;fr=(Math.log(r1)+ir)/(k+1);fi=(t1+ii)/(k+1);}
    result.push({real:p.length*(fr+Math.log(p.length)/(k+1))/(2*Math.PI),imag:p.length*(fi+phi/(k+1))/(2*Math.PI)});
  }
  return result;
}

// Analytic integral of source density a+b*xi+c*xi^2 on a straight segment.
// The complex recurrence I_k=z*I_(k-1)-1/k also gives the endpoint finite
// parts. At wake vertices, common logarithmic singularities cancel between
// adjacent segments; each endpoint uses log(distance/chord) with chord=1.
export function polynomialSourceBasis(point, p, { principal = false,selfMidpoint=false } = {}) {
  // A known self midpoint uses its exact local coordinates. Subtracting large
  // global coordinates can otherwise put it on either side of its own sheet,
  // selecting a spurious half-source normal-velocity jump.
  const x = selfMidpoint ? .5 : ((point.x-p.a.x)*p.tx+(point.y-p.a.y)*p.ty)/p.length;
  const y = selfMidpoint?0:(-(point.x-p.a.x)*p.ty+(point.y-p.a.y)*p.tx)/p.length;
  if(Math.hypot(x,y)>2){
    const zr=x/(x*x+y*y),zi=-y/(x*x+y*y),result=[];
    for(let k=0;k<3;k++){
      let re=0,im=0,pr=1,pi=0;
      for(let n=0;n<60;n++){
        const next=pr*zr-pi*zi;pi=pr*zi+pi*zr;pr=next;
        re+=pr/(k+n+1);im+=pi/(k+n+1);
        if(Math.hypot(pr,pi)/(k+n+1)<Number.EPSILON*1e-2*Math.hypot(re,im))break;
      }
      result.push({u:(re*p.tx+im*p.ty)/(2*Math.PI),v:(re*p.ty-im*p.tx)/(2*Math.PI)});
    }
    return result;
  }
  const atA = Math.hypot(point.x-p.a.x,point.y-p.a.y) < 1e-14*p.length;
  const atB = Math.hypot(point.x-p.b.x,point.y-p.b.y) < 1e-14*p.length;
  let re, im;
  if (atA) { re = -Math.log(p.length); im = 0; }
  else if (atB) { re = Math.log(p.length); im = 0; }
  else {
    re = .5*Math.log((x*x+y*y)/((x-1)**2+y*y));
    im = principal && Math.abs(y)<1e-12 ? 0 : -Math.atan2(y,x*(x-1)+y*y);
  }
  const result = [];
  for (let k = 0; k < 3; k++) {
    result.push({u:(re*p.tx+im*p.ty)/(2*Math.PI),v:(re*p.ty-im*p.tx)/(2*Math.PI)});
    const next = x*re-y*im-1/(k+1); im = x*im+y*re; re = next;
  }
  return result;
}

export function createDisplacementOperator({elements,alpha=0,wakeLength=2,wakeLengths,wakeCount=24,wakePaths,wakeInitialization='straight',computeInfluence=true}) {
  if (!Number.isFinite(alpha) || !Number.isFinite(wakeLength) || wakeLength<=0 || !Number.isInteger(wakeCount) || wakeCount<4) throw new Error('Invalid displacement operator controls.');
  if(!['straight','inviscid'].includes(wakeInitialization))throw new Error('Invalid wake initialization.');
  if(typeof computeInfluence!=='boolean')throw new Error('Invalid displacement influence control.');
  elements=elements.map(prepareAirfoilElement);
  const topologies=elements.map(e=>e.trailingEdge?.kind==='finite-base'?createContourTopology(e.points,e):null);
  const solids=elements.map((e,k)=>topologies[k]?.points??prepareContour(e.points)); validateAssembly(solids);
  const contours=solids.map((p,k)=>topologies[k]?.surface.points??p),hasFiniteBase=topologies.some(Boolean);
  if(wakeLengths&&(!Array.isArray(wakeLengths)||wakeLengths.length!==elements.length||wakeLengths.some(v=>!Number.isFinite(v)||v<=0)))throw new Error('Supply one positive dimensional wake length per element.');
  const panels=[], bodies=[],basePanels=[]; let ns=0;
  for (let e=0;e<contours.length;e++) {
    const points=contours[e], start=ns, first=panels.length, s=[0];
    for(let j=0;j<points.length-1;j++) { const p={...makePanel(points[j],points[j+1],e),node:start+j}; panels.push(p); s.push(s.at(-1)+p.length); }
    ns+=points.length; const body={start,end:ns-1,first,last:panels.length-1,points,s};
    if(topologies[e]){
      const baseGeometry=panelTrailingEdgeGeometry(points),t=baseGeometry.tangentDerivative,length=Math.hypot(t.x,t.y);
      if(!(baseGeometry.width>0&&length>0))throw new Error('A finite trailing edge requires positive projected width and a resolved downstream direction.');
      const chart=assemblyFiniteBaseSourceChart(topologies[e],{x:t.x/length,y:t.y/length},solids.filter((_,k)=>k!==e)),cutDirection=chart.direction;
      Object.assign(body,{baseGeometry,solidPoints:solids[e],trailingEdge:structuredClone(elements[e].trailingEdge)});
      body.basePanels=topologies[e].base.panels.map(b=>{
        const p=makePanel(b.start,b.end,e);
        return{...p,sourcePanelIndex:b.sourcePanelIndex,upperNode:start,lowerNode:body.end,
          sourceCoefficient:-.5*(t.x*p.ty-t.y*p.tx),vortexCoefficient:-.5*(t.x*p.tx+t.y*p.ty),
          cutDirection,cutOrigin:{...chart.origin}};
      });
      basePanels.push(...body.basePanels);
    }
    bodies.push(body);
  }
  const nodes=contours.flat(), n=panels.length, size=ns+bodies.length;
  // The external BL assembly needs geometry and conservative source weights,
  // but its potential/Euler closure supplies velocity. Do not factor and
  // solve a dense panel-displacement operator that it will never evaluate.
  // Cold inviscid wake tracing still requires the ordinary inviscid field.
  const needInviscid=computeInfluence||(!wakePaths&&wakeInitialization==='inviscid');
  const a=needInviscid?new Float64Array(size*size):null,rad=alpha*Math.PI/180,u=Math.cos(rad),v=Math.sin(rad),rhs=needInviscid?new Float64Array(size):null;
  const special=new Set(bodies.filter(b=>!b.baseGeometry).map(b=>b.end));
  const baseVelocity=(point,p)=>{
    const q=vortexBasis(point,p),vortex={u:q[0].u+q[1].u,v:q[0].v+q[1].v};
    return{u:p.sourceCoefficient*vortex.v+p.vortexCoefficient*vortex.u,
      v:-p.sourceCoefficient*vortex.u+p.vortexCoefficient*vortex.v};
  };
  if(needInviscid)for(let i=0;i<ns;i++)if(!special.has(i)){
    const p=nodes[i];rhs[i]=-u*p.y+v*p.x;
    for(const q of panels){
      if(hasFiniteBase){const f=vortexStreamfunctionBasis(p,q);a[i*size+q.node]+=f[0];a[i*size+q.node+1]+=f[1];}
      else{const f=polynomialPotential(p,q);a[i*size+q.node]-=f[0].real-f[1].real;a[i*size+q.node+1]-=f[1].real;}
    }
    for(const q of basePanels){
      const f=basePanelStreamfunctionBasis(p,q);
      if(f.onSourceCut&&q.sourceCoefficient!==0)throw new Error('A surface node lies on a finite-base source cut.');
      const value=q.sourceCoefficient*f.source+q.vortexCoefficient*f.vortex;
      a[i*size+q.upperNode]+=value;a[i*size+q.lowerNode]-=value;
    }
  }
  bodies.forEach((b,e)=>{
    b.centroid={x:b.points.reduce((s,p)=>s+p.x,0)/b.points.length,y:b.points.reduce((s,p)=>s+p.y,0)/b.points.length};
    if(needInviscid)for(let i=b.start;i<b.end+Number(Boolean(b.baseGeometry));i++)a[i*size+ns+e]=-1;
    if(b.baseGeometry){
      if(needInviscid){a[(ns+e)*size+b.start]=1;a[(ns+e)*size+b.end]=1;}
      return;
    }
    // Replace the redundant repeated-TE equation by zero interior velocity
    // along its bisector, including the source-sheet contribution below.
    let tx=-panels[b.first].tx+panels[b.last].tx,ty=-panels[b.first].ty+panels[b.last].ty;
    if(hasFiniteBase){
      const curve=createContourCurve(b.points),upper=curve.evaluate(0).derivative,lower=curve.evaluate(curve.length).derivative;
      tx=-upper.x/Math.hypot(upper.x,upper.y)+lower.x/Math.hypot(lower.x,lower.y);
      ty=-upper.y/Math.hypot(upper.x,upper.y)+lower.y/Math.hypot(lower.x,lower.y);
    }
    const length=Math.hypot(tx,ty);tx/=length;ty/=length;
    const distance=.1*Math.min(panels[b.first].length,panels[b.last].length);
    b.bisector={x:b.points[0].x-distance*tx,y:b.points[0].y-distance*ty,tx,ty};
    if(needInviscid){
      rhs[b.end]=-u*tx-v*ty;
      for(const q of panels)vortexBasis(b.bisector,q).forEach((v,k)=>{a[b.end*size+q.node+k]+=dot(v,b.bisector);});
      for(const q of basePanels){const value=dot(baseVelocity(b.bisector,q),b.bisector);a[b.end*size+q.upperNode]+=value;a[b.end*size+q.lowerNode]-=value;}
      a[(ns+e)*size+b.start]=1;a[(ns+e)*size+b.end]=1;
    }
  });
  const solve=needInviscid?factorLinear(a):null,inviscid=solve?.(rhs),gamma=inviscid?.slice(0,ns);
  if(needInviscid)for(const p of basePanels){const d=gamma[p.upperNode]-gamma[p.lowerNode];p.sourceStrength=p.sourceCoefficient*d;p.vortexStrength=p.vortexCoefficient*d;}
  const field=needInviscid?{panels,gamma,u,v,...(hasFiniteBase?{basePanels}: {})}:null;
  let total=ns;
  // Closing the viscous wake shape remains a separate equation. Either a
  // straight or inviscid-streamline path can supply its starting geometry.
  const wakes=bodies.map((b,e)=>{
    const center=b.baseGeometry?.center??b.points[0],chord=Math.max(...b.points.map(p=>Math.hypot(p.x-center.x,p.y-center.y)));
    const extent=wakeLengths?.[e]??wakeLength*chord;
    let tx=-panels[b.first].tx+panels[b.last].tx,ty=-panels[b.first].ty+panels[b.last].ty;
    if(b.baseGeometry){tx=b.baseGeometry.tangentDerivative.x;ty=b.baseGeometry.tangentDerivative.y;}
    const dl=Math.hypot(tx,ty);tx/=dl;ty/=dl;
    // XYWAKE starts just downstream of the TE midpoint but assigns zero
    // wake arc there. Scale its 1e-4 offset with chord for reference covariance.
    const offset=b.baseGeometry?1e-4*chord:0,p0=offset?{x:center.x+offset*tx,y:center.y+offset*ty}:center;
    if(b.baseGeometry&&pointInside(p0,b.solidPoints))throw new Error('The native-style first wake point lies inside the retained base.');
    const points=[],s=[],start=total;
    const initial=.5*(panels[b.first].length+panels[b.last].length);
    let lo=1,hi=4;
    const length=r=>initial*(r===1?wakeCount-1:(r**(wakeCount-1)-1)/(r-1));
    if (length(lo)>extent) throw new Error('Wake count is too large for its specified initial spacing.');
    while(length(hi)<extent)hi*=2;
    for(let k=0;k<70;k++){const mid=(lo+hi)/2;if(length(mid)>extent)hi=mid;else lo=mid;}
    const ratio=(lo+hi)/2;let distance=0;
    for(let j=0;j<wakeCount;j++){
      s.push(distance);
      if(j===0||wakeInitialization==='straight'||wakePaths)points.push({x:p0.x+tx*distance,y:p0.y+ty*distance});
      else{
        const step=distance-s[j-1],previous=points.at(-1),v0=j===1?{u:tx,v:ty}:velocityAt(previous,field),speed0=Math.hypot(v0.u,v0.v);
        if(speed0<1e-10)throw new Error('Inviscid wake initialization reaches stagnation.');
        const mid={x:previous.x+.5*step*v0.u/speed0,y:previous.y+.5*step*v0.v/speed0},vm=velocityAt(mid,field),speed=Math.hypot(vm.u,vm.v);
        if(speed<1e-10)throw new Error('Inviscid wake initialization reaches stagnation.');
        points.push({x:previous.x+step*vm.u/speed,y:previous.y+step*vm.v/speed});
      }
      distance+=initial*ratio**j;
    }
    total+=wakeCount;return{points,s,tx,ty,start,end:total-1,body:b};
  });
  if(wakePaths&&(!Array.isArray(wakePaths)||wakePaths.length!==wakes.length))throw new Error('Supply one wake path per element.');
  wakes.forEach((w,e)=>{
    if(wakePaths){
      const path=wakePaths[e];
      if(!Array.isArray(path)||path.length!==wakeCount||!path.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.y)))throw new Error('Invalid wake path.');
      if(Math.hypot(path[0].x-w.points[0].x,path[0].y-w.points[0].y)>1e-10)throw new Error('A wake must start at its element trailing edge.');
      w.points=path.map(p=>({...p}));w.s=[0];
      for(let j=1;j<wakeCount;j++)w.s.push(w.s.at(-1)+makePanel(path[j-1],path[j]).length);
    }
    w.segments=w.points.slice(1).map((p,j)=>makePanel(w.points[j],p));
    w.s=[0];for(const p of w.segments)w.s.push(w.s.at(-1)+p.length);
    w.tangents=w.points.map((p,j)=>{
      const left=w.segments[Math.max(0,j-1)],right=w.segments[Math.min(j,wakeCount-2)];
      const tx=left.tx+right.tx,ty=left.ty+right.ty,length=Math.hypot(tx,ty);
      if(length<.1)throw new Error('Unresolved reversal in wake geometry.');
      return{tx:tx/length,ty:ty/length};
    });
    for(let j=0;j<w.segments.length;j++){
      const segment=w.segments[j];
      for(const body of bodies){
        if(pointInside(segment,body.solidPoints??body.points))throw new Error('A wake intersects an element.');
        for(let k=body.first;k<=body.last;k++){
          if(j===0&&body===w.body&&(k===body.first||k===body.last))continue;
          if(segmentsTouch(segment.a,segment.b,panels[k].a,panels[k].b,1e-11))throw new Error('A wake intersects an element.');
        }
        for(const p of body.basePanels??[])if(segmentsTouch(segment.a,segment.b,p.a,p.b,1e-11))throw new Error('A wake intersects an element base.');
      }
    }
  });
  for(let e=0;e<wakes.length;e++)for(let f=e+1;f<wakes.length;f++)for(const a of wakes[e].segments)for(const b of wakes[f].segments){
    if(segmentsTouch(a.a,a.b,b.a,b.b,1e-11))throw new Error('Wake confluence needs an explicit merging closure.');
  }
  const sourcePanels=panels.map(p=>({...p,degree:0}));
  for(const w of wakes){w.sourceStart=sourcePanels.length;for(let j=0;j<wakeCount-1;j++)sourcePanels.push({...makePanel(w.points[j],w.points[j+1]),degree:2,midpoint:w.segments[j]});}
  const nq=sourcePanels.length, nc=3*nq;
  // Maps nodal displacement flux to polynomial source coefficients. Hermite
  // interpolation of mass makes sigma=dm/ds continuous in each wake and
  // exactly conservative over every segment, including nonuniform spacing.
  const sourceMatrix=new Float64Array(nc*total);
  const set=(row,col,value)=>{sourceMatrix[row*total+col]+=value;};
  for(let j=0;j<n;j++){const p=panels[j];set(3*j,p.node,-1/p.length);set(3*j,p.node+1,1/p.length);}
  for(const w of wakes){
    const d=Array.from({length:wakeCount},()=>new Float64Array(total));
    // TE source matching. The two surface mass derivatives add in downstream
    // coordinates, while their mass fluxes already carry the contour sign.
    if(w.body.baseGeometry){
      // PSWLIN uses the first wake mass secant at its upstream endpoint.
      // Distinct finite-TE branches must not share a forced source slope.
      const h=w.s[1]-w.s[0];d[0][w.start]=-1/h;d[0][w.start+1]=1/h;
    }else for(let col=0;col<total;col++)d[0][col]=sourceMatrix[3*w.body.first*total+col]+sourceMatrix[3*w.body.last*total+col];
    for(let j=1;j<wakeCount-1;j++){
      const hl=w.s[j]-w.s[j-1],hr=w.s[j+1]-w.s[j];
      d[j][w.start+j-1]=-hr/(hl*(hl+hr));d[j][w.start+j]=hr/(hl*(hl+hr))-hl/(hr*(hl+hr));d[j][w.start+j+1]=hl/(hr*(hl+hr));
    }
    // The downstream continuation has constant displacement flux: d_end=0.
    for(let j=0;j<wakeCount-1;j++){
      const row=3*(w.sourceStart+j),h=w.s[j+1]-w.s[j];
      for(let col=0;col<total;col++){set(row,col,d[j][col]);set(row+1,col,-4*d[j][col]-2*d[j+1][col]);set(row+2,col,3*d[j][col]+3*d[j+1][col]);}
      set(row+1,w.start+j,-6/h);set(row+1,w.start+j+1,6/h);set(row+2,w.start+j,6/h);set(row+2,w.start+j+1,-6/h);
    }
  }
  if(!computeInfluence)return{bodies,wakes,panels,sourcePanels,sourceMatrix,total,surfaceCount:ns,field,
    diagnostics:{inviscidResidual:needInviscid?normInf(linearResidual(a,inviscid,rhs)):null}};
  const streamSource=new Float64Array(ns*nc);
  for(const body of bodies)for(let i=body.start;i<body.end+Number(Boolean(body.baseGeometry));i++)for(let j=0;j<nq;j++){
    const p=nodes[i],q=sourcePanels[j];
    const dx=body.centroid.x-q.a.x,dy=body.centroid.y-q.a.y;
    const chart=j<n&&q.element===bodies.indexOf(body)?Math.PI/2:Math.atan2(-dx*q.ty+dy*q.tx,dx*q.tx+dy*q.ty);
    const basis=polynomialPotential(p,q,chart);
    for(let k=0;k<=(j<n?0:2);k++)streamSource[i*nc+3*j+k]=basis[k].imag;
  }
  for(const body of bodies.filter(b=>!b.baseGeometry))for(let j=0;j<nq;j++){
    const basis=polynomialSourceBasis(body.bisector,sourcePanels[j]);
    for(let k=0;k<=(j<n?0:2);k++)streamSource[body.end*nc+3*j+k]=dot(basis[k],body.bisector);
  }
  const wakeVortex=new Float64Array((total-ns)*ns),wakeSource=new Float64Array((total-ns)*nc),q0=new Float64Array(total);q0.set(gamma);
  for(const w of wakes)for(let j=0;j<wakeCount;j++){
    const row=w.start+j-ns,p=w.points[j],tangent=w.tangents[j];
    if(j===0){wakeVortex[row*ns+w.body.end]=1;q0[w.start]=gamma[w.body.end];continue;}
    q0[w.start+j]=u*tangent.tx+v*tangent.ty;
    for(const q of panels)vortexBasis(p,q).forEach((b,k)=>{const coeff=dot(b,tangent);wakeVortex[row*ns+q.node+k]+=coeff;q0[w.start+j]+=coeff*gamma[q.node+k];});
    for(const q of basePanels){
      const coeff=dot(baseVelocity(p,q),tangent);
      wakeVortex[row*ns+q.upperNode]+=coeff;wakeVortex[row*ns+q.lowerNode]-=coeff;
      q0[w.start+j]+=coeff*(gamma[q.upperNode]-gamma[q.lowerNode]);
    }
    for(let k=0;k<nq;k++){
      const q=sourcePanels[k],basis=k<n?[sourceVelocity(p,q)]:polynomialSourceBasis(p,q,{principal:true});
      // At a curved wake vertex the common bisector gives equal projections
      // of the adjacent tangents. Their continuous-source logarithmic terms
      // cancel in this tangential principal value. Normal velocity is sampled
      // at segment interiors, where no corner finite part is needed.
      basis.forEach((b,l)=>{wakeSource[row*nc+3*k+l]=dot(b,tangent);});
    }
  }
  const influence=new Float64Array(total*total);
  const bodyMass=new Float64Array(ns*total),wakeMass=new Float64Array((total-ns)*total);
  // A source coefficient depends on at most a handful of neighboring mass
  // values. Use that sparsity before solving the dense body influence system.
  for(let k=0;k<nc;k++)for(let col=0;col<total;col++){
    const coefficient=sourceMatrix[k*total+col];if(coefficient===0)continue;
    for(let i=0;i<ns;i++)bodyMass[i*total+col]-=streamSource[i*nc+k]*coefficient;
    for(let i=0;i<total-ns;i++)wakeMass[i*total+col]+=wakeSource[i*nc+k]*coefficient;
  }
  for(let col=0;col<total;col++){
    const b=new Float64Array(size);
    for(let i=0;i<ns;i++)b[i]=bodyMass[i*total+col];
    const dg=solve(b);
    for(let i=0;i<ns;i++)influence[i*total+col]=dg[i];
    for(let i=ns;i<total;i++){
      let value=wakeMass[(i-ns)*total+col];for(let k=0;k<ns;k++)value+=wakeVortex[(i-ns)*ns+k]*dg[k];
      influence[i*total+col]=value;
    }
  }
  const evaluate=mass=>{
    if(mass.length!==total||!mass.every(Number.isFinite))throw new Error('Invalid displacement flux vector.');
    return q0.map((q,i)=>{for(let j=0;j<total;j++)q+=influence[i*total+j]*mass[j];return q;});
  };
  const sources=mass=>Float64Array.from({length:nc},(_,i)=>{let value=0;for(let j=0;j<total;j++)value+=sourceMatrix[i*total+j]*mass[j];return value;});
  const velocityField=mass=>{
    const q=evaluate(mass),sigma=sources(mass);
    return point=>{
      let vx=u,vy=v;
      for(const p of panels)vortexBasis(point,p).forEach((b,k)=>{vx+=b.u*q[p.node+k];vy+=b.v*q[p.node+k];});
      for(const p of basePanels){const b=baseVelocity(point,p),d=q[p.upperNode]-q[p.lowerNode];vx+=b.u*d;vy+=b.v*d;}
      for(let j=0;j<nq;j++){
        const basis=polynomialSourceBasis(point,sourcePanels[j],{principal:true,selfMidpoint:point===sourcePanels[j].midpoint});
        for(let k=0;k<=(j<n?0:2);k++){vx+=basis[k].u*sigma[3*j+k];vy+=basis[k].v*sigma[3*j+k];}
      }
      return{u:vx,v:vy};
    };
  };
  return{bodies,wakes,panels,sourcePanels,sourceMatrix,total,surfaceCount:ns,q0,influence,evaluate,sources,field,
    velocityField,
    diagnostics:{inviscidResidual:normInf(linearResidual(a,inviscid,rhs))}};
}
