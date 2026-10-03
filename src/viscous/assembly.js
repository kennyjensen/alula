// SPDX-License-Identifier: GPL-2.0-or-later
// Augmented simultaneous BL/outer-flow system, shared by all lifting elements.
import { createDisplacementOperator } from '../inviscid/displacement.js';
import { createIntegralKernel } from './integral.js';
import { runCoupled } from './context.js';
import { geometricTripArc,tripDistance } from './trips.js';
import { solveActiveNewton as solveNewton } from './active-newton.js';
import { solveNewton as solveSmoothNewton } from '../numerics/newton.js';
import { solveBlSchur } from '../numerics/bl-schur.js';
import { remapInitialMach } from './seed.js';
import { createXfoilDeadAirGap } from './xfoil-dead-air-gap.js';
import { initializeSharedPanelBoundaryLayers } from './shared-field-initializer.js';

function createAssembly({ elements, alpha=0, mach=0,reynolds=1e6, ncrit=9,
  trips=[1,1], wakeCount=24, wakeLength=2,wakeLengths,wakePaths,wakeInitialization='straight',initialization='native',seed,jacobianMode='analytic' }={},externalVelocity=false) {
  if(!['auto','native','march','shared-mrchue'].includes(initialization))throw new Error('Invalid BL initialization.');
  if(!['analytic','finite-difference'].includes(jacobianMode))throw new Error('Invalid BL Jacobian mode.');
  if(mach!==0&&!externalVelocity)throw new Error('The displacement coupling backend requires Mach zero; compressible Euler/BL coupling is separate.');
  const outer=createDisplacementOperator({elements,alpha,wakeCount,wakeLength,wakeLengths,wakePaths,wakeInitialization,computeInfluence:!externalVelocity});
  const elementTrips=elements.map(e=>e.trips??trips);
  if(elementTrips.some(pair=>!Array.isArray(pair)||pair.length!==2||pair.some(v=>!Number.isFinite(v)||v<=0||v>1)))throw new Error('Each element needs upper/lower trip fractions in (0, 1].');
  const kernel=createIntegralKernel({reynolds,mach,ncrit,exactJacobian:true});
  const similarityRoot=solveSmoothNewton({initial:[.29,.65],admissible:z=>z[0]>0&&z[1]>z[0],
    residual:z=>kernel.interval({regime:'similarity',downstream:{s:1,ue:1,aux:0,theta:z[0]/Math.sqrt(reynolds),deltaStar:z[1]/Math.sqrt(reynolds)}}).residual.slice(1)});
  if(!similarityRoot.converged)throw new Error('Stagnation similarity initialization failed.');
  const {total,surfaceCount,bodies,wakes}=outer;
  const hasFiniteBase=bodies.some(b=>b.baseGeometry!=null);
  // Geometry is frozen within this BL solve. w.s is distance from the solid
  // TE midpoint, independent of the moving stagnation-origin BL coordinate.
  // Total displacement enters outer mass; only BLPRV subtracts this gap.
  const deadAir=hasFiniteBase?bodies.map(b=>b.baseGeometry?createXfoilDeadAirGap({
    normalGap:b.baseGeometry.width,upperDerivative:b.baseGeometry.upperDerivative,
    lowerDerivative:b.baseGeometry.lowerDerivative,sharp:b.baseGeometry.sharp}):null):null;
  const wakeGaps=hasFiniteBase?new Float64Array(total):null;
  if(hasFiniteBase)wakes.forEach((w,e)=>w.s.forEach((s,j)=>{wakeGaps[w.start+j]=deadAir[e]?.at(s).gap??0;}));
  const thicknessScale=1/Math.sqrt(reynolds), initial=new Float64Array(4*total);
  const surfaces=[], turbulent=new Uint8Array(total),zeroNodes=new Int32Array(bodies.length).fill(-1);
  const sharedRequested=initialization==='shared-mrchue'||(initialization==='auto'&&!externalVelocity&&hasFiniteBase&&bodies.length>1);
  const shared=!seed&&sharedRequested?initializeSharedPanelBoundaryLayers(outer,{reynolds,mach,ncrit,elementTrips,wakeGaps}):null;
  let autoRaw;
  if(initialization==='auto'&&!seed&&!shared){
    autoRaw=bodies.map((body,e)=>runCoupled(body.points,{alpha,reynolds,ncrit,trips:elementTrips[e],maxIterations:100}));
    const admissibleSeed=raw=>{
      if(!raw.converged)return false;
      for(let side=1;side<=2;side++)for(let i=2;i<=raw.bl.IBLTE[side];i++){
        const d=raw.bl.DSTR[i][side]/raw.panel.CHORD;
        if(!Number.isFinite(d)||d<=0||d>.1)return false;
      }
      return true;
    };
    // Keep the same initialization stage on every element. A highly loaded
    // isolated element may have an invalid solution even when neighboring
    // elements make the assembly flow admissible. In that case discard the
    // isolated roots and start from the native first-step states throughout.
    if(!autoRaw.every(admissibleSeed))autoRaw=bodies.map((body,e)=>runCoupled(body.points,{alpha,reynolds,ncrit,trips:elementTrips[e],maxIterations:1}));
  }
  let guesses=seed||shared?null:autoRaw??bodies.map((body,e)=>runCoupled(body.points,{alpha,reynolds,ncrit,trips:elementTrips[e],maxIterations:initialization==='march'?1:100}));
  const finiteGuess=raw=>{
    const b=raw.bl;
    for(let side=1;side<=2;side++)for(let i=2;i<=b.NBL[side];i++)if(![b.CTAU[i][side],b.THET[i][side],b.DSTR[i][side],b.UEDG[i][side]].every(Number.isFinite)
      ||b.THET[i][side]<=0||b.DSTR[i][side]-(hasFiniteBase&&i>b.IBLTE[side]?(b.WGAP[i-b.IBLTE[side]]??0):0)<=b.THET[i][side]||b.UEDG[i][side]<0)return false;
    return true;
  };
  // Native TRCHEK cannot initialize a trip before its first laminar station.
  // An untripped first-step state is a finite starting guess; the requested
  // material trips and all simultaneous equations are still used below.
  if(guesses&&!guesses.every(finiteGuess))guesses=bodies.map(body=>runCoupled(body.points,{alpha,reynolds,ncrit,trips:[1,1],maxIterations:1}));
  if(guesses&&!guesses.every(finiteGuess))throw new Error('No finite boundary-layer initialization for this geometry.');
  const seedKey=JSON.stringify({points:bodies.map(b=>b.points),alpha,mach,reynolds,ncrit,trips,wakeCount,...(elements.some(e=>e.trips)?{elementTrips}:{}),
    ...(hasFiniteBase?{finiteBases:bodies.map(b=>b.baseGeometry?{solidPoints:b.solidPoints,trailingEdge:b.trailingEdge,baseGeometry:b.baseGeometry}:null)}:{})});
  if(shared){initial.set(shared.initial);turbulent.set(shared.turbulent);zeroNodes.set(shared.zeroNodes);surfaces.push(...shared.surfaces);}
  else if(seed){
    if(seed.key!==seedKey||seed.x.length!==initial.length)throw new Error('Coupled warm-start geometry or conditions changed.');
    initial.set(seed.x);turbulent.set(seed.turbulent);zeroNodes.set(seed.zeroNodes);
    if(hasFiniteBase){
      if(!seed.wakeGaps||seed.wakeGaps.length!==total||!Array.from(seed.wakeGaps).every(g=>Number.isFinite(g)&&g>=0))
        throw new Error('Finite-base BL warm start requires its previous prescribed wake gaps.');
      // Wake relaxation changes physical arc distances. Transfer the fluid
      // thickness and restore current dead air, never scale the solid base.
      for(let id=surfaceCount;id<total;id++)if(seed.wakeGaps[id]!==wakeGaps[id])
        initial[4*id+2]=(initial[4*id+2]*thicknessScale-seed.wakeGaps[id]+wakeGaps[id])/thicknessScale;
    }
    surfaces.push(...seed.surfaces.map(s=>({...s,ids:[...s.ids],body:bodies[s.element]})));
  }else for(let e=0;e<bodies.length;e++){
    const body=bodies[e],wake=wakes[e];
    const raw=guesses[e];
    const {bl}=raw;
    for(let side=1;side<=2;side++){
      const ids=[];
      for(let i=2;i<=bl.IBLTE[side];i++){
        const id=body.start+bl.IPAN[i][side]-1;ids.push(id);
        initial[4*id]=bl.CTAU[i][side];initial[4*id+1]=bl.THET[i][side]/thicknessScale;initial[4*id+2]=bl.DSTR[i][side]/thicknessScale;
        initial[4*id+3]=(side===1?-1:1)*bl.UEDG[i][side];
        turbulent[id]=i>=bl.ITRAN[side]?1:0;
      }
      const transition=bl.ITRAN[side]-2;
      // Store trip as a geometric arc coordinate, so its distance from
      // stagnation moves together with every other BL station.
      const tripContext={...bl,XSTRIP:[...bl.XSTRIP]};tripContext.XSTRIP[side]=elementTrips[e][side-1];
      const tripArc=geometricTripArc(tripContext,side);
      surfaces.push({body,side,ids,transition,transitionId:ids[transition],tripArc,forced:elementTrips[e][side-1]<1});
    }
    // An exact symmetric stagnation node carries no displacement flux. Its
    // BL thickness is the common limit of the two adjacent similarity states;
    // do not pass Ue=0 into logarithmic interval formulae or invent a floor.
    const nearest=body.s.reduce((best,s,i)=>Math.abs(s-bl.SST)<Math.abs(body.s[best]-bl.SST)?i:best,0),id=body.start+nearest;
    // This is an initialization hypothesis based on the isolated native
    // stagnation point. Neighboring elements may give this node a nonzero
    // velocity; every solver releases it and solves the ordinary equations
    // before returning success. Requiring the assembly's initial inviscid
    // velocity to be zero leaves tiny native U/s states ill-conditioned.
    if(elementTrips[e][0]===elementTrips[e][1]&&Math.abs(body.s[nearest]-bl.SST)<1e-8){
      zeroNodes[e]=id;initial[4*id+3]=0;
      initial[4*id]=0;for(let k=1;k<=2;k++)initial[4*id+k]=.5*(initial[4*(id-1)+k]+initial[4*(id+1)+k]);
      for(const surf of surfaces.filter(s=>s.body===body))if(surf.ids[0]===id){surf.ids.shift();surf.transition--;}
    }
    const w=wakes[e],te=bl.IBLTE[2];
    const nativeS=[];for(let i=te+1;i<=bl.NBL[2];i++)nativeS.push(bl.XSSI[i][2]-bl.XSSI[te+1][2]);
    for(let j=0;j<wakeCount;j++){
      const id=w.start+j,s=w.s[j];let k=0;while(k<nativeS.length-2&&nativeS[k+1]<s)k++;
      const f=Math.min(1,(s-nativeS[k])/(nativeS[k+1]-nativeS[k]));
      const interp=key=>bl[key][te+1+k][2]*(1-f)+bl[key][te+2+k][2]*f;
      // XICALC supplies total native DSTR. Interpolate its fluid part, then
      // restore the prescribed gap at this assembly's actual wake distance.
      const delta=body.baseGeometry?((bl.DSTR[te+1+k][2]-bl.WGAP[1+k])*(1-f)
        +(bl.DSTR[te+2+k][2]-bl.WGAP[2+k])*f+ wakeGaps[id]):interp('DSTR');
      initial[4*id]=interp('CTAU');initial[4*id+1]=interp('THET')/thicknessScale;initial[4*id+2]=delta/thicknessScale;initial[4*id+3]=interp('UEDG');turbulent[id]=1;
    }
  }
  if(!seed&&!shared&&mach!==0){
    if(hasFiniteBase){
      const fluid=initial.slice();for(let id=surfaceCount;id<total;id++)fluid[4*id+2]-=wakeGaps[id]/thicknessScale;
      initial.set(remapInitialMach(fluid,0,mach));
      for(let id=surfaceCount;id<total;id++)initial[4*id+2]+=wakeGaps[id]/thicknessScale;
    }else initial.set(remapInitialMach(initial,0,mach));
  }
  const topology=x=>bodies.map((body,e)=>{
    if(zeroNodes[e]>=0){const id=zeroNodes[e];return{left:id-1,right:id+1,s:body.s[id-body.start],derivatives:[0,0],zero:id};}
    let j=body.start;while(j<body.end&&x[4*(j+1)+3]<0)j++;
    if(j===body.start||j>=body.end)throw new Error('Unresolved leading-edge stagnation point.');
    const ua=-x[4*j+3],ub=x[4*(j+1)+3],a=body.s[j-body.start],b=body.s[j+1-body.start];
    return{left:j,right:j+1,s:(ub*a+ua*b)/(ua+ub),distances:[ua*(b-a)/(ua+ub),ub*(b-a)/(ua+ub)],derivatives:[-ub*(b-a)/(ua+ub)**2,-ua*(b-a)/(ua+ub)**2]};
  });
  const activeSurfaces=st=>surfaces.map(surf=>{
    const e=bodies.indexOf(surf.body),ids=[];
    if(surf.side===1)for(let id=st[e].left;id>=surf.body.start;id--)ids.push(id);
    else for(let id=st[e].right;id<=surf.body.end;id++)ids.push(id);
    return{...surf,ids,transition:ids.indexOf(surf.transitionId)};
  });
  const similarity=(state,id,e,st,x)=>{
    if(mach!==0)return state;
    const b=bodies[e],t=st[e];
    // At Mach zero the native similarity residual depends on Ue/s only.
    // Evaluate that exact ratio without subtracting nearly equal arc lengths
    // or cancelling its two enormous partial derivatives near stagnation.
    const slope=t.zero!==undefined?state.ue/state.s:(x[4*t.right+3]-x[4*t.left+3])/(b.s[t.right-b.start]-b.s[t.left-b.start]);
    return{...state,s:1,ue:slope};
  };
  const decode=x=>{
    const st=topology(x),states=Array(total);
    for(let e=0;e<bodies.length;e++){
      const b=bodies[e];for(let id=b.start;id<=b.end;id++)states[id]={s:Math.abs(b.s[id-b.start]-st[e].s),aux:x[4*id],theta:x[4*id+1]*thicknessScale,deltaStar:x[4*id+2]*thicknessScale,ue:Math.abs(x[4*id+3])};
      if(st[e].distances){states[st[e].left].s=st[e].distances[0];states[st[e].right].s=st[e].distances[1];}
      const w=wakes[e];for(let id=w.start;id<=w.end;id++)states[id]={s:b.s.at(-1)-st[e].s+w.s[id-w.start],aux:x[4*id],theta:x[4*id+1]*thicknessScale,deltaStar:x[4*id+2]*thicknessScale,ue:x[4*id+3],...(b.baseGeometry?{wakeGap:wakeGaps[id]}:{})};
    }
    return{states,st};
  };
  // A material trip may lie between stagnation and the first contour node.
  // There is then no resolved laminar upstream station for TRDIF. Supply a
  // virtual station from the same leading-edge similarity solution, halfway
  // to the actual trip. Its velocity follows the linear stagnation segment.
  // The trip is not moved to a mesh node or extrapolated upstream of X1.
  const leadingTransitionInput=(downstream,tripS)=>{
    const s=.5*tripS,slope=downstream.ue/downstream.s,ue=s*slope;
    if(!(tripS>0&&tripS<=downstream.s*(1+1e-12)))throw new Error('Leading transition must precede the first BL station.');
    let z=similarityRoot.x.map(v=>v/Math.sqrt(slope));
    if(mach!==0){
      const root=solveSmoothNewton({initial:z,admissible:a=>a[0]>0&&a[1]>a[0],
        residual:a=>kernel.interval({regime:'similarity',downstream:{s,ue,aux:0,theta:a[0]*thicknessScale,deltaStar:a[1]*thicknessScale}}).residual.slice(1)});
      if(!root.converged)throw new Error('Leading transition similarity initialization failed.');
      z=root.x;
    }
    return{upstream:{s,ue,aux:0,theta:z[0]*thicknessScale,deltaStar:z[1]*thicknessScale},downstream,regime:'transition',tripS};
  };
  let continuation=1;
  const initialMass=Float64Array.from({length:total},(_,i)=>initial[4*i+3]*initial[4*i+2]*thicknessScale);
  const initialOuter=externalVelocity?null:outer.evaluate(initialMass),continuationOffset=initialOuter?.map((q,i)=>initial[4*i+3]-q);
  const setContinuation=value=>{if(!Number.isFinite(value)||value<0||value>1)throw new Error('Invalid continuation parameter.');continuation=value;};
  const admissible=x=>{
    if(!x.every(Number.isFinite))return false;
    for(let i=0;i<total;i++)if(!(x[4*i+1]>0&&x[4*i+2]-(wakeGaps?.[i]??0)/thicknessScale>x[4*i+1])||(turbulent[i]&&x[4*i]<=0)||(i>=surfaceCount&&x[4*i+3]<=0))return false;
    try{
      const {states,st}=decode(x);
      for(let e=0;e<bodies.length;e++)for(let id=bodies[e].start;id<=bodies[e].end;id++){
        if(id===zeroNodes[e])continue;
        if(!(x[4*id+3]*(id<=st[e].left?-1:1)>0))return false;
      }
      return states.every((s,i)=>zeroNodes.includes(i)||(s.s>0&&s.ue>0));
    }catch{return false;}
  };
  const residual=x=>{
    const {states,st}=decode(x),r=new Float64Array(4*total);
    for(const surf of activeSurfaces(st)){
      const e=bodies.indexOf(surf.body),tripS=tripDistance(surf,st[e].s);
      for(let k=0;k<surf.ids.length;k++){
        const id=surf.ids[k],downstream=k===0?similarity(states[id],id,e,st,x):states[id],upstream=k===0?downstream:states[surf.ids[k-1]];
        const regime=k===0?'similarity':k<surf.transition?'laminar':k===surf.transition?'transition':'turbulent';
        const result=kernel.interval(k===0&&surf.transition===0?leadingTransitionInput(states[id],tripS):{upstream,downstream,regime,tripS});
        for(let row=0;row<3;row++)r[4*id+row]=result.residual[row]*(row===0&&turbulent[id]?20:1);
      }
    }
    for(const id of zeroNodes)if(id>=0){
      r[4*id]=x[4*id];
      for(let k=1;k<=2;k++)r[4*id+k]=x[4*id+k]-.5*(x[4*(id-1)+k]+x[4*(id+1)+k]);
    }
    for(const w of wakes)for(let j=0;j<wakeCount;j++){
      const id=w.start+j,downstream=states[id];
      const result=j===0?kernel.trailingEdge(states[w.body.start],states[w.body.end],downstream,w.body.baseGeometry?.width??0):kernel.interval({upstream:states[id-1],downstream,regime:'wake'});
      for(let row=0;row<3;row++)r[4*id+row]=result.residual[row]*(row===0?20:j===0?1/thicknessScale:1);
    }
    if(!externalVelocity){
      const mass=Float64Array.from(states,(s,i)=>x[4*i+3]*s.deltaStar),q=outer.evaluate(mass);
      for(let i=0;i<total;i++)r[4*i+3]=x[4*i+3]-q[i]-(1-continuation)*continuationOffset[i];
    }
    return r;
  };
  const jacobian=x=>{
    const {states,st}=decode(x),size=4*total,jac=new Float64Array(size*size);
    const keys=['aux','theta','deltaStar','ue','s'];
    // Native analytic blocks include the optional exact lag correction.
    // Transition retains local numerical derivatives. The full numerical
    // mode remains an independent audit of these local blocks and chains.
    const block=(id,upId,regime,tripS,e)=>{
      const input={upstream:states[upId],downstream:regime==='similarity'?similarity(states[id],id,e,st,x):states[id],regime,tripS};
      const shapeLog=hasFiniteBase?Math.log(((states[id].deltaStar-(states[id].wakeGap??0))/states[id].theta-1)
        /((states[upId].deltaStar-(states[upId].wakeGap??0))/states[upId].theta-1))
        :Math.log((states[id].deltaStar/states[id].theta-1)/(states[upId].deltaStar/states[upId].theta-1));
      // Keep independent numerical differentiation for transition and for the
      // native extreme-shape upwind cap, whose derivative is approximate.
      const analytic=jacobianMode==='analytic'&&regime!=='transition'&&shapeLog**2<15?kernel.interval(input):null;
      const stableSimilarity=regime==='similarity'&&mach!==0&&st[e].zero===undefined;
      if(stableSimilarity){
        // F depends on the stagnation slope a=U/s and on U^2. Assemble
        // derivatives in those coordinates before applying the two velocity
        // chains. Separate F_U and F_s terms diverge near a stagnation node
        // and lose precision when cancelled in the global matrix.
        const t=st[e],b=bodies[e],ds=b.s[t.right-b.start]-b.s[t.left-b.start],u=states[id].ue,a=(x[4*t.right+3]-x[4*t.left+3])/ds;
        const base=analytic??kernel.interval(input),at=(velocity,slope)=>kernel.interval({...input,downstream:{...states[id],ue:velocity,s:velocity/slope}}).residual;
        let slopeDerivative=base.downstream.map(row=>-row[4]*states[id].s/a);
        if(!analytic){const h=Math.cbrt(Number.EPSILON)*a,p=at(u,a+h),m=at(u,a-h);slopeDerivative=p.map((v,i)=>(v-m[i])/(2*h));}
        let velocityDerivative=base.downstream.map(row=>row[3]+row[4]/a);
        if(u<1e-3||!analytic){
          const z=u*u,h=Math.cbrt(Number.EPSILON)*Math.max(z,1e-3),lo=Math.max(z-h,z*.5),hi=z+h;
          const p=at(Math.sqrt(hi),a),m=at(Math.sqrt(lo),a);
          velocityDerivative=p.map((v,i)=>2*u*(v-m[i])/(hi-lo));
        }
        for(let row=0;row<3;row++){
          const scale=row===0&&turbulent[id]?20:1,d=scale*slopeDerivative[row]/ds;
          jac[(4*id+row)*size+4*t.left+3]-=d;jac[(4*id+row)*size+4*t.right+3]+=d;
          jac[(4*id+row)*size+4*id+3]+=scale*velocityDerivative[row]*Math.sign(x[4*id+3]);
        }
      }
      for(const side of ['upstream','downstream']){
        if(regime==='similarity'&&side==='upstream')continue;
        const colId=side==='upstream'?upId:id;
        for(let k=0;k<5;k++){
          if(stableSimilarity&&(k===3||k===4))continue;
          if(regime==='similarity'&&k===4&&mach===0)continue;
          const key=keys[k],value=input[side][key],h=Math.cbrt(Number.EPSILON)*Math.max(Math.abs(value),k===0?.01:k<3?1e-7:1e-6);
          const plus={...input,[side]:{...input[side],[key]:value+h}},minus={...input,[side]:{...input[side],[key]:value-h}};
          let rp,rm,width=2*h;
          if(!analytic){
            try{rp=kernel.interval(plus).residual;}catch{}
            try{rm=kernel.interval(minus).residual;}catch{}
            if(!rp&&!rm)throw new Error('No admissible local BL derivative.');
            if(!rp){rp=kernel.interval(input).residual;width=h;}
            if(!rm){rm=kernel.interval(input).residual;width=h;}
          }
          for(let row=0;row<3;row++){
            const scale=row===0&&turbulent[id]?20:1,d=(analytic?analytic[side][row][k]:(rp[row]-rm[row])/width)*scale;
            if(regime==='similarity'&&k===3&&mach===0){
              const t=st[e],b=bodies[e];
              if(t.zero!==undefined)jac[(4*id+row)*size+4*id+3]+=d*Math.sign(x[4*id+3])/states[id].s;
              else{const ds=b.s[t.right-b.start]-b.s[t.left-b.start];jac[(4*id+row)*size+4*t.left+3]-=d/ds;jac[(4*id+row)*size+4*t.right+3]+=d/ds;}
            }
            else if(k<4){const factor=k===1||k===2?thicknessScale:k===3?Math.sign(x[4*colId+3]):1;jac[(4*id+row)*size+4*colId+k]+=d*factor;}
            else{
              const direction=colId<surfaceCount?(colId<=st[e].left?-1:1):1;
              for(let a=0;a<2;a++)jac[(4*id+row)*size+4*(a===0?st[e].left:st[e].right)+3]-=d*direction*st[e].derivatives[a];
            }
          }
        }
      }
      if(regime==='transition'&&tripS<Number.MAX_VALUE){
        const h=Math.cbrt(Number.EPSILON)*Math.max(tripS,1e-6);
        let rp,rm,width=2*h;
        try{rp=kernel.interval({...input,tripS:tripS+h}).residual;}catch{}
        try{rm=kernel.interval({...input,tripS:tripS-h}).residual;}catch{}
        if(!rp&&!rm)throw new Error('No admissible transition-trip derivative.');
        if(!rp){rp=kernel.interval(input).residual;width=h;}
        if(!rm){rm=kernel.interval(input).residual;width=h;}
        const d=rp.map((v,i)=>(v-rm[i])/width);
        const direction=id<=st[e].left?-1:1;
        for(let row=0;row<3;row++)for(let a=0;a<2;a++)jac[(4*id+row)*size+4*(a===0?st[e].left:st[e].right)+3]-=d[row]*(row===0?20:1)*direction*st[e].derivatives[a];
      }
    };
    for(const surf of activeSurfaces(st)){
      const e=bodies.indexOf(surf.body),tripS=tripDistance(surf,st[e].s);
      surf.ids.forEach((id,k)=>{
        if(k!==0||surf.transition!==0){block(id,k===0?id:surf.ids[k-1],k===0?'similarity':k<surf.transition?'laminar':k===surf.transition?'transition':'turbulent',tripS,e);return;}
        // Differentiate the composite TRDIF + virtual similarity state. This
        // includes velocity, stagnation-distance and material-trip chains.
        const evaluate=(state,trip)=>kernel.interval(leadingTransitionInput(state,trip)).residual;
        const base=evaluate(states[id],tripS);
        for(const [l,key] of [...keys,'tripS'].entries()){
          const value=key==='tripS'?tripS:states[id][key],h=Math.cbrt(Number.EPSILON)*Math.max(Math.abs(value),l===0?.01:l<3?1e-7:1e-12);
          let p,m,width=2*h;
          try{p=key==='tripS'?evaluate(states[id],value+h):evaluate({...states[id],[key]:value+h},tripS);}catch{}
          try{m=key==='tripS'?evaluate(states[id],value-h):evaluate({...states[id],[key]:value-h},tripS);}catch{}
          if(!p&&!m)throw new Error('No admissible leading-transition derivative.');
          if(!p){p=base;width=h;}if(!m){m=base;width=h;}
          for(let row=0;row<3;row++){
            const d=(p[row]-m[row])/width*(row===0?20:1);
            if(l<4)jac[(4*id+row)*size+4*id+l]+=d*(l===1||l===2?thicknessScale:l===3?Math.sign(x[4*id+3]):1);
            else for(let a=0;a<2;a++)jac[(4*id+row)*size+4*(a===0?st[e].left:st[e].right)+3]-=d*(surf.side===1?-1:1)*st[e].derivatives[a];
          }
        }
      });
    }
    for(const id of zeroNodes)if(id>=0){
      jac[4*id*size+4*id]=1;
      for(let k=1;k<=2;k++){jac[(4*id+k)*size+4*id+k]=1;jac[(4*id+k)*size+4*(id-1)+k]=-.5;jac[(4*id+k)*size+4*(id+1)+k]=-.5;}
    }
    wakes.forEach((w,e)=>{
      const upper=states[w.body.start],lower=states[w.body.end],sum=upper.theta+lower.theta,aux=(upper.theta*upper.aux+lower.theta*lower.aux)/sum;
      for(const id of[w.body.start,w.body.end]){
        jac[4*w.start*size+4*id]=-20*states[id].theta/sum;
        jac[4*w.start*size+4*id+1]=-20*(states[id].aux-aux)/sum*thicknessScale;
        jac[(4*w.start+1)*size+4*id+1]=-1;jac[(4*w.start+2)*size+4*id+2]=-1;
      }
      jac[4*w.start*size+4*w.start]=20;jac[(4*w.start+1)*size+4*w.start+1]=1;jac[(4*w.start+2)*size+4*w.start+2]=1;
      for(let id=w.start+1;id<=w.end;id++)block(id,id-1,'wake',Number.MAX_VALUE,e);
    });
    if(!externalVelocity)for(let i=0;i<total;i++){
      for(let j=0;j<total;j++){
        const d=outer.influence[i*total+j];jac[(4*i+3)*size+4*j+2]=-d*x[4*j+3]*thicknessScale;jac[(4*i+3)*size+4*j+3]=-d*states[j].deltaStar;
      }
      jac[(4*i+3)*size+4*i+3]+=1;
    }
    return jac;
  };
  const stagnationKey=t=>2*t.left+Number(t.zero!==undefined);
  let activeStagnation=topology(initial).map(stagnationKey);
  const snapshot=()=>({transitions:surfaces.map(s=>s.transitionId),turbulent:turbulent.slice(),stagnation:activeStagnation.slice(),zeroNodes:zeroNodes.slice()});
  const restore=saved=>{surfaces.forEach((s,i)=>{s.transitionId=saved.transitions[i];});turbulent.set(saved.turbulent);activeStagnation=saved.stagnation.slice();zeroNodes.set(saved.zeroNodes);};
  const updateActive=x=>{
    const {states,st}=decode(x);let changed=false;
    st.forEach((t,e)=>{
      if(stagnationKey(t)===activeStagnation[e])return;
      changed=true;activeStagnation[e]=stagnationKey(t);
      for(const id of[t.left,t.right]){
        const slope=mach===0?similarity(states[id],id,e,st,x).ue:states[id].ue/states[id].s;
        x[4*id]=0;turbulent[id]=0;x[4*id+1]=similarityRoot.x[0]/Math.sqrt(slope);x[4*id+2]=similarityRoot.x[1]/Math.sqrt(slope);
        if(mach!==0){
          const root=solveSmoothNewton({initial:[x[4*id+1],x[4*id+2]],admissible:z=>z[0]>0&&z[1]>z[0],
            residual:z=>kernel.interval({regime:'similarity',downstream:{...states[id],aux:0,theta:z[0]*thicknessScale,deltaStar:z[1]*thicknessScale}}).residual.slice(1)});
          if(!root.converged)throw new Error('Compressible stagnation similarity initialization failed.');
          x[4*id+1]=root.x[0];x[4*id+2]=root.x[1];
        }
        states[id]={...states[id],aux:0,theta:x[4*id+1]*thicknessScale,deltaStar:x[4*id+2]*thicknessScale};
      }
    });
    for(const [index,surf] of activeSurfaces(st).entries()){
      const e=bodies.indexOf(surf.body),tripS=tripDistance(surf,st[e].s),amplification=[0];
      let transition=tripS>0&&tripS<=states[surf.ids[0]].s?0:surf.ids.length-1;
      for(let k=1;transition!==0&&k<surf.ids.length;k++){
        const interval={upstream:{...states[surf.ids[k-1]],aux:amplification[k-1]},downstream:states[surf.ids[k]],tripS};
        let check;
        try{check=kernel.transitionCheck(interval);}catch(error){throw Object.assign(error,{interval});}
        amplification[k]=check.amplification;
        if(check.transition){transition=k;break;}
      }
      const target=surf.ids[transition];
      if(target!==surfaces[index].transitionId)changed=true;
      surfaces[index].transitionId=target;
      surf.ids.forEach((id,k)=>{
        const next=k>=transition?1:0;
        if(next!==turbulent[id]){
          changed=true;
          x[4*id]=next?kernel.station({...states[id],aux:.03},'turbulent').transitionShear:amplification[k];
          if(next&&k===0&&transition===0){
            const root=solveSmoothNewton({initial:[x[4*id],x[4*id+1],x[4*id+2]],maxIterations:30,
              admissible:z=>z[0]>0&&z[1]>0&&z[2]>z[1],
              residual:z=>kernel.interval(leadingTransitionInput({...states[id],aux:z[0],theta:z[1]*thicknessScale,deltaStar:z[2]*thicknessScale},tripS)).residual});
            if(root.converged){x.set(root.x,4*id);states[id]={...states[id],aux:root.x[0],theta:root.x[1]*thicknessScale,deltaStar:root.x[2]*thicknessScale};}
          }
        }
        // Keep the amplification state consistent with the criterion used to
        // choose the transition interval. A stale N can place TRCHEK's root in
        // a different interval and make every otherwise valid trial fail.
        if(!next)x[4*id]=amplification[k];
        turbulent[id]=next;
      });
    }
    return changed;
  };
  const releaseNonzeroStagnation=x=>{
    let changed=false;
    for(let e=0;e<zeroNodes.length;e++)if(zeroNodes[e]>=0&&Math.abs(x[4*zeroNodes[e]+3])>1e-8){zeroNodes[e]=-1;changed=true;}
    if(changed)updateActive(x);
    return changed;
  };
  const linearHistory=[];
  const exportSeed=x=>({...(shared||seed?.initialization?{initialization:{method:shared?.diagnostics.method??seed.initialization.method,flowSolved:false,isolatedSolves:0}}:{}),...(hasFiniteBase?{wakeGaps:wakeGaps.slice()}:{}),key:seedKey,x:x.slice(),turbulent:turbulent.slice(),zeroNodes:zeroNodes.slice(),
    surfaces:activeSurfaces(topology(x)).map(s=>({element:bodies.indexOf(s.body),side:s.side,ids:s.ids,
      transition:s.transition,transitionId:s.transitionId,tripArc:s.tripArc,forced:s.forced}))});
  const guardedTrialState=(x,delta,step)=>{
    const candidate=x.map((v,i)=>v+step*delta[i]);
    // Local fraction-to-boundary safeguards affect trial updates only. A wake
    // near H=1 must not freeze every other element's Newton correction. The
    // unmodified residual still decides acceptance and final convergence.
    for(let id=0;id<total;id++){
      const i=4*id,gap=(wakeGaps?.[id]??0)/thicknessScale,excess=x[i+2]-gap-x[i+1];
      candidate[i+2]=Math.max(candidate[i+2],candidate[i+1]+gap+.2*excess);
      if(turbulent[id])candidate[i]=Math.max(candidate[i],.2*x[i]);
    }
    return candidate;
  };
  const linearSolve=(matrix,rhs,x)=>{
    const order=activeSurfaces(topology(x)).flatMap(s=>s.ids);
    for(const id of zeroNodes)if(id>=0)order.push(id);
    for(const wake of wakes)for(let id=wake.start;id<=wake.end;id++)order.push(id);
    const result=solveBlSchur(matrix,rhs,order);
    linearHistory.push({method:result.method,backwardError:result.backwardError,fallbackReason:result.fallbackReason});
    return result.x;
  };
  return{...(shared?{initialization:shared.diagnostics}:seed?.initialization?{initialization:{...seed.initialization,resumed:true}}:{}),...(hasFiniteBase?{baseGeometry:bodies.map(b=>b.baseGeometry??null),wakeGaps}:{}),outer,kernel,total,initial,surfaces,activeSurfaces,turbulent,zeroNodes,thicknessScale,decode,residual,jacobian,admissible,setContinuation,updateActive,snapshot,restore,releaseNonzeroStagnation,linearSolve:externalVelocity?undefined:linearSolve,linearHistory,guardedTrialState,exportSeed,externalVelocity,leadingTransitionInput};
}

export const createCoupledAssembly=input=>createAssembly(input);

// An unclosed BL equation block for an external Euler/potential system. The
// three BL rows per station are populated; every fourth row is deliberately
// zero and must be supplied by the outer-flow velocity closure. This factory
// does not provide a standalone flow solve or publish aerodynamic forces.
export const createBoundaryLayerAssembly=input=>createAssembly(input,true);

export function solveCoupledAssembly(input,{onIteration,maxIterations=40,tolerance=1e-8,continuation=false,initialState,activeSet}={}){
  const system=createCoupledAssembly(input);
  if(initialState){if(initialState.length!==system.initial.length)throw new Error('Invalid coupled initial state dimensions.');system.initial.set(initialState);}
  if(activeSet)system.restore(activeSet);
  let solution;
  if(!continuation){
    const saved=system.snapshot();
    solution=solveNewton({...system,maxIterations,tolerance,onIteration});
    const remaining=maxIterations-solution.history.filter(h=>h.iteration>0).length;
    if(!solution.converged&&remaining>=3){
      system.restore(saved);
      const retry=solveNewton({...system,trialState:system.guardedTrialState,maxIterations:remaining,tolerance,
        onIteration:h=>onIteration?.({...h,attempt:2})});
      solution={...retry,history:[...solution.history,...retry.history.map(h=>({...h,attempt:2}))]};
    }
  }
  else{
    let value=0,step=.25,x=system.initial;const history=[];
    system.setContinuation(0);
    solution=solveNewton({...system,initial:x,maxIterations,tolerance,onIteration:h=>onIteration?.({...h,continuation:0})});
    history.push(...solution.history.map(h=>({...h,continuation:0})));
    x=solution.x;
    while(solution.converged&&value<1){
      const next=Math.min(1,value+step);system.setContinuation(next);
      const active=system.snapshot();
      let attempt;
      try{attempt=solveNewton({...system,initial:x,maxIterations,tolerance,onIteration:h=>onIteration?.({...h,continuation:next})});}
      catch(error){attempt={converged:false,x,history:[],reason:error.message};}
      history.push(...attempt.history.map(h=>({...h,continuation:next})));
      if(attempt.converged){solution=attempt;x=attempt.x;value=next;step=Math.min(.25,step*1.5);}
      else if(step>1/1024){system.restore(active);step/=2;}
      else{solution={...attempt,reason:`Continuation stopped at ${value}: ${attempt.reason}`};break;}
    }
    system.setContinuation(1);
    solution={...solution,history,converged:solution.converged&&value===1,continuation:value};
  }
  for(let restart=0;restart<system.outer.bodies.length&&solution.converged;restart++){
    if(!system.releaseNonzeroStagnation(solution.x))break;
    const next=solveNewton({...system,initial:solution.x,maxIterations,tolerance,onIteration});
    solution={...next,history:[...solution.history,...next.history]};
  }
  return{system,...solution,...system.decode(solution.x)};
}
