// SPDX-License-Identifier: GPL-2.0-or-later
import { createBoundaryLayerAssembly } from '../viscous/assembly.js';
import { solveActiveNewton } from '../viscous/active-newton.js';
import { meshDisplacementSources } from '../viscous/mass.js';
import { multielementMesh } from '../euler/multielement-mesh.js';
import { triangularAirfoilMesh,deformTriangularWakes } from '../euler/triangular-mesh.js';
import { airfoilCirculation } from './airfoil.js';
import { affine,sumAffine,evaluateAffine } from './reconstruction.js';
import { partitionPotentialInterfaces,potentialFluxConstraint } from './interfaces.js';
import { createPotentialSystem,solvePotential } from './system.js';
import { isentropicState } from './isentropic.js';
import { solvePotentialLinear } from './linear.js';
import { solveLinear } from '../numerics/linear.js';
import { remapSeedMach } from '../viscous/seed.js';
import { solveMultielementViscous } from '../viscous/multielement.js';
import { interpolatePanelAverages,integratedPanelVelocity } from './edge-sampling.js';
import { createMultipoleFarfield } from './farfield.js';

// Research assembly: cell potentials, circulations, wall normal velocities,
// shared wake potentials, and all BL variables enter the same residual.
// Public browser promotion additionally requires spatial/physical validation
// and closure of the mean wake geometry; a small algebraic residual alone is
// insufficient. No pressure correction or independently marched final BL.
export function createCoupledPotential(input,{mesh:meshOptions={},sparse=true,order=2,meshSeed,fluxQuadrature=3,edgeVelocitySampling='vertex',edgeVelocityFraction=1,farfield='fixed',bodyMassInterpolation='linear',onMesh}={}){
  if(!['vertex','panel-average','station-average'].includes(edgeVelocitySampling))throw new Error('Unknown BL edge-velocity sampling method.');
  if(!['fixed','multipole'].includes(farfield))throw new Error('Unknown farfield condition.');
  if(!Number.isFinite(edgeVelocityFraction)||edgeVelocityFraction<0||edgeVelocityFraction>1)throw new Error('Invalid edge-velocity continuation fraction.');
  const bl=createBoundaryLayerAssembly(input),{outer}=bl,mach=input.mach??0;
  const {type='cross-line',...controls}=meshOptions,contours=outer.bodies.map(b=>b.points),wakePaths=outer.wakes.map(w=>w.points);
  if(!['cross-line','triangular'].includes(type))throw new Error('Unknown potential mesh topology.');
  const usableSeed=meshSeed&&JSON.stringify(meshSeed.contours)===JSON.stringify(contours);
  const deformed=type==='triangular'&&usableSeed?deformTriangularWakes(meshSeed,wakePaths):null;
  const original=type==='triangular'?(deformed??triangularAirfoilMesh(contours,{...controls,wakePaths}))
    :multielementMesh(contours,{rows:3,padding:4,outerGrowth:1.5,wallStretch:4,...controls,wakePaths});
  const wakeIndices=original.cuts.filter(c=>c.type==='wake-cut').map(c=>c.face);
  const mesh=partitionPotentialInterfaces(original,wakeIndices),count=mesh.cells.length,elements=outer.bodies.length;
  onMesh?.(mesh);
  const wallIndices=mesh.faces.map((f,i)=>f.boundary?.element===undefined?-1:i).filter(i=>i>=0);
  const circulation=airfoilCirculation(mesh,{alpha:input.alpha??0});
  const matching=farfield==='multipole'?createMultipoleFarfield(mesh,{alpha:input.alpha??0,mach,circulation,offset:count+elements}):null;
  const wallOffset=count+elements+(matching?.size??0),wallColumns=new Map(wallIndices.map((i,k)=>[i,wallOffset+k]));
  const wakeOffset=wallOffset+wallIndices.length,blOffset=wakeOffset+mesh.potentialInterfaces.length;
  const sources=meshDisplacementSources(original,outer,{offset:blOffset,thicknessScale:bl.thicknessScale,mach,bodyMassInterpolation});
  const panelFaces=new Map();
  for(const i of wallIndices){const p=sources.bodyPanels.get(i),f=mesh.faces[i];
    if(!panelFaces.has(p))panelFaces.set(p,[]);
    panelFaces.get(p).push({...f,index:i,a:mesh.vertices[f.a],b:mesh.vertices[f.b]});
  }
  const faceIndices=new Map(mesh.faces.map((f,i)=>[f,i]));
  let frozen=false,displacementScale=1,cached,edgeForms;
  const scaledSource=source=>({...source,evaluate:x=>{
    const result=source.evaluate(x);
    return{value:displacementScale*result.value,
      derivatives:new Map([...result.derivatives].map(([i,v])=>[i,displacementScale*v]))};
  }});
  const invalidate=()=>{cached=undefined;};
  const local=x=>x.subarray(blOffset);
  const blocks=(x,jacobian)=>{
    const state=local(x);
    if(!cached||cached.x.some((v,i)=>v!==state[i]))cached={x:state.slice(),residual:bl.residual(state)};
    if(jacobian&&!cached.jacobian)cached.jacobian=bl.jacobian(state);
    return cached;
  };
  const velocityForms=(reconstruction,faceVelocities)=>{
    const edges=new Array(outer.total),same=(p,q)=>Math.hypot(p.x-q.x,p.y-q.y)<1e-10;
    const atVertex=(point,indices,tangent)=>{
      if(!indices.length)throw new Error('BL station has no adjacent potential cells.');
      const terms=[];
      for(const i of indices){const f=mesh.faces[i],[u,v]=reconstruction.velocity(f.owner,point);terms.push([u,tangent.tx/indices.length],[v,tangent.ty/indices.length]);}
      return sumAffine(terms);
    };
    for(const [e,b] of outer.bodies.entries())for(let id=b.start;id<=b.end;id++){
      const point=b.points[id-b.start],left=outer.panels[Math.max(b.first,b.first+id-b.start-1)],right=outer.panels[Math.min(b.last,b.first+id-b.start)];
      let tx=left.tx+right.tx,ty=left.ty+right.ty;const length=Math.hypot(tx,ty);tx/=length;ty/=length;
      const indices=wallIndices.filter(i=>mesh.faces[i].boundary.element===e&&[mesh.vertices[mesh.faces[i].a],mesh.vertices[mesh.faces[i].b]].some(p=>same(p,point))
        &&(id!==b.start||sources.bodyPanels.get(i).node===b.start)&&(id!==b.end||sources.bodyPanels.get(i).node===b.end-1));
      edges[id]=atVertex(point,indices,{tx,ty});
    }
    if(edgeVelocitySampling!=='vertex')for(const b of outer.bodies){
      const samples=[];
      for(let j=b.first;j<=b.last;j++){
        const p=outer.panels[j],terms=[];let covered=0;
        if(edgeVelocitySampling==='station-average'){
          // Half-panel samples place the recovery interval around each BL
          // station. Whole-panel means followed by nodal interpolation can
          // annihilate alternating station-scale velocity/displacement modes.
          for(let half=0;half<2;half++)samples.push({s:b.s[p.node-b.start]+(.25+.5*half)*p.length,
            velocity:integratedPanelVelocity(p,panelFaces.get(p),(f,point)=>reconstruction.velocity(f.owner,point),{lo:.5*half,hi:.5*(half+1)})});
          continue;
        }
        for(const {index:i} of panelFaces.get(p)){
          const weight=mesh.faces[i].length/p.length;covered+=weight;
          terms.push([faceVelocities[i][0],weight*p.tx],[faceVelocities[i][1],weight*p.ty]);
        }
        if(Math.abs(covered-1)>1e-8)throw new Error('Incomplete BL panel velocity coverage.');
        samples.push({s:b.s[p.node-b.start]+.5*p.length,velocity:sumAffine(terms)});
      }
      const velocities=interpolatePanelAverages(samples,b.s);
      velocities.forEach((v,j)=>edges[b.start+j]=sumAffine([[edges[b.start+j],1-edgeVelocityFraction],[v,edgeVelocityFraction]]));
    }
    outer.wakes.forEach((w,e)=>{
      edges[w.start]=edges[w.body.end];
      const interfaces=mesh.potentialInterfaces.filter(p=>original.cuts.find(c=>c.face===p.original)?.element===e);
      for(let j=1;j<w.points.length;j++){
        const point=w.points[j],indices=interfaces.flatMap(p=>p.faces).filter(i=>[mesh.vertices[mesh.faces[i].a],mesh.vertices[mesh.faces[i].b]].some(p=>same(p,point)));
        edges[w.start+j]=atVertex(point,indices,w.tangents[j]);
      }
    });
    return edges;
  };
  const blColumns=id=>{
    const cols=new Set(),add=i=>{for(let k=0;k<4;k++)cols.add(blOffset+4*i+k);};
    const body=outer.bodies.find(b=>id>=b.start&&id<=b.end),wake=outer.wakes.find(w=>id>=w.start&&id<=w.end),b=body??wake.body;
    if(body)for(let i=Math.max(b.start,id-1);i<=Math.min(b.end,id+1);i++)add(i);
    else{add(id);if(id===wake.start){add(b.start);add(b.end);}else add(id-1);}
    for(let i=b.start;i<=b.end;i++)cols.add(blOffset+4*i+3);
    return cols;
  };
  const systemOptions={mach,sparse,order,fluxQuadrature,baseVelocity:circulation.baseVelocity,baseFluxIntegral:circulation.baseFluxIntegral,
    boundaryPotential:f=>f.boundary.type==='potential-interface'?affine(0,[[wakeOffset+f.boundary.interface,1]]):(matching?.potential(f)??affine()),
    boundaryNormalVelocity:f=>f.boundary.element===undefined?affine():affine(0,[[wallColumns.get(faceIndices.get(f)),1]]),
    constraints:({faces,reconstruction,fluxColumns,meanVelocities})=>{
      const edges=velocityForms(reconstruction,meanVelocities);if(!edgeForms)edgeForms=edges;
      const constraints=outer.bodies.map(b=>sumAffine([[edges[b.start],1],[edges[b.end],1]]));
      if(matching)constraints.push(...matching.constraints({faces}));
      for(const i of wallIndices)constraints.push(potentialFluxConstraint(mesh,fluxColumns,[i],{source:scaledSource(sources.bodySources.get(i)),scale:mesh.faces[i].length}));
      for(const p of mesh.potentialInterfaces)constraints.push(potentialFluxConstraint(mesh,fluxColumns,p.faces,{source:scaledSource(sources.wakeSources.get(p.original)),scale:mesh.faces[p.faces[0]].length}));
      for(let row=0;row<bl.initial.length;row++){
        const id=Math.floor(row/4),columns=blColumns(id),velocity=row%4===3?sumAffine([[affine(0,[[blOffset+row,1]]),1],[edges[id],-1]]):null;
        if(velocity)for(const col of velocity.coefficients.keys())columns.add(col);
        constraints.push({columns,evaluate(x,{jacobian}){
          if(frozen)return{value:x[blOffset+row]-bl.initial[row],derivatives:new Map([[blOffset+row,1]])};
          if(velocity){let value=velocity.constant;for(const [col,v] of velocity.coefficients)value+=v*x[col];return{value,derivatives:velocity.coefficients};}
          const b=blocks(x,jacobian),derivatives=new Map();
          if(jacobian)for(const col of columns){const value=b.jacobian[row*bl.initial.length+col-blOffset];if(value!==0)derivatives.set(col,value);}
          return{value:b.residual[row],derivatives};
        }});
      }
      return constraints;
    }};
  const system=createPotentialSystem(mesh,systemOptions);
  const initial=new Float64Array(system.n);initial.set(bl.initial,blOffset);
  for(const [i,col] of wallColumns)initial[col]=-sources.bodySources.get(i).evaluate(initial).value/mesh.faces[i].length;
  const updateActive=x=>{const changed=bl.updateActive(local(x));invalidate();return changed;};
  const restore=snapshot=>{bl.restore(snapshot);invalidate();};
  return{...system,bl,outer,initial,blOffset,wallOffset,wakeOffset,wallIndices,sources,edgeForms,originalMesh:original,meshDeformed:!!deformed,alpha:input.alpha??0,farfield:matching,bodyMassInterpolation,edgeVelocitySampling,edgeVelocityFraction,
    setFrozen:value=>{frozen=value;invalidate();},updateActive,snapshot:bl.snapshot,restore,
    setDisplacementScale:value=>{
      if(!frozen||!Number.isFinite(value)||value<0||value>1)throw new Error('Displacement continuation is only an outer-flow initialization.');
      displacementScale=value;
    },
    incompressibleInitialization:()=>{
      if(!frozen)throw new Error('An incompressible guess is only an outer-flow initialization.');
      return system.incompressibleInitialization();
    },
    admissible:x=>{
      if(!bl.admissible(local(x))||!system.admissible(x))return false;
      try{for(let i=0;i<outer.total;i++)isentropicState(x[blOffset+4*i+3],0,{mach});return true;}catch{return false;}
    },
    releaseNonzeroStagnation:x=>{const changed=bl.releaseNonzeroStagnation(local(x));invalidate();return changed;},
    guardedTrialState:(x,delta,step)=>{
      const candidate=x.map((v,i)=>v+step*delta[i]);
      candidate.set(bl.guardedTrialState(local(x),local(delta),step),blOffset);
      return candidate;
    },
    limitedVariables:()=>({offset:blOffset,length:bl.initial.length})};
}

export function solveCoupledPotential(input,controls={}){
  const target=input.mach??0;
  const useContinuation=controls.machContinuation??(input.elements.length>=3&&!input.seed);
  if(!useContinuation||target===0||controls.initial||controls.previous)return solveCoupledPotentialSampling(input,controls);
  const history=[],report=mach=>h=>controls.onIteration?.({...h,continuationMach:mach});
  let current=solveCoupledPotentialSampling({...input,mach:0,seed:input.seed?remapSeedMach(input.seed,0):undefined},
    {...controls,onIteration:report(0)});
  history.push(...current.history.map(h=>({...h,continuationMach:0})));
  if(!current.converged)return{...current,history,stage:'mach-initialization',reason:`Mach-zero coupled initialization failed: ${current.reason}`};
  let mach=0,step=Math.min(.05,target);
  while(mach<target){
    const next=Math.min(target,mach+step);let attempt;
    try{
      const seed=remapSeedMach(current.system.bl.exportSeed(current.x.subarray(current.system.blOffset)),next);
      const initial=current.x.slice();initial.set(seed.x,current.system.blOffset);
      attempt=solveCoupledPotentialSampling({...input,mach:next,seed,wakePaths:current.system.outer.wakes.map(w=>w.points)},
        {...controls,initial,previous:current,onIteration:report(next)});
    }catch(error){
      if(!/Local sonic flow|stagnation enthalpy|Invalid Newton initial state/.test(error.message))throw error;
      attempt={...current,converged:false,history:[],reason:error.message};
    }
    history.push(...attempt.history.map(h=>({...h,continuationMach:next})));
    if(attempt.converged){current=attempt;mach=next;step=Math.min(.05,step*1.5);}
    else if(step>1e-4)step/=2;
    else return{...attempt,converged:false,history,stage:'mach-initialization',reason:`Mach continuation stopped after M=${mach}: ${attempt.reason}`};
  }
  return{...current,history};
}

function solveCoupledPotentialSampling(input,controls){
  const averaged=['panel-average','station-average'].includes(controls.edgeVelocitySampling);
  const enabled=controls.edgeVelocityContinuation??(input.elements.length>=3&&averaged);
  if(!enabled||controls.initial||controls.previous||!averaged)return solveCoupledPotentialDirect(input,{...controls,edgeVelocityFraction:1});
  // Continue the velocity closure only to obtain a starting point. Every
  // intermediate stage solves all BL and outer variables simultaneously.
  // A result at any fraction below one is never returned as converged.
  const history=[],report=fraction=>h=>controls.onIteration?.({...h,edgeVelocityFraction:fraction});
  let current=solveCoupledPotentialDirect(input,{...controls,edgeVelocityFraction:0,onIteration:report(0)});
  history.push(...current.history.map(h=>({...h,edgeVelocityFraction:0})));
  if(!current.converged)return{...current,history,stage:'sampling-initialization',reason:`Vertex-velocity initialization failed: ${current.reason}`};
  let fraction=0,step=.25;
  while(fraction<1){
    const next=Math.min(1,fraction+step),seed=current.system.bl.exportSeed(current.x.subarray(current.system.blOffset));
    let attempt;
    try{attempt=solveCoupledPotentialDirect({...input,seed,wakePaths:current.system.outer.wakes.map(w=>w.points)},
      {...controls,edgeVelocityFraction:next,initial:current.x.slice(),previous:current,onIteration:report(next)});}
    catch(error){
      if(!/Local sonic flow|stagnation enthalpy|Invalid Newton initial state/.test(error.message))throw error;
      attempt={...current,converged:false,history:[],reason:error.message};
    }
    history.push(...attempt.history.map(h=>({...h,edgeVelocityFraction:next})));
    if(attempt.converged){current=attempt;fraction=next;step=Math.min(.25,step*1.5);}
    else if(step>1e-4)step/=2;
    else return{...attempt,converged:false,history,stage:'sampling-initialization',reason:`Edge-velocity continuation stopped after ${fraction}: ${attempt.reason}`};
  }
  return{...current,history};
}

function solveCoupledPotentialDirect(input,{maxIterations=60,tolerance=1e-8,onIteration,onState,initial,previous,linearTolerance=1e-7,fillLevel=1,linearBackend='klu',eventStepFraction=.25,...options}={}){
  const system=createCoupledPotential(input,{...options,meshSeed:previous?.system.originalMesh??options.meshSeed});
  const sameTopology=previous?.converged&&previous.system.n===system.n&&previous.system.blOffset===system.blOffset
    &&previous.system.mesh.faces.length===system.mesh.faces.length&&system.mesh.faces.every((f,i)=>{
      const p=previous.system.mesh.faces[i];return f.a===p.a&&f.b===p.b&&f.owner===p.owner&&f.neighbor===p.neighbor
        &&f.boundary?.type===p.boundary?.type&&f.boundary?.element===p.boundary?.element&&f.boundary?.interface===p.boundary?.interface;
    });
  const warm=initial??(sameTopology?previous.x.slice():null);
  let x=warm??system.initial,outerInitialization;
  if(!warm){
    system.setFrozen(true);
    const solveStart=(initial,scale)=>solvePotential(system,{initial,tolerance:Math.min(tolerance,1e-10),linearTolerance,fillLevel,linearBackend,
      onIteration:h=>onIteration?.({...h,stage:'outer-initialization',displacementScale:scale})});
    let start;
    try{start=solveStart(x,1);}catch(error){start={converged:false,reason:error.message,x,history:[]};}
    if(!start.converged){
      // Large isolated-BL source gradients can make the very first outer
      // state inadmissible. Initialize circulation/flow without displacement,
      // then continue its strength to one. The joint solve always uses one.
      system.setDisplacementScale(0);
      const zero=x.slice();zero.fill(0,0,system.blOffset);
      try{start=solveStart(zero,0);}catch(error){start={converged:false,reason:error.message,x:zero,history:[]};}
      if(!start.converged){
        // Zero potential is not a uniform velocity near a solid wall: the
        // Neumann reconstruction already enforces impermeability. Its poor
        // first guess can be locally sonic even for a subcritical solution.
        const linear=solvePotential(system.incompressibleInitialization(),{initial:zero,linearTolerance,fillLevel,linearBackend,
          onIteration:h=>onIteration?.({...h,stage:'incompressible-initialization'})});
        if(linear.converged){
          try{start=solveStart(linear.x,0);}catch(error){start={...linear,converged:false,reason:error.message};}
        }else start=linear;
      }
      let scale=0,step=.25;
      while(start.converged&&scale<1){
        const next=Math.min(1,scale+step);system.setDisplacementScale(next);
        let attempt;
        try{attempt=solveStart(start.x,next);}catch(error){attempt={converged:false,reason:error.message};}
        if(attempt.converged){start=attempt;scale=next;step=Math.min(.25,step*1.5);}
        else if(step>1/1024)step/=2;
        else{start={...start,converged:false,reason:`Outer displacement initialization stopped at ${scale}: ${attempt.reason}`};break;}
      }
      system.setDisplacementScale(1);
    }
    system.setFrozen(false);
    outerInitialization={converged:start.converged,reason:start.reason};
    if(!start.converged){
      // A frozen BL is only an initializer. Its prescribed mass gradient can
      // be incompatible with a subcritical outer root even when the joint
      // equations admit a solution after the BL changes. Use its last finite
      // state as a guess, with displacement restored to one and all BL rows
      // released. Only the complete joint residual can certify a result.
      let usable=false;
      try{usable=system.admissible(start.x)&&system.evaluate(start.x).residual.every(Number.isFinite);}catch{}
      if(!usable)return{...start,system,stage:'outer-initialization',outerInitialization};
    }
    x=start.x;
    // Preserve a supplied BL seed as a complete state. Its velocity and
    // thickness already satisfy the same BL equations; replacing velocity
    // alone after remeshing can destroy that consistency. Newton adjusts
    // all variables together to match the new outer discretization.
    // The isolated native edge velocities can be almost zero at a node that
    // is far from stagnation in the assembly. Start the BL on the actual
    // initialized outer field, retaining its finite thickness/shear guesses.
    // This changes only the initial iterate; displacement is recomputed by
    // the unmodified simultaneous residual below.
    if(!input.seed){
      const candidate=x.slice(),saved=system.snapshot();
      for(let i=0;i<system.outer.total;i++)candidate[system.blOffset+4*i+3]=evaluateAffine(system.edgeForms[i],x);
      try{
        system.releaseNonzeroStagnation(candidate);system.updateActive(candidate);
        if(system.admissible(candidate)&&system.evaluate(candidate).residual.every(Number.isFinite))x=candidate;
        else system.restore(saved);
      }catch{system.restore(saved);}
    }
  }
  let lastState,lastResidual;
  const equations={tolerance,eventStepFraction,onIteration:h=>{
      let worstRow=0;for(let i=1;i<lastResidual.length;i++)if(Math.abs(lastResidual[i])>Math.abs(lastResidual[worstRow]))worstRow=i;
      const diagnostics={...h,worstRow,transitions:system.bl.surfaces.map(s=>s.transitionId)};
      onIteration?.(diagnostics);onState?.({system,x:lastState,iteration:diagnostics});
    },
    residual:x=>{lastState=x;lastResidual=system.evaluate(x).residual;return lastResidual;},
    jacobian:x=>system.evaluate(x,{jacobian:true}).jacobian,
    linearSolve:system.sparse?(a,b)=>solvePotentialLinear(a,b,{linearTolerance,fillLevel,linearBackend}):solveLinear,
    admissible:system.admissible,limitedVariables:system.limitedVariables,updateActive:system.updateActive,snapshot:system.snapshot,restore:system.restore};
  const saved=system.snapshot();
  let result=solveActiveNewton({...equations,initial:x,maxIterations});
  const remaining=maxIterations-result.history.filter(h=>h.iteration>0).length;
  if(!result.converged&&remaining>=3){
    system.restore(saved);
    const retry=solveActiveNewton({...equations,initial:x,maxIterations:remaining,trialState:system.guardedTrialState,
      onIteration:h=>equations.onIteration({...h,attempt:2})});
    result={...retry,history:[...result.history,...retry.history.map(h=>({...h,attempt:2}))]};
  }
  while(result.converged&&system.releaseNonzeroStagnation(result.x)){
    const remaining=maxIterations-result.history.filter(h=>h.iteration>0).length;
    const next=solveActiveNewton({...equations,initial:result.x,maxIterations:Math.max(0,remaining)});
    result={...next,history:[...result.history,...next.history]};
  }
  return{...result,...system.evaluate(result.x),system,stage:'coupled',outerInitialization,boundaryLayer:system.bl.decode(result.x.subarray(system.blOffset))};
}

export function potentialWakeVelocity(result,element,point){
  const {system,x}=result,{mesh,reconstruction}=system;
  const candidates=mesh.potentialInterfaces.filter(p=>{
    const cut=mesh.cuts.find(c=>c.face===p.original);if(cut.element!==element)return false;
    const f=mesh.faces[p.faces[0]],a=mesh.vertices[f.a],b=mesh.vertices[f.b];
    return point.x>=Math.min(a.x,b.x)-1e-12&&point.x<=Math.max(a.x,b.x)+1e-12;
  });
  if(!candidates.length)throw new Error('Wake velocity point is outside the mesh interface.');
  const velocity=[0,0];
  for(const p of candidates)for(const i of p.faces){
    const f=mesh.faces[i],atPoint=reconstruction.velocity(f.owner,point),atFace=reconstruction.velocity(f.owner,f);
    for(let k=0;k<2;k++)velocity[k]+=evaluateAffine(sumAffine([[system.faces[i][k],1],[atPoint[k],1],[atFace[k],-1]]),x)/(2*candidates.length);
  }
  return{u:velocity[0],v:velocity[1]};
}

// A separate geometry iteration closes each wake's mean normal velocity.
// Every intermediate shape still solves the full simultaneous potential/BL
// system. Success requires both the equation and wake-shape tolerances.
export function solveMultielementPotentialViscous(input,controls={}){
  const initialize=controls.assemblyInitialization??input.elements.length>=3;
  if(!initialize||input.seed||controls.initial)return solvePotentialWakeGeometry(input,controls);
  // A closed incompressible assembly supplies a consistent wake and BL guess.
  // Isolated airfoil guesses can put a slat wake in a region where the fixed
  // geometry finite-Mach equations have no accessible root. Every final
  // potential/BL equation and wake normal-velocity condition is still solved.
  const start=solveMultielementViscous({...input,mach:0},{...controls,wakeRelaxation:.8,
    onIteration:h=>controls.onIteration?.({...h,stage:'assembly-initialization',continuationMach:0})});
  let prepared=input,seedFailure;
  if(start.converged){
    try{prepared={...input,seed:remapSeedMach(start.system.exportSeed(start.x),input.mach??0),wakePaths:start.system.outer.wakes.map(w=>w.points)};}
    catch(error){if(!/Local sonic flow|stagnation enthalpy/.test(error.message))throw error;seedFailure=error.message;}
  }
  const result=solvePotentialWakeGeometry(prepared,controls);
  return{...result,history:[...start.history.map(h=>({...h,stage:'assembly-initialization',continuationMach:0})),...result.history],
    assemblyInitialization:{converged:start.converged,used:prepared!==input,reason:seedFailure??start.reason,
      equationResidual:start.history.at(-1)?.residual,wakeResidual:start.wakeResidual}};
}

function solvePotentialWakeGeometry(input,{wakeTolerance=1e-6,maxWakeIterations=40,wakeRelaxation=.5,onIteration,initial,...options}={}){
  if(!(wakeTolerance>0)||!Number.isFinite(wakeTolerance)||!Number.isInteger(maxWakeIterations)||maxWakeIterations<0||!(wakeRelaxation>0)||wakeRelaxation>1)throw new Error('Invalid potential wake controls.');
  let wakePaths=input.wakePaths,seed=input.seed,lengths,result,previous;const wakeHistory=[],history=[];
  for(let iteration=0;iteration<=maxWakeIterations;iteration++){
    result=solveCoupledPotential({initialization:input.elements.length>1?'auto':'native',wakeInitialization:'inviscid',...input,wakePaths,seed},
      {...options,previous,initial:iteration===0?initial:undefined,onIteration:h=>onIteration?.({...h,wakeIteration:iteration})});
    history.push(...result.history.map(h=>({...h,wakeIteration:iteration})));
    if(!result.converged&&result.system.meshDeformed){
      // Moving a constrained wake can leave an admissible but poor triangle
      // stencil. Rebuild at the same physical wake geometry and solve all
      // equations again. The failed state never supplies published forces.
      onIteration?.({kind:'remesh',wakeIteration:iteration,reason:result.reason});
      result=solveCoupledPotential({...input,wakePaths:result.system.outer.wakes.map(w=>w.points),seed:result.system.bl.exportSeed(result.x.subarray(result.system.blOffset))},
        {...options,meshSeed:undefined,onIteration:h=>onIteration?.({...h,wakeIteration:iteration,remeshed:true})});
      history.push(...result.history.map(h=>({...h,wakeIteration:iteration,remeshed:true})));
    }
    if(!result.converged)return{...result,history,wakeHistory,wakeConverged:false,wakeResidual:Infinity};
    const {outer,bl,blOffset}=result.system;
    lengths??=outer.wakes.map(w=>w.segments.map(s=>s.length));
    let residual=0;
    const directions=outer.wakes.map((w,e)=>w.segments.map(p=>{
      const {u,v}=potentialWakeVelocity(result,e,p),speed=Math.hypot(u,v);
      if(u*p.tx+v*p.ty<=0)throw new Error('Potential wake flow reverses.');
      residual=Math.max(residual,Math.abs(u*p.nx+v*p.ny));
      return{tx:u/speed,ty:v/speed};
    }));
    wakeHistory.push({iteration,residual,equationResidual:result.diagnostics.residual});
    onIteration?.({kind:'wake',wakeIteration:iteration,wakeResidual:residual});
    if(residual<=wakeTolerance)return{...result,history,wakeHistory,wakeConverged:true,wakeResidual:residual};
    if(iteration===maxWakeIterations)return{...result,converged:false,reason:'wake iteration limit',history,wakeHistory,wakeConverged:false,wakeResidual:residual};
    wakePaths=outer.wakes.map((w,e)=>{
      const points=[{...w.points[0]}];
      for(let j=0;j<w.segments.length;j++){
        const p=w.segments[j],d=directions[e][j],tx=(1-wakeRelaxation)*p.tx+wakeRelaxation*d.tx,ty=(1-wakeRelaxation)*p.ty+wakeRelaxation*d.ty,scale=lengths[e][j]/Math.hypot(tx,ty),a=points.at(-1);
        points.push({x:a.x+scale*tx,y:a.y+scale*ty});
      }
      return points;
    });
    seed=bl.exportSeed(result.x.subarray(blOffset));
    previous=result;
  }
}
