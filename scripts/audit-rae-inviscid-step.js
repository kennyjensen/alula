// Exact next-iteration reproduction; no grid regeneration or BL equations.
import fs from 'node:fs';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {solveStreamtubeIses, continueEulerDissipationFromCheckpoint} from '../src/euler/streamtube-ises-update.js';
export {continueEulerDissipationFromCheckpoint};
import {createStreamtubeBodySystem} from '../src/euler/streamtube-body.js';
import {adjustStreamtubeInlets,reparameterizeStreamtubeInlets,dekinkStreamtubeInteriors} from '../src/geometry/streamtube-grid-maintenance.js';
import {assembleTangentialCoordinate,redistributeStreamtubeTangentially} from '../src/geometry/streamtube-tangential-redistribution.js';
import {streamtubeMeshSnapshot} from '../src/euler/streamtube-mesh-preview.js';
import {proposeDensityNewton} from '../src/euler/streamtube-density-newton.js';
import {solveSparseDirect} from '../src/numerics/klu.js';
import {streamtubeGridConvexity,limitStreamtubeGridStep} from '../src/geometry/streamtube-convex-step.js';
export const directory='docs/solver-reliability/rae-inviscid-audit';
export function reproduceRaeInviscidStep({traceUpdates=true}={}) {
 const path=`${directory}/fixtures/rae32-before-step5.json`, bytes=fs.readFileSync(path);
 const cp=JSON.parse(bytes),start=performance.now();
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:1,tolerance:1e-10,traceUpdates});
 const report={schemaVersion:1,fixtureSHA256:crypto.createHash('sha256').update(bytes).digest('hex'),
  sourceSHA256:crypto.createHash('sha256').update(fs.readFileSync('src/euler/streamtube-ises-update.js')).digest('hex'),
  acceptedColdIteration:4,attemptedColdIteration:5,seconds:(performance.now()-start)/1000,
  converged:result.converged,reason:result.reason,diagnostics:result.diagnostics,
  history:result.history,lastRejectedStep:result.lastRejectedStep,linearDiagnostics:result.linearDiagnostics,updateTrace:result.updateTrace};
 return {cp,result,report};
}

// Read-only discriminator for a finite grid-maintenance jump. Compare the
// centered SMOVE movement with inversion of the *same* solved coordinate
// labels along the existing streamline polylines. Neither candidate is an
// accepted Newton state, and neither changes the production movement policy.
export function auditEulerTangentialMovement(source, {pairs=[5,10,20,50],correctionScale=1}={}) {
 if(!Number.isFinite(correctionScale)||correctionScale<=0||correctionScale>1)
  throw new Error('Invalid tangential audit correction scale.');
 const before=JSON.stringify(source),cp=structuredClone(source);
 const system=createStreamtubeBodySystem(cp.input);
 const x=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes);
 const value=system.evaluate(x),rows=[];
 function measure(stage,nodes,details={}) {
  let maximumNodeMotion=0;
  nodes.forEach((group,g)=>group.forEach((row,i)=>row.forEach((p,j)=>{
   const q=value.nodes[g][i][j];maximumNodeMotion=Math.max(maximumNodeMotion,Math.hypot(p.x-q.x,p.y-q.y));
  })));
  const row={stage,maximumNodeMotion,quality:streamtubeMeshSnapshot({system,nodes}).quality,...details};
  try {
   const target=createStreamtubeBodySystem(cp.input),v=target.evaluate(target.adoptGeometry(x,nodes));
   row.maximumResidual=v.diagnostics.residual;
   row.squaredResidual=v.residual.reduce((a,r)=>a+r*r,0);
  }catch(error){row.error={message:error.message,code:error.code,diagnostics:error.diagnostics};}
  rows.push(row);
 }
 measure('base',value.nodes);
 const inlet=adjustStreamtubeInlets(value.nodes,system.layout.bodies,cp.continuation.fractions,{preserveConvexity:true});
 const {nodes:inletNodes,...inletDetails}=inlet;
 measure('inlet',inletNodes,{inlet:inletDetails});
 const dekink=dekinkStreamtubeInteriors(inlet.nodes,{preserveConvexity:true});
 measure('dekink',dekink.nodes,{repairs:dekink.repairs});
 for(const count of pairs) {
  try {
   const inverse=[],coordinates=[];
   const passages=dekink.nodes.map((nodes,g)=>{
    const options={referenceBank:g===0?nodes[0].length-1:0,
     fixedBanks:[g!==0,g!==dekink.nodes.length-1],pairs:count,correctionScale,quadratureDomain:'sampled-positive'};
    return {nodes,options,ordinary:redistributeStreamtubeTangentially(nodes,options)};
   });
   // Preserve the ordinary result even if the inverse field has no monotone
   // bracket. Failure of the alternative must not hide the control evidence.
   measure('centered',passages.map(p=>p.ordinary.nodes),{pairs:count,coordinates});
   for(let g=0;g<passages.length;g++) {
    const {nodes,options,ordinary}=passages[g];
    const c=assembleTangentialCoordinate(nodes,options),moved=nodes.map(row=>row.map(p=>({...p})));
    let minimumLabelSpan=Infinity;
    for(let j=c.firstStreamline;j<=c.lastStreamline;j++) {
     const labels=Array.from(c.baseline,(v,i)=>v+(i>0&&i<c.nx
      ?correctionScale*ordinary.solution.x[(i-1)*(c.lineDimensions.nt-1)+j-c.firstStreamline]:0));
     for(let i=0;i<c.nx;i++) {
      const span=labels[i+1]-labels[i];minimumLabelSpan=Math.min(minimumLabelSpan,span);
      if(!(span>0))throw new Error(`Nonmonotone tangential coordinate in passage ${g}, station ${i}, streamline ${j}.`);
     }
     let k=0;
     for(let i=1;i<c.nx;i++) {
      const label=c.baseline[i];while(k<c.nx-1&&labels[k+1]<label)k++;
      const fraction=(label-labels[k])/(labels[k+1]-labels[k]);
      if(!(fraction>=0&&fraction<=1))throw new Error('Unbracketed tangential coordinate.');
      const a=nodes[k][j],b=nodes[k+1][j];
      moved[i][j]={x:a.x+fraction*(b.x-a.x),y:a.y+fraction*(b.y-a.y)};
     }
    }
    inverse.push(moved);
    coordinates.push({group:g,relativeResidual:ordinary.solution.relativeResidual,minimumLabelSpan});
   }
   measure('inverse',inverse,{pairs:count,coordinates});
  }catch(error){rows.push({stage:'coordinate-movement',pairs:count,error:{message:error.message,code:error.code}});}
 }
 if(JSON.stringify(source)!==before)throw new Error('Grid audit modified its source checkpoint.');
 return {sourceUnchanged:true,correctionScale,rows,scope:'Zero Newton step, unchanged flow unknowns; grid and residual diagnosis only.'};
}

// Isolate the ordering of the raw-grid gate and inlet re-spacing. The
// proposed polyline may be nonconvex: it is only a geometric intermediate.
// Check the repaired endpoint AND its complete path from the accepted grid;
// evaluate gas/residuals only for convex endpoints. Nothing is committed.
export function auditEulerInletTrialRepair(source, {halvings=16,onTrial}={}) {
 if(!Number.isInteger(halvings)||halvings<0||halvings>30)throw new Error('Invalid inlet trial audit controls.');
 const before=JSON.stringify(source),cp=structuredClone(source),system=createStreamtubeBodySystem(cp.input);
 const x=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes),value=system.evaluate(x);
 const squared=residual=>residual.reduce((sum,r)=>sum+r*r,0),baseSquaredResidual=squared(value.residual);
 const linear=solveSparseDirect(system.jacobian(x,{sparse:true}),value.residual.map(r=>-r));
 function measure(state,nodes) {
  const quality=streamtubeGridConvexity(nodes);
  const row={quality:{...quality,invalidCells:quality.invalidCells.slice(0,4)},invalidCount:quality.invalidCells.length};
  if(!quality.valid)return row;
  try {
   const candidate=createStreamtubeBodySystem(cp.input),v=candidate.evaluate(candidate.adoptGeometry(state,nodes));
   row.squaredResidual=squared(v.residual);row.maximumResidual=v.diagnostics.residual;
   row.residualByFamily=v.diagnostics.residualByFamily;
  }catch(error){row.error={message:error.message,code:error.code,diagnostics:error.diagnostics};}
  return row;
 }
 const rows=[];
 for(let i=0;i<=halvings;i++) {
  const proposal=proposeDensityNewton(system,x,linear.x,{stagnationLimiter:cp.continuation.stagnationLimiter,maximumStep:2**-i});
  const nodes=system.decode(proposal.x).nodes,repaired=reparameterizeStreamtubeInlets(nodes,cp.input.bodies,cp.continuation.fractions).nodes;
  const row={maximumStep:2**-i,step:proposal.step,raw:measure(proposal.x,nodes),inlet:measure(proposal.x,repaired)};
  try {
   const limit=limitStreamtubeGridStep(value.nodes,repaired);
   row.completePath={limited:limit.limited,step:limit.step,limiter:limit.limiter};
  }catch(error){row.pathError={message:error.message,code:error.code};}
  rows.push(row);onTrial?.(structuredClone(row));
 }
 if(JSON.stringify(source)!==before)throw new Error('Inlet trial audit modified its source checkpoint.');
 return {sourceUnchanged:true,baseSquaredResidual,linearRelativeResidual:linear.relativeResidual,rows,
  scope:'One frozen Newton direction; raw and inlet-repaired trials, no accepted updates or changed equations.'};
}
if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){
 if(process.argv[2]==='--inlet-trial-repair'){
  const [source,output]=process.argv.slice(3);
  if(!source||!output)throw new Error('Usage: --inlet-trial-repair checkpoint.json output.json');
  const bytes=fs.readFileSync(source),report=auditEulerInletTrialRepair(JSON.parse(bytes),{
   onTrial:r=>console.log(JSON.stringify({step:r.step,rawInvalid:r.raw.invalidCount,inletInvalid:r.inlet.invalidCount,
    squaredResidual:r.inlet.squaredResidual,completePath:r.completePath,error:r.inlet.error?.message})),
  });
  report.source=source;report.sourceSHA256=crypto.createHash('sha256').update(bytes).digest('hex');
  fs.writeFileSync(output,JSON.stringify(report,null,2));
 }else if(process.argv[2]==='--grid-maintenance'){
  const [source,output]=process.argv.slice(3);
  if(!source||!output)throw new Error('Usage: --grid-maintenance checkpoint.json output.json');
  const bytes=fs.readFileSync(source),report=auditEulerTangentialMovement(JSON.parse(bytes));
  report.source=source;report.sourceSHA256=crypto.createHash('sha256').update(bytes).digest('hex');
  fs.writeFileSync(output,JSON.stringify(report,null,2));
  for(const r of report.rows)console.log(JSON.stringify({stage:r.stage,pairs:r.pairs,valid:r.quality?.valid,
   residual:r.maximumResidual,squaredResidual:r.squaredResidual,error:r.error?.message}));
 }else if(process.argv[2]==='--cold-dissipation'){
  const [source,output]=process.argv.slice(3);
  if(!source||!output)throw new Error('Usage: --cold-dissipation checkpoint.json output.json');
  const result=continueEulerDissipationFromCheckpoint(JSON.parse(fs.readFileSync(source)),{
   onIteration:h=>console.log(JSON.stringify(h)),
   onCheckpoint:cp=>fs.writeFileSync(`${output}.checkpoint.json`,JSON.stringify(cp)),
  });
  const {checkpoint,...report}=result;
  fs.writeFileSync(output,JSON.stringify(report,null,2));
  console.log(JSON.stringify({converged:result.converged,reason:result.reason,seconds:result.seconds}));
 }else{
  const {report}=reproduceRaeInviscidStep();
  fs.writeFileSync(`${directory}/exact-step-baseline.json`,JSON.stringify(report,null,2));
  console.log(JSON.stringify({converged:report.converged,reason:report.reason,trials:report.updateTrace?.length,seconds:report.seconds}));
 }
}
