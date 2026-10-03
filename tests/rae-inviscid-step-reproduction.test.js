// Exact next-iteration regression. This does not claim that the full RAE case
// is converged; it protects the recovered residual-decreasing Newton update.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {reproduceRaeInviscidStep,auditEulerTangentialMovement} from '../scripts/audit-rae-inviscid-step.js';
import {createStreamtubeBodySystem} from '../src/euler/streamtube-body.js';
import {sparseProduct} from '../src/numerics/sparse.js';
import {solveStreamtubeIses} from '../src/euler/streamtube-ises-update.js';
import {adjustStreamtubeInlets,reparameterizeStreamtubeInlets} from '../src/geometry/streamtube-grid-maintenance.js';
import {streamtubeGridConvexity} from '../src/geometry/streamtube-convex-step.js';

test('RAE128x7 harmonic inlet stall replays exactly without a finite zero-step repair jump', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae128x7-harmonic-before-inlet-step18.json'));
 const before=structuredClone(cp), system=createStreamtubeBodySystem(cp.input);
 const x=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes), value=system.evaluate(x);
 assert.equal(cp.input.mach,.74); assert.equal(cp.input.alpha,2.68);
 assert.equal(system.layout.nx,756); assert.deepEqual(system.layout.tubes,[10,10]);
 assert.deepEqual(Array.from(value.residual),cp.residual);
 const repaired=reparameterizeStreamtubeInlets(value.nodes,cp.input.bodies,cp.continuation.fractions);
 const quality=streamtubeGridConvexity(repaired.nodes);
 assert.equal(quality.valid,true);
 const other=createStreamtubeBodySystem(cp.input), candidate=other.evaluate(other.adoptGeometry(x,repaired.nodes));
 const squared=residual=>residual.reduce((sum,r)=>sum+r*r,0);
 assert.ok(Math.abs(squared(candidate.residual)/squared(value.residual)-1)<1e-6,
  'This inlet failure must not be confused with the finite wake-coordinate merit jump.');
 assert.deepEqual(cp,before);
});

test('RAE128x7 partial harmonic inlet candidate preserves boundaries and flow unknowns while decreasing the original residual', () => {
 const root='docs/solver-reliability/rae-inviscid-audit/fixtures/';
 const source=JSON.parse(fs.readFileSync(root+'rae128x7-harmonic-before-inlet-step18.json'));
 const candidate=JSON.parse(fs.readFileSync(root+'rae128x7-partial-harmonic-inlet-recovery.json'));
 const before=structuredClone({source,candidate});
 assert.deepEqual(candidate.input,source.input);
 assert.deepEqual(candidate.continuation,source.continuation);
 const system=createStreamtubeBodySystem(source.input),layout=system.layout;
 assert.equal(layout.nx,756); assert.deepEqual(layout.tubes,[10,10]);
 assert.deepEqual(candidate.initialEuler.x.slice(0,layout.densityCount),source.initialEuler.x.slice(0,layout.densityCount));
 assert.deepEqual(candidate.initialEuler.x.slice(layout.globalOffset),source.initialEuler.x.slice(layout.globalOffset));
 for(let g=0;g<source.initialEuler.nodes.length;g++) {
  const nodes=source.initialEuler.nodes[g],other=candidate.initialEuler.nodes[g];
  for(let i=0;i<nodes.length;i++)for(let j=0;j<nodes[i].length;j++)
   if(!i||i===nodes.length-1||!j||j===nodes[i].length-1)assert.deepEqual(other[i][j],nodes[i][j]);
 }
 const replay=solveStreamtubeIses(undefined,{...candidate.continuation,resume:candidate,maxIterations:0,retainCheckpoint:true});
 assert.equal(replay.finalQuality.valid,true);
 assert.ok(replay.finalQuality.minCornerSine>.03);
 const squared=r=>r.reduce((s,v)=>s+v*v,0);
 assert.ok(squared(replay.checkpoint.residual)<.9*squared(source.residual));
 assert.ok(replay.diagnostics.residual<.8);
 assert.equal(replay.converged,false,'An improved initial state is not a requested-law root.');
 assert.equal(replay.solverInput.upwind.mucon,-2);
 assert.equal(replay.solverInput.upwind.mcrit,.75);
 assert.deepEqual({source,candidate},before);
});

test('RAE128x7 outer-boundary contact admits a decreasing inverse-coordinate correction without changing flow unknowns', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae128x7-after-inlet-correction-farfield-stall.json'));
 const before=structuredClone(cp);
 const replay=solveStreamtubeIses(undefined,{...cp.continuation,resume:cp,maxIterations:0});
 assert.equal(replay.converged,false);
 assert.equal(replay.finalQuality.valid,true);
 assert.ok(replay.finalQuality.minCornerSine<1e-10);
 assert.ok(replay.diagnostics.residual>.7);
 const audit=auditEulerTangentialMovement(cp,{pairs:[5],correctionScale:1/512});
 const base=audit.rows.find(r=>r.stage==='base'),inverse=audit.rows.find(r=>r.stage==='inverse');
 assert.ok(inverse,JSON.stringify(audit.rows));
 assert.equal(inverse.quality.valid,true);
 assert.ok(inverse.quality.minCornerSine>.03);
 assert.ok(inverse.squaredResidual<.96*base.squaredResidual);
 assert.ok(inverse.maximumNodeMotion<.002);
 assert.ok(inverse.coordinates.every(c=>c.minimumLabelSpan>0));
 assert.equal(audit.sourceUnchanged,true);
 assert.deepEqual(cp,before);
});

test('RAE64x11 first harmonic-grid backtrack preserves the full Euler Jacobian and maintained descent', t => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae64x11-harmonic-before-step5.json'));
 const before=structuredClone(cp),system=createStreamtubeBodySystem(cp.input);
 const x=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes),value=system.evaluate(x);
 assert.equal(cp.input.mach,.74); assert.equal(cp.input.alpha,2.68);
 assert.deepEqual(system.layout.tubes,[14,14]); assert.equal(system.layout.nx,380);
 assert.ok(value.residual.every((v,i)=>v===cp.residual[i]));
 assert.deepEqual(value.nodes,cp.initialEuler.nodes);
 const matrix=system.jacobian(x,{sparse:true});
 const direction=x.map((_,i)=>Math.sin(1.7*i+.4)*(i<system.layout.densityCount?.05:.001));
 const product=sparseProduct(matrix,direction),peak=a=>a.reduce((m,v)=>Math.max(m,Math.abs(v)),0);
 const errors=[1e-4,1e-5,1e-6].map(h=>{
  const plus=system.residual(x.map((v,i)=>v+h*direction[i]));
  const minus=system.residual(x.map((v,i)=>v-h*direction[i]));
  return peak(product.map((v,i)=>v-(plus[i]-minus[i])/(2*h)))/peak(product);
 });
 assert.ok(Math.min(...errors)<1e-6,`First shock backtrack Jv: ${errors}`);
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:1,tolerance:1e-10,traceUpdates:true});
 assert.equal(result.history.length,2,result.reason);
 const accepted=result.history[1],merit=accepted.residualDecrease;
 assert.ok(accepted.step>0&&accepted.step<=1);
 assert.ok(merit.afterSquaredNorm<=merit.allowedSquaredNorm);
 assert.ok(result.linearDiagnostics.maxRelativeResidual<=1e-10);
 assert.equal(streamtubeGridConvexity(result.nodes).valid,true);
 assert.deepEqual(cp,before);
 t.diagnostic(JSON.stringify({relativeDirectionalErrors:errors,step:accepted.step,backtracks:accepted.backtracks,
  beforeSquaredNorm:merit.beforeSquaredNorm,afterSquaredNorm:merit.afterSquaredNorm,
  rawTrial:result.updateTrace[0].stages.find(s=>s.name==='raw-newton')}));
});

test('RAE8x7 standalone outlet-stall checkpoint replays with a valid retained grid and consistent Jacobian', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae8x7-standalone-before-outlet-rejection.json'));
 const before=structuredClone(cp),system=createStreamtubeBodySystem(cp.input);
 const state=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes);
 const value=system.evaluate(state);
 assert.ok(value.residual.every((v,i)=>v===cp.residual[i]));
 assert.equal(streamtubeGridConvexity(value.nodes).valid,true);
 assert.deepEqual(system.layout.tubes,[10,10]);
 assert.equal(system.layout.nx,63);
 assert.equal(cp.input.mach,.74); assert.equal(cp.input.alpha,2.68);
 const matrix=system.jacobian(state,{sparse:true});
 const direction=state.map((_,i)=>Math.sin(1.7*i+.4)*(i < system.layout.densityCount ? .02 : .0001));
 const product=sparseProduct(matrix,direction),peak=a=>a.reduce((m,v)=>Math.max(m,Math.abs(v)),0);
 const errors=[1e-5,1e-6,1e-7,1e-8].map(h=>{
  try {
   const plus=system.residual(state.map((v,i)=>v+h*direction[i]));
   const minus=system.residual(state.map((v,i)=>v-h*direction[i]));
   return peak(product.map((v,i)=>v-(plus[i]-minus[i])/(2*h)))/peak(product);
  } catch(error) {
   // The stalled outlet is close to Mach 1. A difference crossing that
   // physical boundary is not a derivative of the admissible residual.
   if(error.code !== 'streamtube-boundary-subsonic') throw error;
   return Infinity;
  }
 });
 assert.ok(Math.min(...errors)<1e-6,`Outlet-stall Jv: ${errors}`);
 assert.deepEqual(cp,before);
});

for (const [label, fixture, tubes, stagnationMotion] of [
 ['RAE64x11', 'rae64x11-interpolated-euler-root', 14, 'interpolated'],
 ['RAE64x9', 'rae64x9-harmonic-cold-euler-root', 12, 'walls-only'],
]) test(`${label} restored second-order root replays on the complete fine grid`, t => {
 const cp=JSON.parse(fs.readFileSync(`docs/solver-reliability/rae-inviscid-audit/fixtures/${fixture}.json`));
 const before=structuredClone(cp),system=createStreamtubeBodySystem(cp.input);
 assert.equal(cp.input.mach,.74); assert.equal(cp.input.alpha,2.68);
 assert.equal(cp.input.upwind.mucon,1); assert.equal(cp.input.upwind.mcrit,.99);
 assert.equal(cp.input.stagnationMotion,stagnationMotion);
 assert.deepEqual(system.layout.tubes,[tubes,tubes]); assert.equal(system.layout.nx,380);
 const state=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes);
 const value=system.evaluate(state);
 // The production checkpoint guard uses numeric equality: JSON normalizes -0.
 assert.equal(value.residual.length,cp.residual.length);
 assert.ok(value.residual.every((r,i)=>r===cp.residual[i]));
 assert.ok(value.diagnostics.residual<=1e-10);
 assert.deepEqual(value.nodes,cp.initialEuler.nodes);
 const quality=streamtubeGridConvexity(value.nodes);
 assert.equal(quality.valid,true); assert.ok(quality.minCornerSine>.1);
 assert.ok(value.sections.flat(2).every(s=>s.rho>0 && s.p>0
  && Number.isFinite(s.rho) && Number.isFinite(s.p)));
 const matrix=system.jacobian(state,{sparse:true});
 const direction=state.map((_,i)=>Math.sin(1.7*i+.4)*(i<system.layout.densityCount ? .05 : .001));
 const product=sparseProduct(matrix,direction),maximum=a=>a.reduce((m,v)=>Math.max(m,Math.abs(v)),0);
 const errors=[1e-3,1e-4,1e-5].map(h=>{
  const plus=system.residual(state.map((v,i)=>v+h*direction[i]));
  const minus=system.residual(state.map((v,i)=>v-h*direction[i]));
  return maximum(product.map((v,i)=>v-(plus[i]-minus[i])/(2*h)))/maximum(product);
 });
 assert.ok(Math.min(...errors)<1e-6,`restored fine-grid root Jv: ${errors}`);
 assert.deepEqual(cp,before);
 t.diagnostic(JSON.stringify({residual:value.diagnostics.residual,relativeDirectionalErrors:errors}));
});

for (const [label, fixture, mucon] of [
 ['RAE64 before the SMOVE trigger', 'rae64x7-before-smove-trigger', -1],
 ['RAE64x7 cold before the DEKINK pressure failure', 'rae64x7-cold-before-dekink-pressure', 1],
]) test(`${label} has a consistent fine-grid Euler Jacobian`, t => {
 const cp=JSON.parse(fs.readFileSync(`docs/solver-reliability/rae-inviscid-audit/fixtures/${fixture}.json`));
 const before=structuredClone(cp),system=createStreamtubeBodySystem(cp.input);
 assert.equal(cp.input.mach,.74); assert.equal(cp.input.alpha,2.68);
 assert.equal(cp.input.upwind.mucon,mucon);
 assert.deepEqual(system.layout.tubes,[10,10]); assert.equal(system.layout.nx,380);
 const state=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes);
 assert.ok(system.residual(state).every((r,i)=>r===cp.residual[i]));
 assert.equal(streamtubeGridConvexity(system.decode(state).nodes).valid,true);
 const matrix=system.jacobian(state,{sparse:true});
 const maximum=values=>values.reduce((m,v)=>Math.max(m,Math.abs(v)),0),errors={};
 for (const block of ['all','density','geometry']) {
  const direction=state.map((_,i)=>{
   const density=i<system.layout.densityCount;
   if (block==='density' && !density || block==='geometry' && density) return 0;
   return Math.sin(1.7*i+.4)*(density ? .05 : .001);
  });
  const product=sparseProduct(matrix,direction),scale=maximum(product);
  assert.ok(scale>0);
  errors[block]=[1e-3,1e-4,1e-5,1e-6,1e-7].map(h=>{
   const plus=system.residual(state.map((v,i)=>v+h*direction[i]));
   const minus=system.residual(state.map((v,i)=>v-h*direction[i]));
   return maximum(product.map((v,i)=>v-(plus[i]-minus[i])/(2*h)))/scale;
  });
  assert.ok(Math.min(...errors[block])<1e-6,`${block}: ${errors[block]}`);
 }
 assert.deepEqual(cp,before);
 t.diagnostic(JSON.stringify({relativeDirectionalErrors:errors}));
});

test('RAE128 restores inlet spacing along its curve instead of flattening the pinched cell', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae128-inlet-pinch.json'));
 const before=structuredClone(cp),nodes=cp.initialEuler.nodes;
 const moved=adjustStreamtubeInlets(nodes,cp.input.bodies,cp.continuation.fractions,{preserveConvexity:true});
 assert.equal(moved.reparameterization?.method,'polyline-arclength');
 assert.ok(moved.reparameterization.replacedTangentScale<1e-3);
 const quality=streamtubeGridConvexity(moved.nodes);
 assert.equal(quality.valid,true);
 assert.ok(quality.minCornerSine>1e-3,'do not keep approaching the same flat inlet corner');
 const system=createStreamtubeBodySystem(cp.input);
 const value=system.evaluate(system.adoptGeometry(Float64Array.from(cp.initialEuler.x),moved.nodes));
 assert.ok(value.residual.reduce((sum,r)=>sum+r*r,0)<cp.residual.reduce((sum,r)=>sum+r*r,0));
 for(let body=0;body<cp.input.bodies.length;body++) {
  const leading=cp.input.bodies[body].leadingIndex;
  for(let i=0;i<=leading;i++)assert.deepEqual(moved.nodes[body][i].at(-1),moved.nodes[body+1][i][0]);
  for(let i=leading;i<nodes[body].length;i++) {
   assert.deepEqual(moved.nodes[body][i].at(-1),nodes[body][i].at(-1));
   assert.deepEqual(moved.nodes[body+1][i][0],nodes[body+1][i][0]);
  }
 }
 assert.deepEqual(cp,before);
});

test('RAE128 checks the complete inlet-repaired Newton trial before clipping its raw corner', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae128-inlet-pinch.json'));
 const before=structuredClone(cp);
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:1,tolerance:1e-10,traceUpdates:true});
 assert.equal(result.history.length,2,result.reason);
 const h=result.history[1],accepted=result.updateTrace.at(-1);
 assert.equal(accepted.accepted,true);
 assert.equal(accepted.stages.find(s=>s.name==='raw-newton').admissible,false);
 assert.equal(accepted.stages.find(s=>s.name==='inlet-adjustment').admissible,true);
 assert.equal(h.maintenance.inlet.reparameterization.beforeGridGate,true);
 assert.ok(h.step>10000*4.523597214975224e-9,'do not clip to a raw corner which inlet re-spacing repairs');
 assert.ok(h.residualDecrease.afterSquaredNorm<h.residualDecrease.beforeSquaredNorm);
 assert.ok(h.residual<result.history[0].residual);
 assert.equal(result.finalQuality.valid,true);
 assert.deepEqual(result.checkpoint.input,cp.input);
 assert.deepEqual(cp,before);
});

test('RAE128 keeps a valid Newton trial when an earlier trial requested a repair that now fails', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae128-before-step21.json'));
 const before=structuredClone(cp);
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:1,tolerance:1e-10,traceUpdates:true});
 assert.equal(result.history.length,2,result.reason);
 const trials=result.updateTrace,h=result.history[1];
 assert.equal(trials[0].rejection.code,'streamtube-grid-step');
 // The next grid is valid but actual stagnation motion still needs SMOVE.
 // Preserve that trigger; only discard the obsolete geometry-repair request.
 assert.equal(trials[1].stages[0].admissible,true);
 assert.equal(trials[1].redistributionAttempted,true);
 assert.equal(trials.at(-1).accepted,true);
 assert.equal(trials.at(-1).redistributionAttempted,false);
 assert.equal(h.maintenance.discardedGridRepair.code,'streamtube-grid-nonconvex');
 assert.ok(h.step>100*.00016311047670829006,'do not discard a useful Newton step because an earlier grid was invalid');
 assert.ok(h.residual<1.6);
 assert.ok(h.residualDecrease.afterSquaredNorm<h.residualDecrease.beforeSquaredNorm);
 assert.equal(result.finalQuality.valid,true);
 assert.deepEqual(result.checkpoint.input,cp.input);
 assert.deepEqual(cp,before);
});

test('RAE128 backtracks SMOVE without unnecessarily shrinking its decreasing Newton step', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae128-before-step7.json'));
 const before=structuredClone(cp);
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:1,tolerance:1e-10});
 assert.equal(result.history.length,2,result.reason);
 assert.equal(result.reason,'iteration limit');
 const h=result.history[1];
 assert.equal(h.maintenance.residualBacktracking,true);
 assert.ok(h.maintenance.gridCorrectionBacktracks>0);
 assert.ok(h.step>.01,'do not couple every grid backtrack to a smaller Newton step');
 assert.ok(h.residualDecrease.afterSquaredNorm<h.residualDecrease.beforeSquaredNorm);
 assert.equal(result.finalQuality.valid,true);
 assert.deepEqual(result.checkpoint.input,cp.input,'keep the grid and target equations');
 assert.deepEqual(cp,before);
});

test('RAE64 repairs the rejected leading-edge trial before testing its complete grid update', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae64-before-step14.json'));
 const before=structuredClone(cp);
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:1,tolerance:1e-10,traceUpdates:true});
 assert.equal(result.history.length,2,result.reason);
 assert.equal(result.reason,'iteration limit','one accepted update is not a converged solution');
 assert.ok(result.updateTrace.slice(0,13).every(t=>!t.accepted),'preserve the ordinary search first');
 const ordinary=result.updateTrace[12];
 assert.ok(ordinary.stages.find(s=>s.name==='raw-newton').squaredNorm<ordinary.beforeSquaredNorm);
 assert.ok(ordinary.stages.find(s=>s.name==='dekink').squaredNorm>ordinary.beforeSquaredNorm);
 const accepted=result.updateTrace.at(-1);
 assert.equal(accepted.accepted,true);
 assert.equal(accepted.gridRepairRecovery,true);
 assert.equal(accepted.stages.find(s=>s.name==='raw-newton').admissible,false);
 assert.equal(accepted.stages.find(s=>s.name==='dekink').admissible,true);
 assert.ok(accepted.stages.at(-1).squaredNorm<accepted.beforeSquaredNorm);
 assert.ok(accepted.gridCorrectionScale>0&&accepted.gridCorrectionScale<1);
 assert.ok(result.history[1].step>100*ordinary.step);
 assert.equal(result.finalQuality.valid,true);
 assert.equal(result.checkpoint.input.mach,.74);
 assert.equal(result.checkpoint.input.alpha,2.68);
 assert.deepEqual(result.checkpoint.input.upwind,cp.input.upwind);
 assert.deepEqual(result.checkpoint.input.hybrid,cp.input.hybrid);
 assert.deepEqual(cp,before);
});

test('RAE8 backtracks triggered redistribution and reaches the unchanged Euler root', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae8-before-step14.json'));
 const before=structuredClone(cp);
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:8,tolerance:1e-10,traceUpdates:true});
 const trials=result.updateTrace.filter(t=>t.iteration===1);
 assert.equal(trials[0].gridCorrectionScale,1);
 assert.equal(trials[0].rejection.code,'EULER_RESIDUAL_DECREASE');
 const raw=trials[0].stages.find(s=>s.name==='raw-newton');
 assert.ok(raw.squaredNorm<trials[0].beforeSquaredNorm,'Newton alone already decreases the residual');
 assert.ok(trials.at(-1).accepted);
 assert.ok(trials.at(-1).gridCorrectionScale<1,'the rejected grid correction must also shrink');
 assert.equal(result.converged,true,result.reason);
 assert.ok(result.diagnostics.residual<=1e-10);
 assert.equal(result.finalQuality.valid,true);
 assert.equal(result.checkpoint.input.mach,.74);
 assert.equal(result.checkpoint.input.alpha,2.68);
 assert.deepEqual(result.checkpoint.input.upwind,cp.input.upwind);
 assert.deepEqual(result.checkpoint.input.hybrid,cp.input.hybrid);
 assert.deepEqual(cp,before);
});

test('RAE8 with retries disabled retains the rejected source state', () => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae8-before-step14.json'));
 const result=solveStreamtubeIses(undefined,{resume:cp,...cp.continuation,maxIterations:1,maxBacktracks:0,tolerance:1e-10,traceUpdates:true});
 assert.equal(result.converged,false);
 assert.equal(result.updateTrace.length,1);
 assert.equal(result.lastRejectedStep.code,'EULER_RESIDUAL_DECREASE');
 assert.deepEqual(result.checkpoint,cp);
});

test('RAE32 accepts the fifth update after a trial-local grid repair',()=>{
 const traced=reproduceRaeInviscidStep(),plain=reproduceRaeInviscidStep({traceUpdates:false});
 assert.equal(traced.result.converged,false);
 assert.equal(traced.result.reason,'iteration limit');
 assert.equal(traced.result.history.length,2,'the residual-decreasing fifth update must be committed');
 assert.ok(traced.result.history[1].residual<traced.result.history[0].residual);
 assert.deepEqual(traced.result.checkpoint,plain.result.checkpoint);
 assert.equal(traced.report.updateTrace.length,5);
 assert.equal(traced.report.updateTrace[0].rejection.code,'streamtube-grid-step');
 assert.equal(traced.report.updateTrace[1].redistributionAttempted,true);
 assert.ok(traced.report.updateTrace.slice(2).every(t=>t.redistributionAttempted===false));
 assert.equal(traced.report.updateTrace.at(-1).accepted,true);
 assert.ok(Math.abs(traced.result.diagnostics.residual-.2976181014030512)<1e-10);
});

test('RAE32 stalled ninth iterate retains consistent full, density and geometry Jacobians', t => {
 const cp=JSON.parse(fs.readFileSync('docs/solver-reliability/rae-inviscid-audit/fixtures/rae32-before-step10.json'));
 assert.equal(cp.input.mach,.74); assert.equal(cp.input.alpha,2.68);
 const system=createStreamtubeBodySystem(cp.input);
 const x=system.adoptGeometry(Float64Array.from(cp.initialEuler.x),cp.initialEuler.nodes);
 assert.deepEqual(Array.from(system.residual(x)),cp.residual);
 const jacobian=system.jacobian(x,{sparse:true});
 const maximum=values=>values.reduce((m,v)=>Math.max(m,Math.abs(v)),0);
 const errors={};
 for(const block of ['all','density','geometry']) {
  const direction=x.map((_,i)=>{
   const density=i<system.layout.densityCount;
   if ((block === 'density' && !density) || (block === 'geometry' && density)) return 0;
   return Math.sin(1.7 * i + .4) * (density ? .05 : .001);
  });
  const product=sparseProduct(jacobian,direction),scale=maximum(product);
  assert.ok(scale>0);
  errors[block]=[1e-4,1e-5,1e-6].map(h=>{
   const plus=system.residual(x.map((v,i)=>v+h*direction[i]));
   const minus=system.residual(x.map((v,i)=>v-h*direction[i]));
   return maximum(product.map((v,i)=>v-(plus[i]-minus[i])/(2*h)))/scale;
  });
  assert.ok(Math.min(...errors[block])<1e-6,`${block}: ${errors[block]}`);
 }
 t.diagnostic(JSON.stringify({relativeDirectionalErrors:errors}));
});
