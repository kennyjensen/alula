import {test,expect} from '@playwright/test';
import fs from 'node:fs';
import {numericalSourceHashes,changedSources} from '../../scripts/validation/provenance.js';
import {directStreamtubeVolumeGeometry} from '../oracles/streamtube-control-volume-geometry.js';
import {directChannelConservation} from '../oracles/streamtube.js';
for(const {smooth,tubes} of [{smooth:false,tubes:7},{smooth:true,tubes:7},{smooth:false,tubes:9}])test(`cold two-element Euler/BL assembly in a browser worker, ${tubes} tubes, SLOR ${smooth?'on':'off'}`,async({page})=>{
 test.setTimeout(240000);
 const sourceHashes=numericalSourceHashes(['tests/browser/streamtube-cold-coupled.spec.js','index.html','src/ui/app.js']);
 await page.goto('/');
 const report=await page.evaluate(settings=>new Promise((resolve,reject)=>{
  const code=`self.onmessage=async({data:{smooth,tubes}})=>{try{
   const {naca4,transform}=await import('/src/geometry/airfoil.js');
   const {solveCoupledStreamtubeAssembly}=await import('/src/euler/streamtube-coupled-assembly.js');
   const input={elements:[{name:'Main element',points:naca4('2412',160)},
    {name:'Flap',points:transform(naca4('0012',160),{chord:.3,x:.94,y:-.08,angle:-15})}],
    alpha:4,mach:.2,referenceChord:1,reynolds:1e6,ncrit:9,materialTrips:[[.05,.05],[.05,.05]],
    gridIntervals:16,gridTubes:tubes,gridEllipticSmoothing:smooth};
   const phases=[],iterations=[],frames=[],checkpoints={};let previous=null;
   const started=performance.now();
   const result=solveCoupledStreamtubeAssembly(input,{
    onStage:s=>{phases.push(s.stage);previous=null;},onIteration:h=>iterations.push(h),
    onMesh:(mesh,phase,stage)=>{let movement=null;if(previous?.length===mesh.vertices.length)
      movement=Math.max(...mesh.vertices.map((p,i)=>Math.hypot(p.x-previous[i].x,p.y-previous[i].y)));
      frames.push({stage,phase,iteration:mesh.iteration?.iteration,startupAttempt:mesh.iteration?.startupAttempt,quality:mesh.quality,movement,flow:!!mesh.flow});previous=mesh.vertices;},
    onCheckpoint:q=>{if(q.stage==='boundary-layer-initialization'){checkpoints.seed=q;checkpoints['seed'+q.startupAttempt]=q;}}
   });
   self.postMessage({input,phases,iterations,frames,checkpoints,result,seconds:(performance.now()-started)/1000});
  }catch(error){self.postMessage({error:error.stack,stage:error.stage});}};`;
  const url=URL.createObjectURL(new Blob([code.replaceAll("import('/",`import('${location.origin}/`)],{type:'text/javascript'})),w=new Worker(url,{type:'module'});
  const close=()=>{w.terminate();URL.revokeObjectURL(url);};w.onerror=e=>{close();reject(new Error(e.message));};
  w.onmessage=({data})=>{close();data.error?reject(new Error(data.stage+': '+data.error)):resolve(data);};w.postMessage(settings);
 }),{smooth,tubes});
 const r=report.result;expect(r.converged,r.reason).toBe(true);
 expect(r.model).toBe('research-streamtube-euler-bl');expect(r.physicalAcceptance).toBe(false);
 const attempts=r.initialization.attempts;
 expect(report.phases).toEqual(['euler',...attempts.flatMap(()=>['boundary-layer-initialization','coupled'])]);
 expect(attempts).toHaveLength(tubes===9?2:1);expect(attempts.at(-1).converged).toBe(true);
 if(tubes===9){
  expect(attempts[0].converged).toBe(false);expect(attempts[0].quality.valid).toBe(false);
  expect(attempts[1].thicknessFactor).toBe(.25*attempts[0].thicknessFactor);
  const first=report.checkpoints.seed1.restart,second=report.checkpoints.seed2.restart;
  expect(second.options).toEqual(first.options);expect(second.input).toEqual(first.input);
  for(let k=0;k<first.initialBL.length;k++)expect(second.initialBL[k]).toBe(first.initialBL[k]*(k%4===1||k%4===2?.25:1));
 }
 expect(r.x).toHaveLength(tubes===9?12184:10510);expect(r.boundaryLayer.surfaces).toHaveLength(4);expect(r.boundaryLayer.wakes).toHaveLength(2);
 expect(r.elementOrder).toEqual([1,0]);expect(r.kernelReynolds/r.solverLength).toBeCloseTo(r.referenceReynolds/r.referenceChord,7);
 expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-10);
 expect(r.cl).toBeNull();expect(r.cd).toBeNull();expect(r.coefficientStatus).toBe('unavailable');
 expect(r.solverSettings.wakeOutlet).toBe('banks');expect(r.restart.input.wakeGeometry).toBe('independent-banks');
 expect(Boolean(r.initialization.euler.gridSmoothing?.converged)).toBe(smooth);
 const preview=report.frames.find(f=>f.stage==='euler'&&!f.flow);expect(preview.quality.valid).toBe(true);
 const moving=report.frames.filter(f=>f.stage==='coupled'&&f.iteration>0);
 expect(moving.length).toBe(attempts.reduce((n,a)=>n+a.iterations,0));
 expect(moving.filter(f=>f.startupAttempt===attempts.length)).toHaveLength(r.history.length-1);
 expect(moving.some(f=>f.movement>1e-7)).toBe(true);
 const geometry=directStreamtubeVolumeGeometry(r.flow.nodes);expect(geometry.valid).toBe(true);expect(geometry.concavePrimal).toHaveLength(0);
 const conservation=r.flow.nodes.map((nodes,g)=>directChannelConservation({nodes,sections:r.flow.sections.map(row=>row[g]),cells:r.flow.cells.map(row=>row[g])}));
 for(const c of conservation)for(const k of [0,3])for(const key of ['total','maxLocal'])expect(Math.abs(c[key][k])).toBeLessThan(1e-7);
 for(const c of conservation)expect(Math.max(...c.internalCancellation.map(Math.abs))).toBeLessThan(2e-9);
 expect(Math.max(...r.flow.outletBankTangency.flat().map(Math.abs))).toBeLessThan(2e-10);
 expect(changedSources(sourceHashes)).toEqual([]);
 const suffix=tubes===9?'nine-tubes':smooth?'smoothed':'default';
 fs.writeFileSync('docs/current-coupled-startup-'+suffix+'-browser.json',JSON.stringify({date:new Date().toISOString(),physicalAcceptance:false,sourceHashes,geometry,conservation,
  ...report,restart:r.restart,seed:report.checkpoints.seed.restart,initial:{families:report.checkpoints.seed.families}},(_,v)=>ArrayBuffer.isView(v)?Array.from(v):v)+'\n');
 console.log(JSON.stringify({smooth,tubes,seconds:report.seconds,euler:r.initialization.euler.iterations,coupled:r.history.length-1,attempts:attempts.map(a=>({thicknessFactor:a.thicknessFactor,converged:a.converged,iterations:a.iterations})),
  families:r.families,quality:r.mesh.quality,initialThicknessFactor:r.initialization.boundaryLayer.thicknessFactor}));
});
