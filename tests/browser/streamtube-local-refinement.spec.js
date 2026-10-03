import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';
import { directChannelConservation } from '../oracles/streamtube.js';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

for (const edgeMatching of ['section-velocity', 'section-velocity-distance']) {
test(`a browser worker locally refines and reconverges all four BLs and two independent wakes: ${edgeMatching}`, async ({page}) => {
  const sourceHashes=numericalSourceHashes(['tests/fixtures/intrinsic-body.js','tests/browser/streamtube-local-refinement.spec.js']);
  await page.goto('/');
  const result=await page.evaluate(edgeMatching=>new Promise((resolve,reject)=>{
    const code=`self.onmessage=async({data:{edgeMatching}})=>{try{
      const {intrinsicBodyFixture}=await import('/tests/fixtures/intrinsic-body.js');
      const {solveCoupledStreamtubeIses}=await import('/src/euler/streamtube-coupled-ises.js');
      const {createCoupledStreamtubeBody}=await import('/src/euler/streamtube-coupled.js');
      const {refineCoupledStreamtubeBody}=await import('/src/euler/streamtube-coupled-refinement.js');
      const input={...intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2}),wakeGeometry:'independent-banks',streamwiseMode:'isentropic'};
      const parent=solveCoupledStreamtubeIses(input,{reynolds:1e5,edgeMatching,tolerance:1e-10,stepAcceptance:'admissible'});
      if(!parent.converged)throw new Error(parent.reason);
      const ne=parent.x.length-4*parent.boundaryLayer.stations.length;
      const source=createCoupledStreamtubeBody(parent.solverInput,{...parent.coupledOptions,initialEuler:parent.flow,initialBL:parent.x.slice(ne)});
      const before=source.evaluate(source.initial),counts=Array(source.euler.layout.nx).fill(1),normal=[[1,4],[4,4],[4,1]],positions=[[0,1,5],[0,4,8],[0,4,5]];
      for(const b of source.euler.layout.bodies){counts[b.leadingIndex-1]=2;counts[b.leadingIndex]=2;counts[b.leadingIndex+1]=3;counts[b.trailingIndex]=2;}
      const refined=refineCoupledStreamtubeBody(parent.solverInput,source,{streamwiseSubdivisions:counts,normalSubdivisions:normal,normalInterpolation:'streamfunction-quadratic'});
      const after=refined.system.evaluate(refined.system.initial),retained=refined.diagnostics.retainedStreamwiseStations;
      let nodeError=0,massError=0;
      before.outer.nodes.forEach((grid,g)=>grid.forEach((row,i)=>row.forEach((p,j)=>{
        const q=after.outer.nodes[g][retained[i]][positions[g][j]];nodeError=Math.max(nodeError,Math.hypot(q.x-p.x,q.y-p.y));
      })));
      before.outer.allocation.groups.forEach((group,g)=>group.forEach((p,j)=>{
        const mass=after.outer.allocation.groups[g].slice(positions[g][j],positions[g][j+1]).reduce((a,t)=>a+t.massFlow,0);
        massError=Math.max(massError,Math.abs(mass/p.massFlow-1));
      }));
      const r=solveCoupledStreamtubeIses({...refined.input,geometryDomain:'positive-simple'},{...refined.options,
        initialEuler:refined.initialEuler,initialBL:refined.initialBL,maxIterations:16,tolerance:1e-10,stepAcceptance:'admissible',iterationGeometry:'ises-sampled'});
      self.postMessage({nodeError,massError,unknowns:refined.system.n,diagnostics:refined.diagnostics,result:r});
    }catch(error){self.postMessage({error:error.stack});}};`;
    const source=code.replaceAll("import('/",`import('${location.origin}/`);
    const url=URL.createObjectURL(new Blob([source],{type:'text/javascript'})),worker=new Worker(url,{type:'module'});
    const close=()=>{worker.terminate();URL.revokeObjectURL(url);};
    worker.onerror=e=>{close();reject(new Error(e.message));};
    worker.onmessage=({data})=>{close();data.error?reject(new Error(data.error)):resolve(data);};worker.postMessage({edgeMatching});
  }),edgeMatching);
  const r=result.result;expect(r.converged).toBe(true);expect(result.unknowns).toBe(1187);
  expect(result.nodeError).toBeLessThan(2e-12);expect(result.massError).toBeLessThan(1e-14);
  expect(r.boundaryLayer.surfaces).toHaveLength(4);expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.coupledOptions.edgeMatching).toBe(edgeMatching);
  expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-10);
  const geometry=directStreamtubeVolumeGeometry(r.flow.nodes);expect(geometry.valid).toBe(true);expect(geometry.concavePrimal).toHaveLength(0);
  const conservation=r.flow.nodes.map((nodes,g)=>directChannelConservation({nodes,sections:r.flow.sections.map(row=>row[g]),cells:r.flow.cells.map(row=>row[g])}));
  for(const c of conservation) for(const key of ['maxLocal','total','internalCancellation'])
    for(const k of key==='internalCancellation'?[0,1,2,3]:[0,3])expect(Math.abs(c[key][k])).toBeLessThan(2e-9);
  expect(changedSources(sourceHashes)).toEqual([]);
  const suffix=edgeMatching==='section-velocity'?'':'-distance';
  fs.writeFileSync(`/tmp/mses-local-coupled-worker${suffix}.json`,JSON.stringify({date:new Date().toISOString(),physicalAcceptance:false,sourceHashes,geometry,conservation,...result},
    (_,v)=>ArrayBuffer.isView(v)?Array.from(v):v)+'\n');
  console.log(JSON.stringify({edgeMatching,unknowns:result.unknowns,families:r.families,quality:r.mesh.quality,nodeError:result.nodeError,massError:result.massError}));
});
}
