import {test,expect} from '@playwright/test';
import fs from 'node:fs';
import {numericalSourceHashes,changedSources} from '../../scripts/validation/provenance.js';
import {directStreamtubeVolumeGeometry} from '../oracles/streamtube-control-volume-geometry.js';
import {directChannelConservation} from '../oracles/streamtube.js';
test('a browser worker solves four surface BLs and two wakes with both outlet banks constrained',async({page})=>{
 const sourceHashes=numericalSourceHashes(['tests/fixtures/intrinsic-body.js','tests/browser/streamtube-wake-outlet.spec.js']);
 await page.goto('/');
 const result=await page.evaluate(()=>new Promise((resolve,reject)=>{
  const code=`self.onmessage=async()=>{try{
   const {intrinsicBodyFixture}=await import('/tests/fixtures/intrinsic-body.js');
   const {createCoupledStreamtubeBody,solveCoupledStreamtubeBody}=await import('/src/euler/streamtube-coupled.js');
   const input={...intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2}),streamwiseMode:'isentropic',wakeGeometry:'independent-banks',wakeOutlet:'banks'};
   const options={reynolds:1e5,edgeMatching:'section-velocity-distance'};
   const s=createCoupledStreamtubeBody(input,options),updates=[];
   const r=solveCoupledStreamtubeBody(s,{maxIterations:10,tolerance:1e-10,onIteration:h=>updates.push(h)});
   self.postMessage({input,options,unknowns:s.n,updates,result:r});
  }catch(error){self.postMessage({error:error.stack});}};`;
  const url=URL.createObjectURL(new Blob([code.replaceAll("import('/",`import('${location.origin}/`)],{type:'text/javascript'}));
  const w=new Worker(url,{type:'module'}),close=()=>{w.terminate();URL.revokeObjectURL(url);};
  w.onerror=e=>{close();reject(new Error(e.message));};w.onmessage=({data})=>{close();data.error?reject(new Error(data.error)):resolve(data);};w.postMessage({});
 }));
 const r=result.result;expect(r.converged).toBe(true);expect(result.unknowns).toBe(362);
 expect(r.boundaryLayer.surfaces).toHaveLength(4);expect(r.boundaryLayer.wakes).toHaveLength(2);
 expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-10);
 expect(Math.max(...r.flow.outletBankTangency.flat().map(Math.abs))).toBeLessThan(2e-10);
 const geometry=directStreamtubeVolumeGeometry(r.flow.nodes);expect(geometry.valid).toBe(true);expect(geometry.concavePrimal).toHaveLength(0);
 const conservation=r.flow.nodes.map((nodes,g)=>directChannelConservation({nodes,sections:r.flow.sections.map(row=>row[g]),cells:r.flow.cells.map(row=>row[g])}));
 for(const c of conservation)for(const key of ['maxLocal','total','internalCancellation'])for(const k of key==='internalCancellation'?[0,1,2,3]:[0,3])expect(Math.abs(c[key][k])).toBeLessThan(2e-9);
 expect(changedSources(sourceHashes)).toEqual([]);
 fs.writeFileSync('docs/current-wake-outlet-browser.json',JSON.stringify({date:new Date().toISOString(),physicalAcceptance:false,sourceHashes,geometry,conservation,...result},(_,v)=>ArrayBuffer.isView(v)?Array.from(v):v)+'\n');
 console.log(JSON.stringify({unknowns:result.unknowns,iterations:result.updates.length-1,families:r.families,quality:r.mesh.quality}));
});
