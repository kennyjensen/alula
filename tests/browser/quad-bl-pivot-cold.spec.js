import {test,expect} from '@playwright/test';import fs from 'node:fs';import assert from 'node:assert/strict';
import {numericalSourceHashes,changedSources,sha256} from '../../scripts/validation/provenance.js';
test('public Euler/BL cold solve preserves its default flow after thin-cell arithmetic improvement',async({page})=>{
 test.setTimeout(180000);
 const sourceHashes=numericalSourceHashes(['tests/browser/quad-bl-pivot-cold.spec.js','src/ui/app.js','src/ui/quad-coupled-result.js',
  'src/ui/quad-coupled-coefficients.js','src/ui/quad-mesh-progress.js','src/worker/solver.js','src/ui/style.css','index.html']);
 await page.addInitScript(()=>{const OriginalWorker=window.Worker;window.Worker=class extends OriginalWorker{
  constructor(...args){super(...args);this.addEventListener('message',({data})=>{if(data.type==='result')window.lastResult=data.result;});}
 };});
 const errors=[];page.on('pageerror',e=>errors.push(e.message));
 await page.goto('/');await expect(page.locator('#status')).toContainText('Solved',{timeout:45000});
 await page.locator('#flow-model').selectOption('streamtube-bl');await page.locator('#solve-button').click();
 await expect(page.locator('#status')).toHaveText('Euler/BL converged · research',{timeout:120000});
 const result=await page.evaluate(()=>({cl:window.lastResult.cl,cd:window.lastResult.cd,cm:window.lastResult.cm,coefficients:window.lastResult.coefficients,
  coefficientStatus:window.lastResult.coefficientStatus,physicalAcceptance:window.lastResult.physicalAcceptance,
  families:window.lastResult.families,restart:window.lastResult.restart,
  surfaces:window.lastResult.boundaryLayer.surfaces.length,wakes:window.lastResult.boundaryLayer.wakes.length}));
 for(const id of ['cl','cd','cm']){expect(Number.isFinite(result[id])).toBe(true);await expect(page.locator('#'+id)).toHaveText(result[id].toFixed(id==='cd'?6:4));}
 expect(result.cd).toBeGreaterThan(0);expect(result.coefficientStatus).toBe('research-unvalidated');expect(result.physicalAcceptance).toBe(false);
 await expect(page.locator('#coefficient-warning')).toContainText('Unvalidated');
 const evidencePath='docs/current-coupled-checkpoint-cold.json',evidence=JSON.parse(fs.readFileSync(evidencePath));
 const path=evidence.baseline.path,baseline=JSON.parse(fs.readFileSync(path));
 assert.equal(evidence.completeRestartExact,true);assert.equal(sha256(path),evidence.baseline.sha256);
 const restart=JSON.parse(JSON.stringify(result.restart,(_,v)=>ArrayBuffer.isView(v)?Array.from(v):v));
 assert.deepEqual(restart.input,baseline.restart.input);assert.deepEqual(restart.options,baseline.restart.options);
 let maximumNodeDifference=0,maximumStateDifference=0;
 restart.initialEuler.nodes.forEach((grid,g)=>grid.forEach((row,i)=>row.forEach((p,j)=>{
  const q=baseline.restart.initialEuler.nodes[g][i][j];maximumNodeDifference=Math.max(maximumNodeDifference,Math.hypot(p.x-q.x,p.y-q.y));
 })));
 const old=[...baseline.restart.initialEuler.x,...baseline.restart.initialBL],current=[...restart.initialEuler.x,...restart.initialBL];
 expect(current.length).toBe(old.length);current.forEach((v,i)=>{maximumStateDifference=Math.max(maximumStateDifference,Math.abs(v-old[i]));});
 expect(maximumNodeDifference).toBeLessThan(1e-9);expect(maximumStateDifference).toBeLessThan(1e-9);
 expect(Math.max(...Object.values(result.families))).toBeLessThan(1e-10);expect(result.surfaces).toBe(4);expect(result.wakes).toBe(2);
 for(const key of ['cl','cd','cm'])expect(Math.abs(result[key]-evidence.result[key])).toBeLessThan(1e-9);
 delete result.restart;
 expect(errors).toEqual([]);expect(changedSources(sourceHashes)).toEqual([]);
 await page.screenshot({path:'docs/current-coupled-relative-cold.png',fullPage:true});
 fs.writeFileSync('docs/current-coupled-relative-cold.json',JSON.stringify({date:new Date().toISOString(),sourceHashes,
  scope:'Actual cold public worker/GUI solve after relative-vector geometry arithmetic, with all four BLs and both wakes. Full state/geometry and coefficients agree with the retained pre-edit default within stated numerical bounds; physical accuracy remains unvalidated.',
  baseline:{path,sha256:sha256(path),evidencePath,evidenceHash:sha256(evidencePath)},
  maximumNodeDifference,maximumStateDifference,physicalAcceptance:false,result,restart})+'\n');
 console.log(JSON.stringify(result));
});
