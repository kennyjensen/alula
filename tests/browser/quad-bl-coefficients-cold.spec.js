import {test,expect} from '@playwright/test';import fs from 'node:fs';import assert from 'node:assert/strict';
import {numericalSourceHashes,changedSources,sha256} from '../../scripts/validation/provenance.js';
test('real public Euler/BL worker returns coefficient estimates without changing its numerical root',async({page})=>{
 test.setTimeout(180000);
 const sourceHashes=numericalSourceHashes(['tests/browser/quad-bl-coefficients-cold.spec.js','src/ui/app.js','src/ui/quad-coupled-result.js',
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
  families:window.lastResult.families,restart:window.lastResult.restart}));
 for(const id of ['cl','cd','cm']){expect(Number.isFinite(result[id])).toBe(true);await expect(page.locator('#'+id)).toHaveText(result[id].toFixed(id==='cd'?6:4));}
 expect(result.cd).toBeGreaterThan(0);expect(result.coefficientStatus).toBe('research-unvalidated');expect(result.physicalAcceptance).toBe(false);
 await expect(page.locator('#coefficient-warning')).toContainText('Unvalidated');
 const path='docs/current-coupled-startup-default-browser.json',baseline=JSON.parse(fs.readFileSync(path));
 assert.deepEqual(JSON.parse(JSON.stringify(result.restart,(_,v)=>ArrayBuffer.isView(v)?Array.from(v):v)),baseline.restart);
 assert.deepEqual(result.families,baseline.result.families);delete result.restart;
 expect(errors).toEqual([]);expect(changedSources(sourceHashes)).toEqual([]);
 await page.screenshot({path:'docs/current-quad-bl-coefficients-cold.png',fullPage:true});
 fs.writeFileSync('docs/current-quad-bl-coefficients-cold.json',JSON.stringify({date:new Date().toISOString(),sourceHashes,
  scope:'Actual cold public worker/GUI solve with research coefficients. Complete numerical restart and residual families equal the previously audited root exactly.',
  baseline:{path,sha256:sha256(path)},completeRestartExact:true,physicalAcceptance:false,result})+'\n');
 console.log(JSON.stringify(result));
});
