import {test,expect} from '@playwright/test';import fs from 'node:fs';
import {numericalSourceHashes,changedSources} from '../../scripts/validation/provenance.js';
test('a bounded real coupled solve displays provisional BL/pressure profiles and exports its unfinished state',async({page})=>{
 test.setTimeout(180000);const errors=[];page.on('pageerror',e=>errors.push(e.message));
 const sourceHashes=numericalSourceHashes(['tests/browser/quad-bl-unconverged.spec.js','src/ui/app.js','src/ui/quad-coupled-result.js','src/ui/quad-coupled-coefficients.js','src/ui/quad-mesh-progress.js','src/worker/solver.js','index.html','src/ui/style.css']);
 await page.addInitScript(()=>{
  const OriginalWorker=window.Worker;
  window.Worker=class extends OriginalWorker{
   constructor(...args){super(...args);this.addEventListener('message',({data})=>{if(data.type==='result')window.quadResult=data.result;});}
   postMessage(data,...rest){if(data.caseData?.quadBoundaryLayers)data.caseData.maxIterations=1;return super.postMessage(data,...rest);}
  };
 });
 await page.goto('/');await expect(page.locator('#status')).toContainText('Solved',{timeout:45000});
 await page.locator('#flow-model').selectOption('streamtube-bl');
 await page.locator('#solve-button').click();await expect(page.locator('#status')).toHaveText('Euler/BL unconverged',{timeout:120000});
 await expect(page.locator('#error-message')).toContainText('profiles are provisional');await expect(page.locator('#coefficient-warning')).toContainText('Unconverged');
 await expect(page.locator('#bl-panel')).toBeVisible();await expect(page.locator('#bl-panel')).toHaveClass(/provisional/);
 await expect(page.locator('#mesh-status')).toContainText('unconverged');await expect(page.locator('#export-button')).toBeEnabled();await expect(page.locator('#solve-button')).toBeEnabled();
 const first=await page.locator('#bl-canvas').evaluate(c=>c.toDataURL());await page.locator('#bl-quantity').selectOption('theta');expect(await page.locator('#bl-canvas').evaluate(c=>c.toDataURL())).not.toBe(first);
 await page.locator('#pressure-canvas').hover();await expect(page.locator('#cp-tooltip')).toContainText('Cp');
 const result=await page.evaluate(()=>({model:window.quadResult.model,converged:window.quadResult.converged,physicalAcceptance:window.quadResult.physicalAcceptance,families:window.quadResult.families,
  surfaces:window.quadResult.boundaryLayer.surfaces.length,wakes:window.quadResult.boundaryLayer.wakes.length,attempts:window.quadResult.initialization.attempts.map(a=>({converged:a.converged,iterations:a.iterations,thicknessFactor:a.thicknessFactor})),
  finiteProfiles:window.quadResult.boundaryLayer.surfaces.every(s=>s.stations.every(p=>[p.theta,p.deltaStar,p.cf,p.cp].every(Number.isFinite)))}));
 expect(result.converged).toBe(false);expect(result.physicalAcceptance).toBe(false);expect(result.surfaces).toBe(4);expect(result.wakes).toBe(2);expect(result.finiteProfiles).toBe(true);
 expect(result.attempts).toHaveLength(2);expect(result.attempts.every(a=>a.iterations===1&&!a.converged)).toBe(true);
 await page.screenshot({path:'docs/current-quad-bl-gui-unconverged.png',fullPage:true});
 const download=page.waitForEvent('download');await page.locator('#export-button').click();const stream=await(await download).createReadStream();let text='';for await(const chunk of stream)text+=chunk;
 const exported=JSON.parse(text);expect(exported.result.converged).toBe(false);expect(exported.result.restart.initialBL).toHaveLength(1180);expect(exported.result.boundaryLayer.surfaces).toHaveLength(4);expect(exported.result.coefficientStatus).toBe('unconverged');
 expect(errors).toEqual([]);expect(changedSources(sourceHashes)).toEqual([]);
 fs.writeFileSync('docs/current-quad-bl-gui-unconverged.json',JSON.stringify({date:new Date().toISOString(),sourceHashes,physicalAcceptance:false,
  scope:'One real Newton update per coupled attempt, deliberately unfinished, through the public worker/UI. No equation or result replacement.',result})+'\n');
 console.log(JSON.stringify(result));
});
