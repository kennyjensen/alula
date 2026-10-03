import {test,expect} from '@playwright/test';import fs from 'node:fs';
import {numericalSourceHashes,changedSources,sha256} from '../../scripts/validation/provenance.js';
// Presentation regression from an audited actual solve. The unfinished
// variant changes only status to test warning/color handling, not physics.
for(const converged of [true,false])test(`saved coupled coefficients: ${converged?'research':'unconverged amber'} display and export`,async({page})=>{
 const sourceHashes=numericalSourceHashes(['src/ui/app.js','src/ui/quad-coupled-result.js','src/ui/quad-coupled-coefficients.js',
  'src/ui/style.css','index.html','tests/browser/quad-bl-coefficients.spec.js']);
 const path='docs/current-coupled-startup-default-browser.json';
 await page.addInitScript(({converged,path})=>{
  const OriginalWorker=window.Worker;
  window.Worker=class extends OriginalWorker{
   postMessage(data,...rest){
    if(!data.caseData?.quadBoundaryLayers)return super.postMessage(data,...rest);
    Promise.all([fetch('/'+path).then(r=>r.json()),import('/src/ui/quad-coupled-result.js')]).then(([saved,adapter])=>{
     const raw=saved.result;raw.converged=converged;if(!converged)raw.reason='Presentation-only unfinished-state fixture';
     const result=adapter.quadCoupledResultForDisplay(raw,saved.input);window.coefficientFixture=result;
     this.dispatchEvent(new MessageEvent('message',{data:{id:data.id,type:'result',result,elapsed:0}}));
    });
   }
  };
 },{converged,path});
 await page.goto('/');await expect(page.locator('#status')).toContainText('Solved',{timeout:45000});
 await page.locator('#flow-model').selectOption('streamtube-bl');await page.locator('#solve-button').click();
 await expect(page.locator('#status')).toHaveText(converged?'Euler/BL converged · research':'Euler/BL unconverged');
 const result=await page.evaluate(()=>({cl:window.coefficientFixture.cl,cd:window.coefficientFixture.cd,cm:window.coefficientFixture.cm,
  coefficientStatus:window.coefficientFixture.coefficientStatus,method:window.coefficientFixture.coefficients.method}));
 for(const id of ['cl','cd','cm']){
  await expect(page.locator('#'+id)).toHaveText(result[id].toFixed(id==='cd'?6:4));
  const metric=page.locator('#'+id).locator('..');
  if(converged)await expect(metric).not.toHaveClass(/provisional/);else{
   await expect(metric).toHaveClass(/provisional/);await expect(page.locator('#'+id)).toHaveCSS('color','rgb(239, 189, 104)');
  }
 }
 await expect(page.locator('#coefficient-warning')).toContainText(converged?'Unvalidated':'Unconverged');
 await expect(page.locator('#coefficient-warning')).toContainText('solid contour');await expect(page.locator('#drag-description')).toContainText('Wake momentum');
 const download=page.waitForEvent('download');await page.locator('#export-button').click();const stream=await(await download).createReadStream();let text='';for await(const chunk of stream)text+=chunk;
 const exported=JSON.parse(text);for(const id of ['cl','cd','cm'])expect(exported.result[id]).toBe(result[id]);
 expect(exported.result.coefficientStatus).toBe(converged?'research-unvalidated':'unconverged');expect(exported.result.physicalAcceptance).toBe(false);
 const suffix=converged?'research':'unconverged';await page.screenshot({path:`docs/current-quad-bl-coefficients-${suffix}.png`,fullPage:true});
 await page.locator('#quad-reynolds').fill('1100000');await expect(page.locator('#status')).toHaveText('Inputs changed');await expect(page.locator('#export-button')).toBeDisabled();
 expect(changedSources(sourceHashes)).toEqual([]);
 fs.writeFileSync(`docs/current-quad-bl-coefficients-${suffix}.json`,JSON.stringify({date:new Date().toISOString(),sourceHashes,
  fixture:{path,sha256:sha256(path)},scope:'Saved actual coupled solution passed through the current display adapter and GUI. Unconverged case changes status only to check presentation; no new numerical solve.',physicalAcceptance:false,result})+'\n');
});
