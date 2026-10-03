import {test,expect} from '@playwright/test';
import fs from 'node:fs';import assert from 'node:assert/strict';
import {directStreamtubeVolumeGeometry} from '../oracles/streamtube-control-volume-geometry.js';
import {directChannelConservation} from '../oracles/streamtube.js';
import {numericalSourceHashes,changedSources,sha256} from '../../scripts/validation/provenance.js';
const sources=()=>numericalSourceHashes(['tests/browser/quad-bl-gui.spec.js','index.html','src/ui/app.js','src/ui/plots.js','src/ui/style.css','src/ui/quad-coupled-result.js','src/ui/quad-coupled-coefficients.js','src/ui/quad-mesh-progress.js','src/worker/solver.js']);
async function open(page){
 await page.addInitScript(()=>{
  const Worker=window.Worker;window.quadEvents=[];window.quadResult=null;window.lastRequest=null;
  window.Worker=class extends Worker{
   constructor(...args){super(...args);window.latestSolverWorker=this;this.addEventListener('message',({data})=>{
    window.quadEvents.push({id:data.id,type:data.type,stage:data.stage,iteration:data.iteration,
     mesh:data.mesh&&{stage:data.mesh.iteration?.stage,iteration:data.mesh.iteration?.iteration,startupAttempt:data.mesh.iteration?.startupAttempt,
      movement:data.mesh.iteration?.maximumNodeMovement,flow:!!data.mesh.flow,quality:data.mesh.quality,cells:data.mesh.cells.length}});
    if(data.type==='result')window.quadResult=data.result;
    if(window.stopCoupled&&data.type==='mesh'&&data.mesh.iteration?.stage==='coupled'&&data.mesh.iteration.iteration>=1){
     window.stopCoupled=false;setTimeout(()=>document.getElementById('stop-button').click(),0);
    }
   });}
   postMessage(data,...rest){if(data.caseData){if(window.iterationBudget!==undefined&&data.caseData.quadBoundaryLayers)data.caseData.maxIterations=window.iterationBudget;window.lastRequest=data;}return super.postMessage(data,...rest);}
  };
 });
 await page.goto('/');await expect(page.locator('#status')).toContainText('Solved',{timeout:45000});
 await page.locator('#flow-model').selectOption('streamtube-bl');
 // Retain this regression's explicit fixed-trip physical case.
 await page.locator('#quad-transition').selectOption('fixed-trip');
 await page.locator('#grid-elliptic').uncheck();
 await expect(page.locator('#solve-button')).toHaveText('Run quad Euler/BL↗');
 await expect(page.locator('#quad-viscous-conditions')).toBeVisible();await expect(page.locator('#viscous-conditions')).toBeHidden();
 await page.evaluate(()=>{window.quadEvents=[];window.quadResult=null;});
}
for(const smooth of [false,true])test(`public quad Euler/BL cold solve, profiles and export, SLOR ${smooth?'on':'off'}`,async({page})=>{
 test.setTimeout(240000);const errors=[];page.on('pageerror',e=>errors.push(e.message));const sourceHashes=sources();
 await open(page);await page.locator('#grid-elliptic').setChecked(smooth);
 await page.locator('#solve-button').click();
 await expect(page.locator('#status')).toHaveText('Euler/BL converged · research',{timeout:180000});
 await expect(page.locator('#bl-panel')).toBeVisible();await expect(page.locator('.pressure-panel')).toBeVisible();
 await expect(page.locator('#diagnostics')).toContainText('Surface BLs / wakes');await expect(page.locator('#diagnostics')).toContainText('4 / 2');
 await expect(page.locator('#coefficient-warning')).toContainText('Unvalidated');for(const id of ['cl','cd','cm'])await expect(page.locator('#'+id)).not.toHaveText('—');
 await expect(page.locator('#flow-iteration')).toContainText('Euler/BL iteration');await expect(page.locator('#mesh-status')).toContainText('Euler/BL iteration');
 await expect(page.locator('#stop-button')).toBeHidden();await expect(page.locator('#error-message')).toBeHidden();
 const {r,input,events}=await page.evaluate(()=>({r:window.quadResult,input:window.lastRequest.caseData,events:window.quadEvents}));
 expect(input).toMatchObject({quadBoundaryLayers:true,reynolds:1e6,ncrit:9,materialTrips:[[.05,.05],[.05,.05]],gridIntervals:16,gridTubes:7});
 expect(r.converged).toBe(true);expect(r.physicalAcceptance).toBe(false);expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-10);
 const firstMesh=events.findIndex(e=>e.type==='mesh'),firstEuler=events.findIndex(e=>e.type==='iteration');expect(firstMesh).toBeLessThan(firstEuler);
 expect(events.some(e=>e.mesh?.stage==='coupled'&&e.mesh.movement>1e-7)).toBe(true);
 expect(events.filter(e=>e.type==='flow-stage').map(e=>e.stage)).toEqual(['euler','boundary-layer-initialization','coupled']);
 const baselinePath=`docs/current-coupled-startup-${smooth?'smoothed':'default'}-browser.json`,baseline=JSON.parse(fs.readFileSync(baselinePath));
 // The GUI adds presentation/worker plumbing, not equations or a new
 // physical starting state. Compare the complete numerical restart.
 const restart=JSON.parse(JSON.stringify(r.restart,(_,v)=>ArrayBuffer.isView(v)?Array.from(v):v));
 assert.deepEqual(restart,baseline.restart);
 const geometry=directStreamtubeVolumeGeometry(r.flow.nodes);expect(geometry.valid).toBe(true);expect(geometry.concavePrimal).toHaveLength(0);
 const conservation=r.flow.nodes.map((nodes,g)=>directChannelConservation({nodes,sections:r.flow.sections.map(row=>row[g]),cells:r.flow.cells.map(row=>row[g])}));
 for(const c of conservation)for(const k of [0,3])for(const key of ['total','maxLocal'])expect(Math.abs(c[key][k])).toBeLessThan(1e-7);
 for(const c of conservation)expect(Math.max(...c.internalCancellation.map(Math.abs))).toBeLessThan(2e-9);
 const images=new Set();for(const value of ['cf','theta','deltaStar','h']){await page.locator('#bl-quantity').selectOption(value);images.add(await page.locator('#bl-canvas').evaluate(c=>c.toDataURL()));}expect(images.size).toBe(4);
 await page.locator('#flow-color-mode').selectOption('change');await expect(page.locator('#flow-speed-range')).toContainText('Change range');
 await page.locator('#compare-mesh').check();await expect(page.locator('#mesh-update')).toContainText('gold dashed');
 const pressure=page.locator('#pressure-canvas');await pressure.hover();await expect(page.locator('#cp-tooltip')).toContainText('Cp');
 await page.screenshot({path:`docs/current-quad-bl-gui-${smooth?'smoothed':'default'}.png`,fullPage:true});
 const download=page.waitForEvent('download');await page.locator('#export-button').click();const stream=await(await download).createReadStream();let text='';for await(const chunk of stream)text+=chunk;
 const exported=JSON.parse(text);expect(exported.result.restart.initialBL).toBeInstanceOf(Array);expect(exported.result.boundaryLayer.lengthUnit).toBe('case-coordinate');expect(exported.result.numericalBoundaryLayer.stations).toHaveLength(295);
 expect(exported.result.boundaryLayer.surfaces).toHaveLength(4);expect(exported.result.boundaryLayer.wakes).toHaveLength(2);expect([exported.result.cl,exported.result.cd,exported.result.cm].every(Number.isFinite)).toBe(true);expect(exported.result.coefficientStatus).toBe('research-unvalidated');expect(exported.input.materialTrips).toEqual([[.05,.05],[.05,.05]]);
 await page.setViewportSize({width:390,height:844});await page.screenshot({path:`docs/current-quad-bl-gui-${smooth?'smoothed':'default'}-phone.png`,fullPage:true});
 expect(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1)).toBe(true);
 await page.locator('#quad-reynolds').fill('1100000');await expect(page.locator('#status')).toHaveText('Inputs changed');await expect(page.locator('#export-button')).toBeDisabled();await expect(page.locator('#bl-panel')).toBeHidden();
 expect(errors).toEqual([]);expect(changedSources(sourceHashes)).toEqual([]);
 const report={date:new Date().toISOString(),physicalAcceptance:false,sourceHashes,input,events,baseline:{path:baselinePath,sha256:sha256(baselinePath)},numericalRestartExact:true,geometry,conservation,
  result:{converged:r.converged,families:r.families,quality:r.mesh.quality,diagnostics:r.diagnostics,boundaryLayer:r.boundaryLayer,initialization:r.initialization},restart};
 fs.writeFileSync(`docs/current-quad-bl-gui-${smooth?'smoothed':'default'}.json`,JSON.stringify(report)+'\n');
 console.log(JSON.stringify({smooth,families:r.families,coupledIterations:r.history.length-1,profiles:r.boundaryLayer.surfaces.length,wakes:r.boundaryLayer.wakes.length,numericalRestartExact:true}));
});
test('quad BL controls stay distinct, mesh builds before solving, and Stop rejects late coupled results',async({page})=>{
 test.setTimeout(180000);await open(page);
 await page.locator('[data-key="quadTripUpper"]').last().fill('.08');
 await page.locator('#quad-reynolds').fill('');await page.locator('#mesh-button').click();await expect(page.locator('#status')).toHaveText('Mesh ready',{timeout:90000});
 expect(await page.evaluate(()=>window.quadEvents.some(e=>e.type==='iteration'))).toBe(false);
 expect(await page.evaluate(()=>window.lastRequest.caseData.quadBoundaryLayers)).toBe(true);
 await page.locator('#solve-button').click();await expect(page.locator('#error-message')).toContainText('Reynolds');
 await page.locator('#quad-reynolds').fill('1000000');await page.evaluate(()=>{window.stopCoupled=true;window.quadEvents=[];});
 await page.locator('#solve-button').click();await expect(page.locator('#status')).toHaveText('Stopped',{timeout:120000});
 const data=await page.evaluate(()=>({events:window.quadEvents,input:window.lastRequest.caseData}));expect(data.input.materialTrips).toEqual([[.05,.05],[.08,.05]]);
 expect(data.events.some(e=>e.mesh?.stage==='coupled'&&e.mesh.iteration>=1)).toBe(true);
 await expect(page.locator('#mesh-status')).toContainText('stopped');await expect(page.locator('#export-button')).toBeDisabled();
 await page.evaluate(()=>window.latestSolverWorker.dispatchEvent(new MessageEvent('message',{data:{id:window.lastRequest.id,type:'result',result:{model:'research-streamtube-euler-bl'}}})));
 await expect(page.locator('#status')).toHaveText('Stopped');await expect(page.locator('#solve-button')).toBeEnabled();
});
