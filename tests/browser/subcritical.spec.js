import { test,expect } from '@playwright/test';
import { naca4,transform } from '../../src/geometry/airfoil.js';
import {slatMainFlap} from '../fixtures/subcritical-assemblies.js';

test('the workbench builds the default main/flap MSES-style quadrilateral grid with 11 tubes without a flow solve',async({page})=>{
  test.setTimeout(90000);
  const errors=[];page.on('pageerror',error=>errors.push(error.message));
  await page.addInitScript(()=>{
    const OriginalWorker=window.Worker;window.gridEvents=[];
    window.Worker=class extends OriginalWorker{constructor(...args){super(...args);
      this.addEventListener('message',({data})=>{window.gridEvents.push(data.type);if(data.type==='mesh')window.actualGrid=data.mesh;});}};
  });
  await page.goto('/');await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await expect(page.locator('#solve-button')).toBeVisible();await expect(page.locator('#streamtube-grid-conditions')).toBeVisible();
  await page.locator('#grid-tubes').selectOption('11');
  await page.evaluate(()=>{window.gridEvents=[];});
  await page.locator('#mesh-button').click();await expect(page.locator('#status')).toHaveText('Mesh ready',{timeout:60000});
  await expect(page.locator('#mesh-status')).toContainText('MSES-style quadrilateral mesh');
  await expect(page.locator('#mesh-status')).toContainText('initial · unsolved');
  await expect(page.locator('#cl')).toHaveText('—');await expect(page.locator('#export-button')).toBeDisabled();
  const {mesh,events}=await page.evaluate(()=>({mesh:window.actualGrid,events:window.gridEvents}));
  expect(events).toEqual(['mesh','mesh-ready']);expect(mesh.topology).toBe('intrinsic-quadrilateral-streamtubes');
  expect(mesh.cells.length).toBeGreaterThan(500);expect(mesh.cells.every(cell=>cell.length===4)).toBe(true);
  expect(mesh.quality.valid).toBe(true);expect(mesh.initialization.flowSolved).toBe(false);
  expect(mesh.initialization.tubes).toEqual([11,11,11]);
  expect(mesh.cells).toHaveLength(3267);
  expect(mesh.initialization.gridSpacing.inlet.intervals).toBe(32);
  expect(mesh.initialization.gridSpacing.outlet.intervals).toBe(32);
  expect(mesh.initialization.elementOrder).toEqual([1,0]);
  const image=await page.locator('#geometry-canvas').evaluate(canvas=>canvas.toDataURL());
  await page.locator('#show-mesh').uncheck();expect(await page.locator('#geometry-canvas').evaluate(canvas=>canvas.toDataURL())).not.toBe(image);
  await page.locator('#show-mesh').check();
  await page.screenshot({path:'/tmp/mses-quadrilateral-workbench.png',fullPage:true});
  await page.locator('#grid-inlet').selectOption('16');
  await expect(page.locator('#mesh-status')).toHaveText('Build mesh or run quad Euler to view the grid.');
  await page.locator('#grid-outlet').selectOption('64');
  await page.evaluate(()=>{window.gridEvents=[];});
  await page.locator('#mesh-button').click();await expect(page.locator('#status')).toHaveText('Mesh ready',{timeout:60000});
  const custom=await page.evaluate(()=>({mesh:window.actualGrid,events:window.gridEvents}));
  expect(custom.events).toEqual(['mesh','mesh-ready']);
  expect(custom.mesh.initialization.gridSpacing.inlet.intervals).toBe(16);
  expect(custom.mesh.initialization.gridSpacing.outlet.intervals).toBe(64);
  expect(custom.mesh.quality.valid).toBe(true);expect(custom.mesh.initialization.flowSolved).toBe(false);
  expect(custom.mesh.cells.length).toBeGreaterThan(mesh.cells.length);
  await page.locator('#grid-tubes').selectOption('9');
  await expect(page.locator('#mesh-status')).toHaveText('Build mesh or run quad Euler to view the grid.');
  await page.setViewportSize({width:390,height:844});expect(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth)).toBe(true);
  expect(errors).toEqual([]);
});

test('station-centered coupling with the optional linear mass starts cold and closes all six BLs',async({page})=>{
  test.setTimeout(300000);
  await page.goto('/');
  const result=await page.evaluate(input=>new Promise((resolve,reject)=>{
    const moduleUrl=new URL('/src/potential/result.js',location.href).href;
    const source=`self.onmessage=async({data})=>{try{
      const {solveSubcriticalAssembly}=await import(${JSON.stringify(moduleUrl)});
      self.postMessage({result:solveSubcriticalAssembly(data,{edgeVelocitySampling:'station-average',bodyMassInterpolation:'linear'})});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url=URL.createObjectURL(new Blob([source],{type:'text/javascript'})),worker=new Worker(url,{type:'module'});
    const close=()=>{worker.terminate();URL.revokeObjectURL(url);};
    worker.onerror=error=>{close();reject(new Error(error.message));};
    worker.onmessage=({data})=>{close();data.error?reject(new Error(data.error)):resolve(data.result);};
    worker.postMessage(input);
  }),slatMainFlap([60,120,60]));
  expect(result.status).toBe('solved');
  expect(result.assemblyInitialization.used).toBe(true);
  expect(result.solverSettings.edgeVelocitySampling).toBe('station-average');
  expect(result.solverSettings.bodyMassInterpolation).toBe('linear');
  expect(result.solverSettings.linearBackend).toBe('klu');
  expect(result.boundaryLayer.surfaces).toHaveLength(6);
  expect(result.boundaryLayer.wakes).toHaveLength(3);
  expect(result.diagnostics.equationResidual).toBeLessThan(1e-8);
  expect(result.diagnostics.wakeResidual).toBeLessThan(1e-6);
  expect(result.diagnostics.massBudgetError).toBeLessThan(1e-10);
  expect(result.diagnostics.maxMach).toBeLessThan(1);
  expect(Math.abs(result.cd-result.boundaryLayer.wakes.reduce((sum,w)=>sum+w.drag,0))).toBeLessThan(1e-12);
});
