import { test, expect } from '@playwright/test';
import { naca4 } from '../../src/geometry/airfoil.js';
test('refines a selected element, exports its actual spacing, and enforces the assembly panel budget',async({page})=>{
  await page.goto('/');await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#preset').selectOption('three');await expect(page.locator('#status')).toContainText('Solved');
  for(const [name,panels] of [['Slat','80'],['Main element','160'],['Flap','80']]){
    await page.locator('details').filter({has:page.getByLabel(name+' panels',{exact:true})}).evaluate(node=>{node.open=true;});
    await page.getByLabel(name+' panels',{exact:true}).selectOption(panels);
  }
  await expect(page.locator('#export-button')).toBeDisabled();
  await page.locator('#solve-button').click();await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#refine-target').selectOption('0');await page.locator('#refine-button').click();
  await expect(page.locator('#refinement-result')).toContainText('Slat · 320 → 360 panels');
  const downloaded=page.waitForEvent('download');await page.locator('#export-button').click();
  const stream=await(await downloaded).createReadStream();let text='';for await(const chunk of stream)text+=chunk;
  expect(JSON.parse(text).input.elements.map(e=>e.points.length-1)).toEqual([120,160,80]);
  for(const [name,panels]of [['Main element','240'],['Flap','120'],['Slat','320']]){
    await page.locator('details').filter({has:page.getByLabel(name+' panels',{exact:true})}).evaluate(node=>{node.open=true;});
    await page.getByLabel(name+' panels',{exact:true}).selectOption(panels);
  }
  await expect(page.getByLabel('Slat panels',{exact:true}).locator('option[value="400"]')).toHaveJSProperty('disabled',true);
  await expect(page.locator('#resolution')).toBeDisabled();
  await page.locator('#solve-button').click();await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#run-meta')).toContainText('680 panels');
  await expect(page.locator('#refine-button')).toBeDisabled();
  await page.locator('#refine-target').selectOption('all');await expect(page.locator('#refine-button')).toBeDisabled();
});

test('solves, clears stale results, refines, and exports the actual case', async ({ page }) => {
  const errors = []; page.on('pageerror', error => errors.push(error.message));
  await page.goto('/');
  await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#cl')).not.toHaveText('—');
  const original = await page.locator('#cl').textContent();
  await page.locator('#alpha-number').fill('6');
  await expect(page.locator('#cl')).toHaveText('—');
  await expect(page.locator('#export-button')).toBeDisabled();
  await page.getByRole('button', { name: 'Run analysis' }).click();
  await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#cl')).not.toHaveText(original);
  await page.locator('#refine-button').click();
  await expect(page.locator('#refinement-result')).toContainText('320 → 480 panels');
  const downloaded = page.waitForEvent('download'); await page.locator('#export-button').click();
  const download = await downloaded;
  const stream = await download.createReadStream(); let text = ''; for await (const chunk of stream) text += chunk;
  const exported = JSON.parse(text);
  expect(exported.input.alpha).toBe(6); expect(exported.result.panelCount).toBe(480);
  expect(exported.result.cd).toBeNull(); expect(exported.result.cl.toFixed(4)).toBe(await page.locator('#cl').textContent());
  expect(errors).toEqual([]);
});
test('colliding geometry produces an error and cannot export an old result', async ({ page }) => {
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.getByLabel('Flap LE x').fill('0.4'); await page.getByLabel('Flap LE y').fill('0');
  await page.getByRole('button', { name: 'Run analysis' }).click();
  await expect(page.locator('#error-message')).toBeVisible();
  await expect(page.locator('#export-button')).toBeDisabled();
  await expect(page.locator('#cl')).toHaveText('—');
});
test('presets, streamlines, method dialog and mobile layout work', async ({ page }) => {
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#streamline-status')).toHaveText('Potential-flow streamlines');
  await page.locator('#preset').selectOption('single'); await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#element-count')).toHaveText('1 ELEMENT');
  await page.locator('#preset').selectOption('three'); await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#element-count')).toHaveText('3 ELEMENTS');
  await page.locator('#method-button').click(); await expect(page.locator('dialog')).toBeVisible();
  await page.keyboard.press('Escape'); await expect(page.locator('dialog')).not.toBeVisible();
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  await expect(page.locator('#solve-button')).toBeEnabled();
});
test('imports coordinate files, preserves their spacing, and rejects open trailing edges', async ({ page }) => {
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  const points = naca4('0012', 80);
  const coordinates = name => `${name}\n${points.map(p => `${p.x} ${p.y}`).join('\n')}\n`;
  await page.locator('#import-file').setInputFiles({ name: 'test.dat', mimeType: 'text/plain', buffer: Buffer.from(coordinates('Imported NACA')) });
  await expect(page.locator('#case-title')).toHaveText('Imported NACA');
  await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#run-meta')).toContainText('80 panels');
  await expect(page.locator('#preset')).toHaveValue('custom');
  await expect(page.locator('#refine-button')).toBeDisabled();
  points.at(-1).y = 0.002;
  await page.locator('#import-file').setInputFiles({ name: 'blunt.dat', mimeType: 'text/plain', buffer: Buffer.from(coordinates('Blunt section')) });
  await expect(page.locator('#error-message')).toContainText('sharp trailing edge');
  await expect(page.locator('#export-button')).toBeDisabled();
});

test('coupled flow exposes validated drag, BL plots, trips and reproducible exports', async ({ page }) => {
  const errors = []; page.on('pageerror', error => errors.push(error.message));
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#preset').selectOption('single'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('coupled');
  await expect(page.locator('#viscous-conditions')).toBeVisible();
  await page.locator('#solve-button').click(); await expect(page.locator('#status')).toContainText('Solved');
  expect(Number(await page.locator('#cl').textContent())).toBeCloseTo(.4173,3);
  expect(Number(await page.locator('#cd').textContent())).toBeCloseTo(.007212,5);
  await expect(page.locator('#streamline-status')).toContainText('Computed wake');
  await expect(page.locator('#bl-panel')).toBeVisible();
  await expect(page.locator('#diagnostics')).toContainText('transition · upper / lower');
  await page.locator('#bl-quantity').selectOption('deltaStar');
  await page.locator('#trip-upper').fill('0.05'); await page.locator('#trip-lower').fill('0.1');
  await expect(page.locator('#cd')).toHaveText('—'); await expect(page.locator('#bl-panel')).not.toBeVisible();
  await expect(page.locator('#export-button')).toBeDisabled();
  await page.locator('#solve-button').click(); await expect(page.locator('#status')).toContainText('Solved');
  expect(Number(await page.locator('#cd').textContent())).toBeCloseTo(.011328,5);
  const downloaded = page.waitForEvent('download'); await page.locator('#export-button').click();
  const stream = await (await downloaded).createReadStream(); let text = ''; for await (const chunk of stream) text += chunk;
  const data = JSON.parse(text);
  expect(data.input.reynolds).toBe(1e6); expect(data.input.trips).toEqual([0.05, 0.1]);
  expect(data.result.model).toBe('multielement-coupled-incompressible'); expect(data.result.boundaryLayer.wakes[0].stations.length).toBeGreaterThan(10);
  expect(data.result.diagnostics.equationResidual).toBeLessThan(1e-8);
  expect(data.result.diagnostics.wakeResidual).toBeLessThan(1e-6);
  await page.locator('#refine-button').click(); await expect(page.locator('#refinement-result')).toContainText('ΔCD');
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  expect(errors).toEqual([]);
});

test('coupled worker solves multiple elements and labels failed coefficients as provisional', async ({ page }) => {
  test.setTimeout(65000);
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('coupled'); await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toContainText('Solved',{timeout:45000});
  await expect(page.locator('#streamline-status')).toContainText('Computed wakes (2)');
  await expect(page.locator('#diagnostics')).toContainText('Flap transition');
  await expect(page.locator('#export-button')).toBeEnabled();
  expect(Number(await page.locator('#cl').textContent())).toBeCloseTo(1.7098,3);
  const downloaded=page.waitForEvent('download');await page.locator('#export-button').click();
  const stream=await(await downloaded).createReadStream();let text='';for await(const chunk of stream)text+=chunk;
  const result=JSON.parse(text).result;
  expect(result.boundaryLayer.surfaces.length).toBe(4);expect(result.boundaryLayer.wakes.length).toBe(2);
  expect(result.boundaryLayer.wakes.every(w=>w.stations.length>10)).toBe(true);
  const failure = await page.evaluate(async points => {
    const worker = new Worker('/src/worker/solver.js', { type: 'module' });
    return new Promise(resolve => {
      worker.onmessage = ({ data }) => { if (data.type === 'result' || data.type === 'error') { worker.terminate(); resolve(data); } };
      worker.postMessage({ id: 42, caseData: { flowModel: 'coupled', elements: [{ points }], alpha: 4, maxIterations: 1 } });
    });
  }, naca4('0012', 160));
  expect(failure.result.status).toBe('unconverged'); expect(failure.result.coefficientStatus).toBe('unconverged');
  expect(Number.isFinite(failure.result.cd)).toBe(true); expect(Number.isFinite(failure.result.cl)).toBe(true);
});
