import { test, expect } from '@playwright/test';
import { naca4, transform } from '../../src/geometry/airfoil.js';
import { directBodyConservation } from '../oracles/streamtube-body.js';

async function openQuad(page) {
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.quadEvents = [];
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args); window.latestSolverWorker = this;
        this.addEventListener('message', ({ data }) => {
          window.quadEvents.push(data);
          // Stop at the first real Euler iteration, without a timeout race
          // or replacing any solver equation or worker response.
          if (window.stopOnIteration && data.type === 'iteration') document.getElementById('stop-button').click();
        });
      }
    };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await expect(page.locator('#solve-button')).toHaveText('Run quad Euler↗');
  await page.locator('#grid-intervals').selectOption('8');
}

test('Run quad Euler displays the initial and moving grids, converges conservative equations and exports research diagnostics', async ({ page }) => {
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await openQuad(page);
  await page.locator('#preset').selectOption('single');
  await page.locator('#alpha-number').fill('0');
  await page.evaluate(() => { window.quadEvents = []; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler converged · research');
  await expect(page.locator('#show-mesh')).toBeChecked();
  await expect(page.locator('#mesh-status')).toContainText('equations converged · research');
  await expect(page.locator('#diagnostics')).toContainText('Euler equation error');
  await expect(page.locator('#diagnostic-note')).toContainText('No coupled boundary layers');
  expect(Number.isFinite(Number(await page.locator('#cl').textContent()))).toBe(true);
  await expect(page.locator('#coefficient-warning')).toBeHidden();
  await expect(page.locator('#bl-panel')).toBeHidden();
  await expect(page.locator('.pressure-panel')).toBeHidden();
  await expect(page.locator('#stop-button')).toBeHidden();
  const events = await page.evaluate(() => window.quadEvents);
  expect(events[0].type).toBe('mesh'); expect(events[0].stage).toBe('initial');
  const initial = events[0].mesh, r = events.at(-1).result;
  expect(initial.initialization.flowSolved).toBe(false);
  expect(events.some(e => e.type === 'iteration' && e.iteration.stage === 'quad-euler')).toBe(true);
  expect(events.some(e => e.type === 'mesh' && e.stage === 'solving')).toBe(true);
  expect(r.status).toBe('research-converged'); expect(r.mach).toBe(.2);
  expect(r.flow.streamwiseMode).toBe('momentum'); expect(r.flow.linearBackend).toBe('klu');
  expect(r.flow.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(r.diagnostics.equationResidual).toBeLessThanOrEqual(1e-10);
  expect(r.mesh.quality.valid).toBe(true); expect(r.mesh.initialization.flowSolved).toBe(true);
  expect(r.mesh.cells.every(c => c.length === 4)).toBe(true);
  expect(r.mesh.cells).toEqual(initial.cells); expect(r.mesh.vertices).not.toEqual(initial.vertices);
  const conservation = directBodyConservation(r.flow, r.bodies, r.conditions);
  expect(Math.max(...conservation.balance.map(Math.abs))).toBeLessThan(2e-9);
  expect(Math.max(...conservation.cutTraction.map(Math.abs))).toBeLessThan(2e-9);
  const download = page.waitForEvent('download'); await page.locator('#export-button').click();
  const stream = await (await download).createReadStream(); let text = ''; for await (const chunk of stream) text += chunk;
  const exported = JSON.parse(text);
  expect(exported.input).toMatchObject({ mach: .2, gridIntervals: 8, gridTubes: 7, flowModel: 'streamtube-grid' });
  expect(Number.isFinite(exported.result.cl)).toBe(true); expect(Number.isFinite(exported.result.cd)).toBe(true); expect(exported.result.boundaryLayer).toBeNull();
  expect(exported.result.coefficients.dragKind).toBe('pressure');
  expect(exported.result.forceStatus).toContain('Unvalidated');
  expect(exported.result.mesh.vertices).toEqual(r.mesh.vertices);
  await page.locator('#quad-mach').fill('0.15');
  await expect(page.locator('#status')).toHaveText('Inputs changed');
  await expect(page.locator('#export-button')).toBeDisabled();
  await expect(page.locator('#mesh-status')).toHaveText('Build mesh or run quad Euler to view the grid.');
  expect(errors).toEqual([]);
});

test('quad Stop terminates real Euler iteration, keeps the mesh and rejects late worker results', async ({ page }) => {
  await openQuad(page); await page.locator('#preset').selectOption('single');
  await page.locator('#alpha-number').fill('0');
  await page.evaluate(() => { window.quadEvents = []; window.stopOnIteration = true; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Stopped');
  await expect(page.locator('#mesh-status')).toContainText('stopped · unsolved');
  await expect(page.locator('#export-button')).toBeDisabled();
  await expect(page.locator('#solve-button')).toBeEnabled();
  expect(await page.evaluate(() => window.quadEvents[0].type)).toBe('mesh');
  const canvas = page.locator('#geometry-canvas'), image = await canvas.evaluate(c => c.toDataURL());
  await page.locator('#show-mesh').uncheck(); expect(await canvas.evaluate(c => c.toDataURL())).not.toBe(image);
  await page.locator('#show-mesh').check();
  await page.evaluate(() => window.latestSolverWorker.dispatchEvent(new MessageEvent('message', {
    data: { id: window.quadEvents[0].id, type: 'result', result: { model: 'research-streamtube-euler', status: 'research-converged' } }
  })));
  await expect(page.locator('#status')).toHaveText('Stopped');
  await page.locator('#alpha-number').fill('1');
  await expect(page.locator('#mesh-status')).toHaveText('Build mesh or run quad Euler to view the grid.');
});

test('multielement quad nonconvergence remains explicit and displays provisional coefficients with all element grids', async ({ page }) => {
  await openQuad(page);
  const contours = [naca4('0012', 160), transform(naca4('0012', 160), { chord: .3, x: -.4, y: .2 })];
  const blade = 'Two-body quad convergence limit\n' + contours.map(points => points.map(p => `${p.x} ${p.y}`).join('\n')).join('\n999 999\n') + '\n';
  await page.locator('#import-file').setInputFiles({ name: 'blade.quad', mimeType: 'text/plain', buffer: Buffer.from(blade) });
  await page.locator('#alpha-number').fill('2');
  await page.evaluate(() => { window.quadEvents = []; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler unconverged', { timeout: 30000 });
  await expect(page.locator('#error-message')).toContainText('Quad Euler did not converge');
  await expect(page.locator('#mesh-status')).toContainText('unconverged');
  await expect(page.locator('#export-button')).toBeEnabled();
  expect(Number.isFinite(Number(await page.locator('#cl').textContent()))).toBe(true);
  await expect(page.locator('#coefficient-warning')).toContainText('Unconverged');
  const r = await page.evaluate(() => window.quadEvents.at(-1).result);
  expect(r.status).toBe('unconverged'); expect(r.flow.surfaces).toHaveLength(4);
  expect(r.mesh.initialization.groups).toBe(3); expect(r.mesh.initialization.flowSolved).toBe(false);
  expect(r.mesh.cells.every(c => c.length === 4)).toBe(true);
  expect(r.diagnostics.equationResidual).toBeGreaterThan(1e-10);
});

test('quad gas initialization failure retains the initial mesh and never reports a solved flow', async ({ page }) => {
  await openQuad(page); await page.locator('#preset').selectOption('single');
  await page.locator('#alpha-number').fill('0'); await page.locator('#quad-mach').fill('0.9');
  await page.evaluate(() => { window.quadEvents = []; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#error-message')).toContainText('Quad Euler initialization failed');
  await expect(page.locator('#mesh-status')).toContainText('unconverged');
  await expect(page.locator('#export-button')).toBeDisabled(); await expect(page.locator('#solve-button')).toBeEnabled();
  expect(await page.evaluate(() => window.quadEvents.map(e => e.type))).toEqual(['mesh', 'error']);
});

test('128 surface intervals refine both foils while inlet and wake counts stay independently fixed', async ({ page }) => {
  test.setTimeout(90000);
  await openQuad(page);
  await page.locator('#grid-intervals').selectOption('128');
  await expect(page.locator('#grid-inlet option[value="auto"]')).toHaveText('Auto · 256');
  await page.locator('#grid-inlet').selectOption('32'); await page.locator('#grid-outlet').selectOption('32');
  await page.evaluate(() => { window.quadEvents = []; });
  await page.locator('#mesh-button').click();
  await expect(page.locator('#status')).toHaveText('Mesh ready', { timeout: 60000 });
  const events = await page.evaluate(() => window.quadEvents), mesh = events[0].mesh;
  expect(events.map(e => e.type)).toEqual(['mesh', 'mesh-ready']);
  expect(mesh.quality.valid).toBe(true); expect(mesh.initialization.flowSolved).toBe(false);
  const spacing = mesh.initialization.gridSpacing;
  expect(spacing.surfaceIntervals).toBe(128); expect(spacing.surfaceIntervalsByElement).toHaveLength(2);
  expect(spacing.surfaceIntervalsByElement.map(b => b.element).sort()).toEqual([0, 1]);
  expect(spacing.surfaceIntervalsByElement.every(b => b.intervals >= 128)).toBe(true);
  expect(spacing.inlet.intervals).toBe(32); expect(spacing.outlet.intervals).toBe(32);
  await expect(page.locator('#cl')).toHaveText('—'); await expect(page.locator('#export-button')).toBeDisabled();
  await page.locator('#grid-intervals').selectOption('64');
  await expect(page.locator('#status')).toHaveText('Inputs changed');
  await expect(page.locator('#grid-inlet')).toHaveValue('32'); await expect(page.locator('#grid-outlet')).toHaveValue('32');
});

test('coarse-panel Build mesh displays the original grid and then a constrained repair without Euler iteration', async ({ page }) => {
  await openQuad(page);
  await page.locator('#grid-intervals').selectOption('16'); await page.locator('#grid-tubes').selectOption('11');
  await page.locator('#resolution').selectOption('80');
  await page.evaluate(() => { window.quadEvents = []; });
  await page.locator('#mesh-button').click();
  await expect(page.locator('#status')).toHaveText('Mesh ready', { timeout: 30000 });
  await expect(page.locator('#mesh-status')).toContainText('initial · unsolved · repaired');
  const events = await page.evaluate(() => window.quadEvents);
  expect(events.map(e => e.type)).toEqual(['mesh', 'mesh', 'mesh-ready']);
  const raw = events[0].mesh, repaired = events[1].mesh;
  expect(raw.quality.valid).toBe(false); expect(repaired.quality.valid).toBe(true);
  expect(repaired.cells).toEqual(raw.cells); expect(repaired.vertices).not.toEqual(raw.vertices);
  expect(repaired.initialization.gridRepair.converged).toBe(true);
  expect(repaired.initialization.gridRepair.massFluxFraction).toBe(.99);
  expect(repaired.initialization.flowSolved).toBe(false);
  await expect(page.locator('#cl')).toHaveText('—'); await expect(page.locator('#export-button')).toBeDisabled();
});
