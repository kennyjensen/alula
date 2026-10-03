import { test, expect } from '@playwright/test';
import { writeFileSync } from 'node:fs';
import { directBodyConservation } from '../oracles/streamtube-body.js';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

for (const smooth of [false, true]) test(`the actual default GUI ${smooth ? 'with SLOR ' : ''}builds the refined grid before flow and converges with live Newton updates`, async ({ page }) => {
  test.setTimeout(120000);
  const sourceHashes = numericalSourceHashes(['index.html', 'src/ui/app.js', 'src/ui/plots.js', 'src/worker/solver.js']);
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.addInitScript(() => {
    window.quadFrames = [];
    const OriginalWorker = window.Worker;
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => {
          if (data.type === 'mesh') window.quadFrames.push({ mesh: data.mesh, stage: data.stage });
          if (data.type === 'result') window.quadResult = data;
        });
      }
    };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#grid-elliptic').setChecked(smooth);
  await page.locator('#quad-mach').fill('.9');
  await page.locator('#mesh-button').click();
  await expect(page.locator('#status')).toHaveText('Mesh ready', { timeout: 30000 });
  const preview = await page.evaluate(() => window.quadFrames.at(-1).mesh);
  expect(preview.cells).toHaveLength(4587); expect(preview.quality.valid).toBe(true);
  expect(preview.initialization.tubes).toEqual([10, 13, 10]);
  expect(preview.initialization.gridRefinement.normalInterpolation).toBe('streamfunction-quadratic');
  expect(preview.initialization.gridRefinement.flowInitialized).toBe(false);
  expect(Boolean(preview.initialization.gridSmoothing?.converged)).toBe(smooth);
  expect(preview.flow).toBeUndefined(); await expect(page.locator('#flow-legend')).toBeHidden();
  await page.locator('#quad-mach').fill('.2');
  await page.evaluate(() => { window.quadFrames = []; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler converged · research', { timeout: 65000 });
  const { result: r, elapsed, frames } = await page.evaluate(() => ({ ...window.quadResult, frames: window.quadFrames }));
  expect(r.diagnostics.unknowns).toBe(9247); expect(r.diagnostics.equationResidual).toBeLessThan(1e-10);
  expect(r.mesh.quality.valid).toBe(true); expect(r.boundaryLayer).toBeNull();
  expect(Boolean(r.mesh.initialization.gridSmoothing?.converged)).toBe(smooth);
  expect(r.coefficientStatus).toBe('research-unvalidated');
  expect(r.solverSettings.streamwiseMode).toBe('isentropic');
  expect(r.solverSettings.stepMethod).toBe('density-newton');
  expect(r.solverSettings.stepAcceptance).toBe('admissible');
  const built = frames.filter(f => f.stage === 'initial' && !f.mesh.flow).at(-1).mesh;
  expect(built.vertices).toEqual(preview.vertices); expect(built.cells).toEqual(preview.cells);
  const initial = frames.filter(f => f.stage === 'initial' && f.mesh.iteration?.iteration === 0).at(-1).mesh;
  const moving = frames.filter(f => f.mesh.iteration?.iteration > 0).map(f => f.mesh);
  expect(moving.map(m => m.iteration.iteration)).toEqual(r.flow.history.slice(1).map(h => h.iteration));
  let previous = initial;
  for (const mesh of moving) {
    const maxMove = Math.max(...mesh.vertices.map((p, i) => Math.hypot(p.x - previous.vertices[i].x, p.y - previous.vertices[i].y)));
    expect(mesh.iteration.maximumNodeMovement).toBe(maxMove);
    expect(mesh.flow.iteration).toBe(mesh.iteration.iteration);
    expect(mesh.flow.residual).toBe(mesh.iteration.residual);
    previous = mesh;
  }
  expect(r.mesh.vertices).toEqual(moving.at(-1).vertices); expect(r.mesh.flow).toEqual(moving.at(-1).flow);
  expect(directStreamtubeVolumeGeometry(r.flow.nodes).valid).toBe(true);
  const conservation = directBodyConservation(r.flow, r.bodies, r.conditions);
  expect(Math.abs(conservation.balance[0])).toBeLessThan(2e-12);
  expect(Math.abs(conservation.balance[3])).toBeLessThan(2e-11);
  expect(Math.max(...conservation.cutTraction.map(Math.abs))).toBeLessThan(2e-9);
  await expect(page.locator('#flow-iteration')).toContainText('equations converged · research, unvalidated');
  await expect(page.locator('#export-button')).toBeEnabled();
  const download = page.waitForEvent('download'); await page.locator('#export-button').click();
  const stream = await (await download).createReadStream(); let exportedText = ''; for await (const part of stream) exportedText += part;
  const exported = JSON.parse(exportedText).result;
  expect(exported.solverInput).toEqual(r.solverInput);
  expect(Array.isArray(exported.flow.x)).toBe(true); expect(exported.flow.x.length).toBe(r.flow.x.length);
  expect(exported.flow.x.every((x, i) => x === r.flow.x[i])).toBe(true);
  expect(exported.coefficientStatus).toBe('research-unvalidated');
  await page.setViewportSize({ width: 390, height: 844 });
  await page.locator('#show-mesh').uncheck(); await page.locator('#flow-color-mode').selectOption('change');
  await expect(page.locator('#flow-legend')).toBeVisible();
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  const artifact = `/tmp/mses-default-newton${smooth ? '-smoothed' : ''}`;
  await page.locator('.geometry-panel').screenshot({ path: `${artifact}-mobile.png` });
  expect(errors).toEqual([]); expect(changedSources(sourceHashes)).toEqual([]);
  writeFileSync(`${artifact}-gui.json`, JSON.stringify({ date: new Date().toISOString(), physicalAcceptance: false,
    sourceHashes, elapsed, conservation, result: r,
    startup: { nodes: initial.vertices, flow: initial.flow },
    updates: moving.map(m => ({ ...m.iteration, maximumSpeedChangeFromPrevious: m.flow.maximumSpeedChangeFromPrevious })) },
    (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v) + '\n');
  console.log(JSON.stringify({ unknowns: r.diagnostics.unknowns, iterations: r.diagnostics.iterations,
    residual: r.diagnostics.equationResidual, quality: r.mesh.quality, conservation: conservation.balance, elapsed }));
});
