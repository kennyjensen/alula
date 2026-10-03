import { test, expect } from '@playwright/test';
import { writeFileSync } from 'node:fs';

test('unchecked Smooth grid displays the unmodified initial mesh from the real worker', async ({ page }) => {
  test.setTimeout(60000);
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.meshEvents = []; window.meshRequests = [];
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => window.meshEvents.push(data));
      }
      postMessage(data, ...args) { window.meshRequests.push(data); super.postMessage(data, ...args); }
    };
  });
  await page.goto('/');
  await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  // Explicitly select the unsmoothed control from the current smooth default.
  await page.locator('#grid-elliptic').uncheck();
  await page.evaluate(() => { window.meshEvents = []; window.meshRequests = []; });
  await page.locator('#mesh-button').click();
  await page.waitForFunction(() => window.meshEvents.some(e => e.type === 'mesh-ready' || e.type === 'error'),
    undefined, { timeout: 45000 });
  const outcome = await page.evaluate(() => ({ requests: window.meshRequests, events: window.meshEvents,
    status: document.querySelector('#status').textContent,
    meshStatus: document.querySelector('#mesh-status').textContent }));
  writeFileSync('/tmp/mses-unsmoothed-browser.json', JSON.stringify(outcome, null, 2) + '\n');
  await page.locator('.geometry-panel').screenshot({ path: '/tmp/mses-unsmoothed-browser.png' });
  const box = await page.locator('#geometry-canvas').boundingBox();
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.wheel(0, 750);
  await page.locator('.geometry-panel').screenshot({ path: '/tmp/mses-unsmoothed-browser-domain.png' });
  expect(outcome.requests).toHaveLength(1);
  expect(outcome.requests[0]).toMatchObject({ meshOnly: true, caseData: { gridEllipticSmoothing: false } });
  expect(outcome.events.map(e => [e.type, e.stage])).toEqual([['mesh', 'initial'], ['mesh-ready', undefined]]);
  const mesh = outcome.events[0].mesh;
  expect(mesh.quality.valid).toBe(true);
  expect(mesh.initialization.gridSmoothing).toBeUndefined();
  expect(mesh.initialization.gridRepair).toBeUndefined();
  expect(mesh.initialization.flowSolved).toBe(false);
  expect(mesh.initialization.cutStationSpacing).toBe('physical-x');
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#mesh-status')).not.toContainText('SLOR smoothed');
  await expect(page.locator('#mesh-status')).toContainText('unsmoothed');
  await expect(page.locator('#mesh-status')).toContainText('SLOR not run');
  await expect(page.locator('#cl')).toHaveText('—');
  expect(errors).toEqual([]);
});
