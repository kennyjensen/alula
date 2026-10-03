import { test, expect } from '@playwright/test';
import { writeFileSync } from 'node:fs';

test('RAE2822 builds the maximum standard smoothed mesh through the real GUI worker', async ({ page }) => {
  test.setTimeout(240000);
  const errors = []; page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.raeMesh = { requests: [], publications: 0 };
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => {
          if (data.type === 'mesh') {
            window.raeMesh.mesh = data.mesh; window.raeMesh.publications++;
          }
          if (data.type === 'mesh-ready' || data.type === 'error') window.raeMesh.terminal = data;
        });
      }
      postMessage(data, ...args) { window.raeMesh.requests.push(data); super.postMessage(data, ...args); }
    };
  });
  await page.goto('/');
  await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#preset').selectOption('rae2822-mses');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#alpha-number').fill('2.68');
  await page.locator('#quad-mach').fill('0.74');
  await page.locator('#grid-intervals').selectOption('128');
  await page.locator('#grid-tubes').selectOption('11');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await page.evaluate(() => { window.raeMesh = { requests: [], publications: 0 }; });
  await page.locator('#mesh-button').click();
  await page.waitForFunction(() => window.raeMesh.terminal, undefined, { timeout: 220000 });
  const result = await page.evaluate(() => window.raeMesh);
  writeFileSync('/tmp/mses-rae2822-max-mesh-browser.json', JSON.stringify(result) + '\n');
  expect(result.requests).toHaveLength(1);
  expect(result.requests[0]).toMatchObject({ meshOnly: true, caseData: {
    gridIntervals: 128, gridTubes: 11, gridEllipticSmoothing: true, mach: .74, alpha: 2.68,
  } });
  expect(result.requests[0].caseData.elements[0].points).toHaveLength(129);
  expect(result.terminal.type, result.terminal.message).toBe('mesh-ready');
  expect(result.mesh.quality.valid).toBe(true);
  expect(result.mesh.quality.invalidCells).toEqual([]);
  expect(result.mesh.initialization.gridSmoothing.converged).toBe(true);
  expect(result.mesh.initialization.gridRefinement.flowInitialized).toBe(false);
  expect(result.mesh.flow).toBeUndefined();
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#mesh-status')).toContainText('SLOR');
  await page.locator('.geometry-panel').screenshot({ path: '/tmp/mses-rae2822-max-mesh-browser.png' });
  expect(errors).toEqual([]);
});
