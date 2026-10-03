import { test, expect } from '@playwright/test';
import { readFileSync } from 'node:fs';

test('real GUI worker checks the displayed SLOR mesh, exposes failed checks and exports matching highlights', async ({ page }) => {
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker; window.auditEvents = [];
    window.Worker = class extends OriginalWorker { constructor(...args) { super(...args); this.addEventListener('message', ({ data }) => window.auditEvents.push(data)); } };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await expect(page.locator('#grid-audit-button')).toBeDisabled();
  await page.locator('#grid-elliptic').check(); await page.locator('#mesh-button').click();
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#grid-audit-button')).toBeEnabled();
  const before = await page.evaluate(() => window.auditEvents.filter(e => e.type === 'mesh').at(-1).mesh);
  await page.locator('#grid-audit-button').click();
  await expect(page.locator('#grid-audit-summary')).toContainText('grid checks need attention', { timeout: 60000 });
  await expect(page.locator('#grid-audit-overlay')).toHaveValue('spacing');
  const audit = await page.evaluate(() => window.auditEvents.find(e => e.type === 'grid-audit').audit);
  expect(audit.accepted).toBe(false); expect(audit.harmonic).toHaveLength(3);
  // Physical-distance preparation removes the former severe cut jumps.
  // The remaining ratio still exceeds the explicit engineering screen.
  expect(audit.spacing.maximumAdjacentRatio).toBeGreaterThan(audit.spacing.limits.adjacentRatioLimit);
  expect(audit.spacing.maximumAdjacentRatio).toBeLessThan(4);
  expect(audit.spacing.spacingLocations.some(p => p.region === 0 && p.j === 0 || p.region === 2 && p.j === 7)).toBe(false);
  expect(audit.fieldErrors).toHaveLength(before.cells.length);
  expect(await page.evaluate(() => window.auditEvents.some(e => e.type === 'iteration'))).toBe(false);
  const beforeOverlay = await page.locator('#geometry-canvas').screenshot();
  await page.locator('#grid-audit-overlay').selectOption('eta');
  expect((await page.locator('#geometry-canvas').screenshot()).equals(beforeOverlay)).toBe(false);
  await page.locator('#grid-audit-details summary').click();
  await expect(page.locator('#grid-audit-table tbody tr')).toHaveCount(3);
  const downloadPromise = page.waitForEvent('download'); await page.locator('#grid-audit-export').click();
  const download = await downloadPromise, exported = JSON.parse(readFileSync(await download.path()));
  expect(exported.mesh.vertices).toEqual(before.vertices); expect(exported.mesh.cells).toEqual(before.cells);
  expect(exported.audit.spacing).toEqual(audit.spacing); expect(exported.audit.harmonic).toEqual(audit.harmonic);
  await page.locator('.geometry-panel').screenshot({ path: '/tmp/mses-grid-audit-browser.png' }).catch(async () => page.locator('#grid-audit').screenshot({ path: '/tmp/mses-grid-audit-browser.png' }));
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.locator('#grid-audit').screenshot({ path: '/tmp/mses-grid-audit-phone.png' });
  await page.locator('#grid-audit-button').click(); await page.locator('#stop-button').click();
  await expect(page.locator('#grid-audit-summary')).toContainText('stopped'); await expect(page.locator('#grid-audit-export')).toBeDisabled();
  await page.locator('#grid-audit-button').click(); await page.locator('#alpha-number').fill('5');
  await expect(page.locator('#grid-audit-button')).toBeDisabled(); await expect(page.locator('#grid-audit-details')).toBeHidden();
  await expect(page.locator('#status')).toHaveText('Inputs changed'); expect(errors).toEqual([]);
});
