import { test, expect } from '@playwright/test';

async function openGrid(page) {
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.gridEvents = []; window.gridRequests = [];
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => window.gridEvents.push(data));
      }
      postMessage(data, ...args) { window.gridRequests.push(data); super.postMessage(data, ...args); }
    };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#preset').selectOption('single');
  await page.locator('#alpha-number').fill('0'); await page.locator('#resolution').selectOption('80');
  await page.locator('#more-mesh-settings > summary').click();
}

test('MSET-style controls build and display the configured quad grid through the real GUI worker', async ({ page }) => {
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await openGrid(page);
  await page.locator('#grid-chord-exponent').fill('0.4');
  await page.locator('#grid-upper-streamlines').fill('10');
  await page.locator('#grid-lower-streamlines').fill('8');
  await page.locator('#grid-gap-streamlines').fill('6');
  await page.locator('#grid-inlet').selectOption('16'); await page.locator('#grid-outlet').selectOption('32');
  await expect(page.locator('#grid-le-ratio')).toBeDisabled();
  await page.locator('#grid-surface-spacing').selectOption('curvature');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await page.locator('#grid-elliptic').check();
  await expect(page.locator('#grid-le-ratio')).toBeEnabled();
  await page.locator('#grid-match-aspect').check(); await expect(page.locator('#grid-aspect-ratio')).toBeEnabled();
  await page.evaluate(() => { window.gridEvents = []; window.gridRequests = []; });
  await page.locator('#mesh-button').click();
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#mesh-status')).toContainText('SLOR smoothed');
  await expect(page.locator('#mesh-status')).toContainText('unsolved');
  await expect(page.locator('#cl')).toHaveText('—'); await expect(page.locator('#export-button')).toBeDisabled();
  const { events, requests } = await page.evaluate(() => ({ events: window.gridEvents, requests: window.gridRequests }));
  expect(events[0]).toMatchObject({ type: 'mesh', stage: 'initial' });
  expect(events.some(e => e.stage === 'smoothing')).toBe(true);
  expect(events.some(e => e.type === 'iteration')).toBe(false);
  expect(events.at(-1).type).toBe('mesh-ready');
  const mesh = events.filter(e => e.type === 'mesh').at(-1).mesh;
  expect(mesh.cells.every(c => c.length === 4)).toBe(true); expect(mesh.quality.valid).toBe(true);
  expect(mesh.initialization).toMatchObject({ tubes: [8, 9], flowSolved: false,
    normalSpacing: 'stagnation', gridSmoothing: { converged: true },
    gridSpacing: { surfaceChordExponent: .4, inlet: { intervals: 16 }, outlet: { intervals: 32 } } });
  expect(mesh.initialization.normalAllocation.groups.map(g => g.requestedTubes)).toEqual([7, 9]);
  expect(mesh.initialization.normalAllocation.groups.every(g => g.maxAdjacentRatio <= 3)).toBe(true);
  for (const d of mesh.initialization.surfaceDistributions) {
    expect(d.preMatching.actualLeadingSpacingRatio).toBeCloseTo(.2, 8); expect(d.preMatching.actualTrailingSpacingRatio).toBeCloseTo(.4, 8);
    expect(d.actualLeadingSpacingRatio).toBeGreaterThan(0); expect(d.actualTrailingSpacingRatio).toBeGreaterThan(0);
  }
  expect(JSON.stringify(requests)).toContain('"gridGapTubes":5');
  expect(JSON.stringify(requests)).toContain('"gridStagnationAspectRatio":2.5');
  await page.locator('#grid-le-ratio').fill('0.25');
  await expect(page.locator('#status')).toHaveText('Inputs changed');
  await expect(page.locator('#mesh-status')).toHaveText('Build mesh or run quad Euler to view the grid.');
  await page.locator('#grid-upper-streamlines').fill('3');
  await page.locator('#mesh-button').click();
  await expect(page.locator('#error-message')).toContainText('Upper streamlines');
  expect(errors).toEqual([]);
});

test('grid settings remain readable and operable at phone width', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await openGrid(page);
  await page.locator('#grid-surface-spacing').selectOption('curvature');
  await page.locator('#grid-match-aspect').check();
  const bounds = await page.locator('#streamtube-grid-conditions').evaluate(panel => ({
    fits: document.documentElement.scrollWidth <= innerWidth,
    controls: [...panel.querySelectorAll('input, select, summary')].map(el => {
      const r = el.getBoundingClientRect(); return { id: el.id, left: r.left, right: r.right, width: r.width };
    })
  }));
  expect(bounds.fits).toBe(true);
  for (const c of bounds.controls) { expect(c.left, c.id).toBeGreaterThanOrEqual(0); expect(c.right, c.id).toBeLessThanOrEqual(390); expect(c.width).toBeGreaterThan(0); }
  await page.locator('#streamtube-grid-conditions').screenshot({ path: '/tmp/mses-mset-grid-phone.png' });
  await page.locator('#more-mesh-settings > summary').click();
  await expect(page.locator('#grid-le-ratio')).toBeHidden();
  await expect(page.locator('#grid-elliptic')).toBeVisible();
  await page.locator('#more-mesh-settings > summary').click();
  await expect(page.locator('#grid-le-ratio')).toBeVisible();
});

test('main/flap curvature controls refine all four surfaces and build three independently sized regions', async ({ page }) => {
  await openGrid(page); await page.locator('#preset').selectOption('flap');
  await page.locator('#resolution').selectOption('160'); await page.locator('#alpha-number').fill('4');
  await page.locator('#grid-chord-exponent').fill('0.4');
  await page.locator('#grid-upper-streamlines').fill('10'); await page.locator('#grid-lower-streamlines').fill('8');
  await page.locator('#grid-gap-streamlines').fill('6');
  await page.locator('#grid-inlet').selectOption('16'); await page.locator('#grid-outlet').selectOption('32');
  await page.locator('#grid-surface-spacing').selectOption('curvature');
  await page.locator('#grid-elliptic').check();
  await page.evaluate(() => { window.gridEvents = []; }); await page.locator('#mesh-button').click();
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#mesh-status')).toContainText('local LE refinement');
  const events = await page.evaluate(() => window.gridEvents), mesh = events.filter(e => e.type === 'mesh').at(-1).mesh;
  expect(events.some(e => e.type === 'iteration')).toBe(false);
  expect(mesh.quality.valid).toBe(true); expect(mesh.initialization.tubes).toEqual([7, 5, 9]);
  expect(mesh.initialization.gridSmoothing.converged).toBe(true);
  expect(mesh.initialization.gridSpacing.requestedSurfaceIntervalsByElement.map(e => e.intervals).sort((a, b) => b - a)).toEqual([16, 10]);
  expect(mesh.initialization.surfaceDistributions).toHaveLength(4);
  for (const d of mesh.initialization.surfaceDistributions) {
    expect(d.preMatching.actualLeadingSpacingRatio).toBeCloseTo(.2, 8); expect(d.preMatching.actualTrailingSpacingRatio).toBeCloseTo(.4, 8);
    expect(d.actualLeadingSpacingRatio).toBeGreaterThan(0); expect(d.actualTrailingSpacingRatio).toBeGreaterThan(0);
  }
  expect(mesh.initialization.surfaceDistributions.some(d => d.artificialLeadingCurvature)).toBe(true);
  await expect(page.locator('#cl')).toHaveText('—');
});
