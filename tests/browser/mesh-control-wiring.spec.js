import { test, expect } from '@playwright/test';

// Exercise the real form and request assembly without starting any numerical
// worker. The separate ISET browser regression checks the actual mesh path.
async function openControls(page) {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    window.meshControlRequests = [];
    window.Worker = class {
      postMessage(data) {
        window.meshControlRequests.push(structuredClone(data));
        queueMicrotask(() => {
          if (!this.stopped) this.onmessage?.({ data: { id: data.id, type: 'error',
            message: 'Control wiring test: numerical worker was not started.' } });
        });
      }
      terminate() { this.stopped = true; }
    };
  });
  await page.goto('/');
  await page.waitForFunction(() => window.meshControlRequests.length === 1);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.evaluate(() => { window.meshControlRequests = []; });
  return errors;
}

test('simple default mesh controls keep the same request input when only smoothing is toggled', async ({ page }) => {
  const errors = await openControls(page);
  await expect(page.locator('#preset')).toHaveValue('flap');
  await expect(page.locator('#resolution')).toHaveValue('160');
  await expect(page.locator('#alpha-number')).toHaveValue('4');
  await expect(page.locator('#streamtube-grid-conditions details')).toHaveCount(1);
  await expect(page.locator('#more-mesh-settings')).not.toHaveAttribute('open');
  await expect(page.locator('#more-mesh-settings > summary')).toHaveText('More mesh settings');
  for (const id of ['grid-intervals', 'grid-tubes', 'grid-inlet', 'grid-outlet', 'grid-elliptic'])
    await expect(page.locator(`#${id}`)).toBeVisible();
  for (const id of ['grid-upper-streamlines', 'grid-lower-streamlines', 'grid-gap-streamlines',
    'grid-chord-exponent', 'grid-surface-spacing', 'grid-match-aspect'])
    await expect(page.locator(`#${id}`)).toBeHidden();
  await expect(page.locator('#grid-crosslines')).toHaveCount(0);
  await expect(page.locator('#grid-smoothing-method')).toHaveCount(0);
  await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await expect(page.locator('#solve-button')).toBeEnabled();
  await page.locator('#mesh-button').click();
  await expect.poll(() => page.evaluate(() => window.meshControlRequests.length)).toBe(1);
  await page.locator('#grid-elliptic').uncheck();
  await expect(page.locator('#solve-button')).toBeEnabled();
  await page.locator('#mesh-button').click();
  await expect.poll(() => page.evaluate(() => window.meshControlRequests.length)).toBe(2);
  const [smooth, plain] = await page.evaluate(() => window.meshControlRequests);
  expect(plain).toMatchObject({ meshOnly: true, caseData: { flowModel: 'streamtube-grid',
    gridIntervals: 16, gridTubes: 7, gridCrosslinePlacement: 'potential', gridChordExponent: 0,
    gridSurfaceSpacing: 'automatic', gridSmoothingMethod: 'elliptic', gridEllipticSmoothing: false } });
  expect(smooth).toMatchObject({ meshOnly: true, caseData: { gridEllipticSmoothing: true } });
  expect(smooth.caseData).toEqual({ ...plain.caseData, gridEllipticSmoothing: true });
  expect(plain.caseData.gridInletIntervals).toBeUndefined();
  expect(plain.caseData.gridOutletIntervals).toBeUndefined();
  expect(plain.caseData.gridStagnationAspectRatio).toBeUndefined();
  expect(errors).toEqual([]);
});

test('one advanced mesh section wires sizing controls independently and fits on a phone', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  const errors = await openControls(page);
  await page.locator('#more-mesh-settings > summary').click();
  await page.locator('#grid-upper-streamlines').fill('10');
  await page.locator('#grid-lower-streamlines').fill('8');
  await page.locator('#grid-gap-streamlines').fill('6');
  await page.locator('#grid-chord-exponent').fill('0.4');
  await page.locator('#grid-inlet').selectOption('16');
  await page.locator('#grid-outlet').selectOption('32');
  await expect(page.locator('#grid-le-ratio')).toBeDisabled();
  await page.locator('#grid-surface-spacing').selectOption('curvature');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await expect(page.locator('#grid-le-ratio')).toBeEnabled();
  await page.locator('#grid-match-aspect').check();
  await expect(page.locator('#grid-aspect-ratio')).toBeEnabled();
  await page.locator('#mesh-button').click();
  await expect.poll(() => page.evaluate(() => window.meshControlRequests.length)).toBe(1);
  const request = await page.evaluate(() => window.meshControlRequests[0]);
  expect(request.caseData).toMatchObject({ gridCrosslinePlacement: 'potential', gridSmoothingMethod: 'elliptic',
    gridEllipticSmoothing: true, gridUpperTubes: 9, gridLowerTubes: 7, gridGapTubes: 5,
    gridChordExponent: .4, gridInletIntervals: 16, gridOutletIntervals: 32,
    gridSurfaceSpacing: 'curvature', gridStagnationAspectRatio: 2.5,
    gridCurvatureSpacing: { exponent: .5, leadingSpacingRatio: .2, trailingSpacingRatio: .4 } });
  const layout = await page.locator('#streamtube-grid-conditions').evaluate(panel => ({
    fits: document.documentElement.scrollWidth <= innerWidth,
    controls: [...panel.querySelectorAll('input, select, summary')].filter(el => el.getClientRects().length > 0).map(el => {
      const r = el.getBoundingClientRect(); return { id: el.id, left: r.left, right: r.right, width: r.width };
    }),
  }));
  expect(layout.fits).toBe(true);
  for (const control of layout.controls) {
    expect(control.left, control.id).toBeGreaterThanOrEqual(0);
    expect(control.right, control.id).toBeLessThanOrEqual(390);
    expect(control.width, control.id).toBeGreaterThan(0);
  }
  await page.locator('#streamtube-grid-conditions').screenshot({ path: '/tmp/mses-simple-mesh-controls-phone.png' });
  await page.locator('#more-mesh-settings > summary').click();
  await expect(page.locator('#grid-le-ratio')).toBeHidden();
  await expect(page.locator('#grid-elliptic')).toBeVisible();
  expect(errors).toEqual([]);
});
