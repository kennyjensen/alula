import { test, expect } from '@playwright/test';

async function open(page) {
  await page.addInitScript(() => { window.Worker = class { postMessage() {} terminate() {} }; });
  await page.goto('/');
}
async function swipe(page, locator, direction) {
  await expect.poll(() => page.locator('.app-shell').evaluate(node => Math.abs(node.scrollLeft - (document.querySelector('#visualization-tab').getAttribute('aria-selected') === 'true' ? node.clientWidth : 0)))).toBeLessThan(1);
  const box = await locator.boundingBox();
  const x = box.x + box.width * (direction === 'left' ? .8 : .2), y = box.y + box.height / 2;
  const session = await page.context().newCDPSession(page);
  await session.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x, y }] });
  for (let i = 1; i <= 6; i++) {
    await session.send('Input.dispatchTouchEvent', { type: 'touchMove', touchPoints: [{ x: x + (direction === 'left' ? -1 : 1) * box.width * .6 * i / 6, y }] });
    await page.waitForTimeout(20);
  }
  await session.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  await session.detach();
}
test('phone tabs support taps, keyboard and swipes without losing settings', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await open(page);
  await expect(page.locator('#visualization-tab')).toHaveAttribute('aria-selected', 'true');
  await expect(page.locator('#settings-pane')).toHaveAttribute('inert', '');
  await page.locator('#settings-tab').click();
  await expect(page.locator('#settings-pane')).not.toHaveAttribute('inert', '');
  await page.locator('#stop-button').click();
  await page.locator('#alpha-number').fill('6');
  await page.locator('#settings-pane').evaluate(node => node.scrollTop = 0);
  await expect.poll(() => page.locator('.app-shell').evaluate(node => node.scrollLeft)).toBe(0);
  await swipe(page, page.locator('.sidebar h1'), 'left');
  await expect(page.locator('#visualization-tab')).toHaveAttribute('aria-selected', 'true');
  await expect(page.locator('#overlay-alpha')).toHaveText('6.00');
  await expect(page.locator('#flow-speed-range, #flow-view-note')).toHaveCount(0);
  await swipe(page, page.locator('#geometry-canvas'), 'right');
  await expect(page.locator('#settings-tab')).toHaveAttribute('aria-selected', 'true');
  await expect(page.locator('#alpha-number')).toHaveValue('6');
  await page.locator('#settings-tab').focus();
  await page.keyboard.press('ArrowRight');
  await expect(page.locator('#visualization-tab')).toHaveAttribute('aria-selected', 'true');
  await swipe(page, page.locator('#pressure-canvas'), 'right');
  await expect(page.locator('#settings-tab')).toHaveAttribute('aria-selected', 'true');
  await page.locator('#visualization-tab').click();
  await expect.poll(() => page.locator('.app-shell').evaluate(node => Math.abs(node.scrollLeft - node.clientWidth))).toBeLessThan(1);
  const canvas = page.locator('#geometry-canvas');
  const box = await canvas.boundingBox();
  const image = await canvas.evaluate(node => node.toDataURL());
  const session = await page.context().newCDPSession(page);
  const points = distance => [
    { id: 1, x: box.x + box.width / 2 - distance, y: box.y + box.height / 2 },
    { id: 2, x: box.x + box.width / 2 + distance, y: box.y + box.height / 2 },
  ];
  await session.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: points(30) });
  await session.send('Input.dispatchTouchEvent', { type: 'touchMove', touchPoints: points(60) });
  await session.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  await session.detach();
  await expect.poll(() => canvas.evaluate(node => node.toDataURL())).not.toBe(image);
  await swipe(page, canvas, 'right');
  await expect(page.locator('#visualization-tab')).toHaveAttribute('aria-selected', 'true');
  expect(await page.locator('.pressure-panel').evaluate(node => {
    const box = node.getBoundingClientRect();
    return box.left >= 0 && box.right <= innerWidth && box.bottom <= innerHeight + 2;
  })).toBe(true);
});
test('desktop shows both columns and restores mobile selection on resize', async ({ page }) => {
  await open(page);
  await expect(page.locator('#mesh-button, #grid-audit-button')).toHaveCount(0);
  await expect(page.locator('#grid-audit')).toBeHidden();
  await expect(page.locator('#workspace-tabs')).toBeHidden();
  for (const id of ['settings-pane', 'visualization-pane']) {
    await expect(page.locator('#' + id)).not.toHaveAttribute('inert', '');
    expect(await page.locator('#' + id).evaluate(node => {
      const box = node.getBoundingClientRect(); return box.left >= 0 && box.right <= innerWidth;
    })).toBe(true);
  }
  await page.setViewportSize({ width: 390, height: 844 });
  await page.locator('#settings-tab').click();
  await page.setViewportSize({ width: 1440, height: 1100 });
  await expect(page.locator('#visualization-pane')).not.toHaveAttribute('inert', '');
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(page.locator('#settings-tab')).toHaveAttribute('aria-selected', 'true');
});

test('phone analysis action stays available in both tabs and returns to Settings on desktop', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.addInitScript(() => {
    window.analysisRequests = [];
    window.Worker = class {
      postMessage(request) { window.analysisRequests.push(request); }
      terminate() {}
    };
  });
  await page.goto('/');
  const run = page.locator('#solve-button');
  await expect(page.locator('body > #analysis-actions')).toHaveCount(1);
  await expect(run).toBeDisabled();
  await page.locator('#stop-button').click();
  await expect(run).toBeEnabled();
  await expect(run).toHaveText('Run');
  await expect(run.locator('.run-wind-icon')).toHaveCount(1);
  const before = await run.boundingBox();
  await run.click();
  await expect.poll(() => page.evaluate(() => window.analysisRequests.length)).toBe(2);
  await expect(run).toBeDisabled();
  await page.locator('#settings-tab').click();
  await page.locator('#stop-button').click();
  await page.locator('#settings-pane').evaluate(node => node.scrollTop = node.scrollHeight);
  expect(await run.boundingBox()).toEqual(before);
  await run.click();
  await expect.poll(() => page.evaluate(() => window.analysisRequests.length)).toBe(3);
  await page.setViewportSize({ width: 1440, height: 1100 });
  await expect(page.locator('#analysis-actions-slot > #analysis-actions')).toHaveCount(1);
  await expect(run).toHaveCSS('position', 'static');
});
