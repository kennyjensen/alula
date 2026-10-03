import { test, expect } from '@playwright/test';

// Measure rendered airfoil pixels rather than an implementation-only zoom
// variable, so a pan accidentally substituted for a pinch cannot pass.
const outline = page => page.locator('#geometry-canvas').evaluate(canvas => {
  const box = canvas.getBoundingClientRect(), scale = canvas.width / box.width;
  const { data } = canvas.getContext('2d').getImageData(0, 0, canvas.width, canvas.height);
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (let y = 0; y < canvas.height; y++) for (let x = 0; x < canvas.width; x++) {
    const i = 4 * (y * canvas.width + x), r = data[i], g = data[i + 1], b = data[i + 2];
    // Ignore faint antialiasing pixels and the mesh, labels and other elements.
    if (data[i + 3] > 200 && Math.abs(r - 137) < 8 && Math.abs(g - 229) < 8 && Math.abs(b - 208) < 8) {
      minX = Math.min(minX, x); maxX = Math.max(maxX, x); minY = Math.min(minY, y); maxY = Math.max(maxY, y);
    }
  }
  return { minX: minX / scale, maxX: maxX / scale, minY: minY / scale, maxY: maxY / scale, width: (maxX - minX) / scale };
});
const painted = page => page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));

test('phone pinch zoom anchors the airfoil and mesh, transitions to one-finger pan, and preserves page scrolling', async ({ browser, baseURL }) => {
  test.setTimeout(90000);
  const context = await browser.newContext({ baseURL, viewport: { width: 390, height: 844 }, deviceScaleFactor: 2, isMobile: true, hasTouch: true });
  const page = await context.newPage(), errors = []; page.on('pageerror', e => errors.push(e.message));
  try {
    await page.goto('/'); await expect(page.locator('#status')).toHaveText(/^Solved/);
    await page.locator('#flow-model').selectOption('streamtube-grid');
    await page.locator('#mesh-button').click(); await expect(page.locator('#status')).toHaveText('Mesh ready', { timeout: 60000 });
    const canvas = page.locator('#geometry-canvas'); await canvas.scrollIntoViewIfNeeded(); await painted(page);
    const initial = await outline(page), image = await canvas.evaluate(c => c.toDataURL());
    expect(initial.width).toBeGreaterThan(100);
    const box = await canvas.boundingBox(), center = { x: box.x + box.width / 2, y: box.y + box.height / 2 };
    const session = await context.newCDPSession(page);
    const point = (id, x, y) => ({ id, x, y, radiusX: 2, radiusY: 2, force: 1 });
    const touch = (type, touchPoints) => session.send('Input.dispatchTouchEvent', { type, touchPoints });
    let a = point(1, center.x - 40, center.y), b = point(2, center.x + 40, center.y);
    const scrollBefore = await page.evaluate(() => scrollY);
    await touch('touchStart', [a]); await touch('touchStart', [a, b]);
    for (let i = 1; i <= 6; i++) {
      const distance = 40 + 20 * i / 6;
      a = point(1, center.x + 8 * i / 6 - distance, center.y + 4 * i / 6);
      b = point(2, center.x + 8 * i / 6 + distance, center.y + 4 * i / 6);
      await touch('touchMove', [a, b]);
    }
    await painted(page); let zoomed = await outline(page);
    expect(zoomed.width / initial.width).toBeCloseTo(1.5, 1);
    expect(Math.abs(zoomed.minX - (box.width / 2 + 1.5 * (initial.minX - box.width / 2) + 8))).toBeLessThan(2);
    expect(await page.evaluate(() => visualViewport.scale)).toBe(1);
    expect(await page.evaluate(() => scrollY)).toBe(scrollBefore);
    // Bring the same fingers closer together to zoom back out.
    a = point(1, center.x + 8 - 50, center.y + 4); b = point(2, center.x + 8 + 50, center.y + 4);
    await touch('touchMove', [a, b]); await painted(page); zoomed = await outline(page);
    expect(zoomed.width / initial.width).toBeCloseTo(1.25, 1);
    expect(Math.abs(zoomed.minX - (box.width / 2 + 1.25 * (initial.minX - box.width / 2) + 8))).toBeLessThan(2);
    // Lifting one finger must not jump back to the pre-pinch view.
    await touch('touchEnd', [b]); await painted(page);
    expect(await outline(page)).toEqual(zoomed);
    a = { ...a, x: a.x + 18, y: a.y + 12 }; await touch('touchMove', [a]); await painted(page);
    const panned = await outline(page);
    expect(Math.abs(panned.minX - zoomed.minX - 18)).toBeLessThan(2);
    expect(Math.abs(panned.minY - zoomed.minY - 12)).toBeLessThan(2);
    await touch('touchCancel', []);
    await page.locator('#reset-view').click(); await painted(page);
    expect(await canvas.evaluate(c => c.toDataURL())).toBe(image);
    // Cancelled touches leave no stale pointer behind on the next gesture.
    a = point(3, center.x, center.y); await touch('touchStart', [a]);
    await touch('touchMove', [{ ...a, x: a.x + 10 }]); await touch('touchEnd', []); await painted(page);
    expect(Math.abs((await outline(page)).minX - initial.minX - 10)).toBeLessThan(2);
    const beforePageSwipe = await page.evaluate(() => scrollY);
    await touch('touchStart', [point(4, 4, 740)]);
    for (let i = 1; i <= 5; i++) await touch('touchMove', [point(4, 4, 740 - 35 * i)]);
    await touch('touchEnd', []);
    await expect.poll(() => page.evaluate(() => scrollY)).toBeGreaterThan(beforePageSwipe);
    expect(errors).toEqual([]);
  } finally { await context.close(); }
});
