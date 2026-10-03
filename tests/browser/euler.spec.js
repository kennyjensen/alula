import { test, expect } from '@playwright/test';

test('Euler channel draws each accepted moving grid while solving and retains it on cancellation', async ({ page }) => {
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.addInitScript(() => {
    window.channelFrames = [];
    const OriginalWorker = window.Worker;
    window.Worker = class extends OriginalWorker {
      postMessage(data) {
        const receive = this.onmessage;
        this.onmessage = event => {
          receive(event);
          if (event.data.type !== 'mesh') return;
          window.channelFrames.push({ ...event.data,
            image: document.getElementById('field').toDataURL(),
            label: document.getElementById('field-label').textContent,
            busy: !document.getElementById('cancel').hidden });
          if (event.data.iteration.iteration === 2) document.getElementById('cancel').click();
        };
        return super.postMessage(data);
      }
    };
  });
  await page.goto('/euler.html'); await page.locator('#configuration').selectOption('channel');
  await page.locator('#solve').click();
  await expect(page.locator('#status')).toHaveText('Cancelled');
  const frames = await page.evaluate(() => window.channelFrames);
  expect(frames.map(f => f.iteration.iteration)).toEqual([0, 1, 2]);
  for (let i = 1; i < frames.length; i++) {
    expect(frames[i].busy).toBe(true);
    expect(frames[i].label).toContain(`Newton iteration ${i} · unconverged`);
    expect(frames[i].mesh.vertices).not.toEqual(frames[i - 1].mesh.vertices);
    expect(frames[i].image).not.toBe(frames[i - 1].image);
    expect(frames[i].iteration.maximumNodeMovement).toBeGreaterThan(0);
  }
  expect(await page.locator('#field').evaluate(c => c.toDataURL())).toBe(frames.at(-1).image);
  await expect(page.locator('#export')).toBeDisabled();
  await page.locator('#mach').fill('0.25');
  await expect(page.locator('#field-label')).toContainText('preview');
  await page.locator('.lab-advanced summary').click(); await page.locator('#iterations').fill('1');
  await page.evaluate(() => { window.channelFrames = []; });
  await page.locator('#solve').click();
  await expect(page.locator('#status')).toHaveText('Not converged');
  await expect(page.locator('#error')).toContainText('last accepted grid');
  await expect(page.locator('#field-label')).toContainText('Newton iteration 1 · unconverged');
  expect(await page.locator('#field').evaluate(c => c.toDataURL())).toBe(
    await page.evaluate(() => window.channelFrames.at(-1).image));
  await expect(page.locator('#export')).toBeDisabled();
  expect(errors).toEqual([]);
});

test('Euler laboratory solves interacting elements, exports topology, and clears stale states', async ({ page }) => {
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.goto('/euler.html');
  await expect(page.locator('#cells')).toHaveText('192');
  await expect(page.locator('#export')).toBeDisabled();
  await page.locator('#solve').click();
  await expect(page.locator('#status')).toHaveText('Equations converged', { timeout: 30000 });
  expect(Number(await page.locator('#residual').textContent())).toBeLessThan(1e-9);
  expect(Number(await page.locator('#mass').textContent())).toBeLessThan(1e-8);
  await page.locator('#quantity').selectOption('cp');
  await expect(page.locator('#field-label')).toContainText('Cp · cell averages');
  const downloaded = page.waitForEvent('download'); await page.locator('#export').click();
  const stream = await (await downloaded).createReadStream(); let text = ''; for await (const c of stream) text += c;
  const exported = JSON.parse(text);
  expect(exported.result.converged).toBe(true); expect(exported.controls.configuration).toBe('two');
  expect(exported.result.cl).toBeNull(); expect(exported.result.cd).toBeNull();
  expect(exported.result.mesh.cuts.filter(c => c.type === 'wake-cut').length).toBeGreaterThan(1);
  expect(exported.result.states).toHaveLength(192);
  await page.locator('#mach').fill('0.25');
  await expect(page.locator('#mass')).toHaveText('—'); await expect(page.locator('#export')).toBeDisabled();
  await expect(page.locator('#field-label')).toContainText('preview');
  expect(errors).toEqual([]);
});

test('Euler channel solves the moving grid and renders on mobile', async ({ page }) => {
  await page.goto('/euler.html'); await page.locator('#configuration').selectOption('channel');
  await expect(page.locator('#alpha')).toBeDisabled();
  await page.locator('#solve').click(); await expect(page.locator('#status')).toHaveText('Equations converged');
  expect(Number(await page.locator('#streamline').textContent())).toBeLessThan(1e-9);
  await expect(page.locator('#field-label')).toContainText('solved streamlines');
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  await expect(page.locator('#solve')).toBeEnabled();
});

test('Euler failure and cancellation cannot export unconverged or stale fields', async ({ page }) => {
  await page.goto('/euler.html');
  await page.locator('.lab-advanced summary').click(); await page.locator('#iterations').fill('1');
  await page.locator('#solve').click(); await expect(page.locator('#status')).toHaveText('Not converged');
  await expect(page.locator('#error')).toContainText('iteration limit');
  await expect(page.locator('#export')).toBeDisabled(); await expect(page.locator('#mass')).toHaveText('—');
  await page.locator('#iterations').fill('30'); await page.locator('#solve').click(); await page.locator('#cancel').click();
  await expect(page.locator('#status')).toHaveText('Cancelled'); await expect(page.locator('#export')).toBeDisabled();
  await page.locator('#configuration').selectOption('three');
  await expect(page.locator('#solve')).toBeEnabled(); await page.locator('#solve').click();
  await expect(page.locator('#status')).toHaveText('Equations converged', { timeout: 30000 });
  await page.locator('#panels').selectOption('40'); await page.locator('#rows').selectOption('4');
  await expect(page.locator('#error')).toContainText('600'); await expect(page.locator('#solve')).toBeDisabled();
  await expect(page.locator('#export')).toBeDisabled();
});
