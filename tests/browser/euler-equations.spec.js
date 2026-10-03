import { test, expect } from '@playwright/test';

// All requests are held, including page startup. No native Worker or solve.
async function open(page) {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    window.equationWorkers = [];
    window.equationRequests = [];
    window.Worker = class {
      constructor() { window.equationWorkers.push(this); }
      postMessage(request) {
        this.request = structuredClone(request);
        window.equationRequests.push(this.request);
      }
      terminate() { this.stopped = true; }
    };
  });
  await page.goto('/');
  await page.waitForFunction(() => window.equationRequests.length === 1);
  await page.locator('#stop-button').click();
  return errors;
}

async function request(page, { stop = true } = {}) {
  const count = await page.evaluate(() => window.equationRequests.length);
  await page.locator('#solve-button').click();
  await expect.poll(() => page.evaluate(() => window.equationRequests.length)).toBe(count + 1);
  const value = await page.evaluate(() => window.equationRequests.at(-1));
  if (stop) await page.locator('#stop-button').click();
  return value;
}

for (const mode of ['streamtube-grid', 'streamtube-bl']) {
  test(`${mode} defaults to MSES-recommended ISMOM 4 and sends only the four numbered choices`, async ({ page }) => {
    const errors = await open(page);
    await page.locator('#flow-model').selectOption(mode);
    const selector = page.getByLabel('Euler formulation (ISMOM)');
    await expect(selector).toBeVisible();
    await expect(selector).toHaveValue('4');
    const choices = await selector.locator('option').evaluateAll(options => options.map(option => ({
      value: option.value, disabled: option.disabled,
    })));
    expect(choices).toEqual(['1', '2', '3', '4'].map(value => ({ value, disabled: false })));
    const baseline = await request(page);
    expect(baseline.caseData).toMatchObject({ flowModel: 'streamtube-grid', quadBoundaryLayers: mode === 'streamtube-bl', mach: .2 });
    expect(baseline.caseData.eulerIsmom).toBe(4);
    for (const eulerIsmom of [1, 2, 3, 4]) {
      await selector.selectOption(String(eulerIsmom));
      const selected = await request(page);
      expect(selected.caseData).toEqual({ ...baseline.caseData, eulerIsmom });
      expect(selected.task).toBe(baseline.task);
    }
    await selector.selectOption('4');
    expect((await request(page)).caseData).toEqual(baseline.caseData);
    expect(errors).toEqual([]);
  });
}

test('quad formulations are cached independently and panel requests omit the hidden selector', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#euler-ismom').selectOption('2');
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await expect(page.locator('#euler-ismom')).toHaveValue('4');
  await page.locator('#euler-ismom').selectOption('3');
  for (const mode of ['inviscid', 'coupled']) {
    await page.locator('#flow-model').selectOption(mode);
    await expect(page.locator('#euler-ismom')).toBeHidden();
    const panel = await request(page);
    expect(panel.caseData.flowModel).toBe(mode);
    expect(Object.hasOwn(panel.caseData, 'eulerIsmom')).toBe(false);
    expect(Object.hasOwn(panel.caseData, 'hybrid')).toBe(false);
  }
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await expect(page.locator('#euler-ismom')).toHaveValue('2');
  expect((await request(page)).caseData.eulerIsmom).toBe(2);
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await expect(page.locator('#euler-ismom')).toHaveValue('3');
  expect((await request(page)).caseData.eulerIsmom).toBe(3);
  expect(errors).toEqual([]);
});

test('changing equations cancels the old job, clears current values and sends the next selected mode', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await page.locator('#euler-ismom').selectOption('1');
  const previous = await request(page, { stop: false });
  await page.evaluate(() => {
    const worker = window.equationWorkers.at(-1);
    worker.onmessage({ data: { id: worker.request.id, type: 'coefficients', iteration: 2, mach: .2, targetMach: .2,
      coefficients: { cl: .4, cd: .02, cm: -.03 } } });
  });
  await expect(page.locator('#cl')).toHaveText('0.4000');
  await page.locator('#euler-ismom').selectOption('4');
  await expect(page.locator('#status')).toHaveText('Inputs changed');
  await expect(page.locator('#stop-button')).toBeHidden();
  await expect(page.locator('#cl')).toHaveText('—');
  await expect(page.locator('#export-button')).toBeDisabled();
  expect(await page.evaluate(() => window.equationWorkers.at(-1).stopped)).toBe(true);
  const next = await request(page, { stop: false });
  expect(next.id).not.toBe(previous.id);
  expect(next.caseData.eulerIsmom).toBe(4);
  expect(next.caseData).toEqual({ ...previous.caseData, eulerIsmom: 4 });
  await page.evaluate(() => {
    const obsolete = window.equationWorkers.at(-2);
    obsolete.onmessage({ data: { id: obsolete.request.id, type: 'error', message: 'Obsolete equation solve' } });
  });
  await expect(page.locator('#status')).toHaveText('Solving');
  await expect(page.locator('#error-message')).toBeHidden();
  await page.locator('#stop-button').click();
  await expect(page.locator('#euler-ismom')).toHaveValue('4');
  expect(errors).toEqual([]);
});

test('result diagnostics use the returned formulation, including checkpoint fallback, rather than the selected value', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-bl');
  // Match the saved result's target before checking its formulation display.
  // A fixture at Ncrit 4 cannot certify the GUI's default Ncrit 9 request.
  await page.locator('#quad-ncrit').fill('4');
  await page.locator('#grid-intervals').selectOption('32');
  await page.locator('#euler-ismom').selectOption('4');
  for (const source of ['settings', 'checkpoint', 'default']) {
    await request(page, { stop: false });
    await page.evaluate(async source => {
      const response = await fetch('/docs/rae2822/real-browser-mses-coarse-cold32/state.json');
      if (!response.ok) throw new Error('Missing saved display fixture.');
      const { result } = await response.json();
      // Presentation-only metadata variants of the saved display result;
      // these do not represent new numerical qualifications of the modes.
      result.solverSettings = { ...result.solverSettings, streamwiseMode: 'hybrid' };
      if (source === 'settings') result.solverSettings.hybrid = { ismom: 3 };
      else {
        delete result.solverSettings.hybrid;
        result.checkpoint.restart.input.hybrid = source === 'checkpoint' ? { ismom: 2 } : {};
        result.checkpoint.restart.input.streamwiseMode = 'hybrid';
      }
      const worker = window.equationWorkers.at(-1);
      worker.onmessage({ data: { id: worker.request.id, type: 'result', result, elapsed: 0 } });
    }, source);
    await expect(page.locator('#status')).toHaveText('Euler/BL converged · research');
    await expect(page.locator('#diagnostics')).toContainText({
      settings: 'ISMOM 3 · Leading-edge hybrid', checkpoint: 'ISMOM 2 · Entropy',
      default: 'Automatic hybrid · ISMOM 4-style',
    }[source]);
    await expect(page.locator('#diagnostics')).not.toContainText('ISMOM 4 · Automatic hybrid');
    await expect(page.locator('#euler-ismom')).toHaveValue('4');
  }
  expect(errors).toEqual([]);
});
