import { test, expect } from '@playwright/test';

// Hold every Worker, including the initial panel request. These tests exercise
// the real form, message handlers, alert and copy action with zero solver runs.
async function open(page) {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    window.errorWorkers = [];
    window.errorRequests = [];
    window.errorCopies = [];
    window.clipboardDenied = false;
    Object.defineProperty(navigator, 'clipboard', { configurable: true, value: {
      async writeText(text) {
        if (window.clipboardDenied) throw new DOMException('Clipboard permission denied', 'NotAllowedError');
        window.errorCopies.push(text);
      },
    } });
    window.Worker = class {
      constructor() { window.errorWorkers.push(this); }
      postMessage(request) {
        this.request = structuredClone(request);
        window.errorRequests.push(this.request);
      }
      terminate() { this.stopped = true; }
    };
  });
  await page.goto('/');
  await page.waitForFunction(() => window.errorRequests.length === 1);
  await page.locator('#stop-button').click();
  await page.locator('#flow-model').selectOption('streamtube-bl');
  return errors;
}

async function start(page) {
  const before = await page.evaluate(() => window.errorRequests.length);
  await page.locator('#solve-button').click();
  await expect.poll(() => page.evaluate(() => window.errorRequests.length)).toBe(before + 1);
  return page.evaluate(() => ({ index: window.errorWorkers.length - 1, request: window.errorRequests.at(-1) }));
}

async function send(page, data, index = -1) {
  await page.evaluate(({ data, index }) => {
    const worker = window.errorWorkers.at(index);
    worker.onmessage({ data: { id: worker.request.id, ...data } });
  }, { data, index });
}

async function configureRae(page) {
  await page.locator('#preset').selectOption('rae2822-mses');
  await page.locator('#quad-mach').fill('.74');
  await page.locator('#alpha-number').fill('2.68');
  await page.locator('#quad-reynolds').fill('2700000');
  await page.locator('#quad-ncrit').fill('4');
  await page.locator('#quad-transition').selectOption('automatic');
  await page.locator('#quad-trip-upper').fill('1');
  await page.locator('#quad-trip-lower').fill('1');
  await page.locator('#grid-intervals').selectOption('64');
  await page.locator('#grid-tubes').selectOption('11');
  await page.locator('#grid-elliptic').check();
  await page.locator('#euler-ismom').selectOption('3');
  await page.locator('#grid-inlet').selectOption('32');
  await page.locator('#grid-outlet').selectOption('64');
  await page.locator('#more-mesh-settings > summary').click();
  await page.locator('#grid-upper-streamlines').fill('14');
  await page.locator('#grid-lower-streamlines').fill('12');
  await page.locator('#grid-chord-exponent').fill('.4');
  await page.locator('#grid-surface-spacing').selectOption('curvature');
  await page.locator('#grid-le-ratio').fill('.15');
  await page.locator('#grid-te-ratio').fill('.35');
  await page.locator('#grid-curvature-exponent').fill('.6');
  await page.locator('#grid-match-aspect').check();
  await page.locator('#grid-aspect-ratio').fill('3');
}

const rawError = 'Uniform stagnation-density startup failed: nonpositive interface pressure at interval 128, passage 1, tube 1.';
const failure = {
  type: 'error', message: rawError, code: 'streamtube-interface-pressure', stage: 'gas-initialization',
  actualMach: .1, requestedMach: .74,
  diagnostics: { cell: { i: 128, group: 1, tube: 1 }, interfacePressure: { lower: -.002, upper: 1.1 } },
};

test('error includes the submitted case, prioritizes actual failure Mach, and copies exact diagnostics with a phone fallback', async ({ page }, testInfo) => {
  const errors = await open(page);
  await configureRae(page);
  const submitted = await start(page);
  expect(submitted.request.caseData).toMatchObject({ mach: .74, alpha: 2.68, reynolds: 2700000, ncrit: 4,
    geometrySource: { id: 'rae2822-mses' }, gridIntervals: 64, gridTubes: 11, eulerIsmom: 3,
    gridEllipticSmoothing: true, transitionMode: 'automatic', materialTrips: [[1, 1]],
    gridInletIntervals: 32, gridOutletIntervals: 64, gridUpperTubes: 13, gridLowerTubes: 11,
    gridChordExponent: .4, gridSurfaceSpacing: 'curvature', gridStagnationAspectRatio: 3,
    gridCurvatureSpacing: { exponent: .6, leadingSpacingRatio: .15, trailingSpacingRatio: .35 } });
  expect(submitted.request.caseData.elements[0].points).toHaveLength(129);
  await send(page, { type: 'flow-stage', stage: 'euler', mach: .2, targetMach: .74 });
  await send(page, failure);
  await expect(page.locator('#error-message')).toBeVisible();
  await expect(page.locator('#error-description')).toHaveText(rawError);
  await expect(page.locator('#error-context')).toContainText(/RAE\s*2822/);
  await expect(page.locator('#error-context')).toContainText('Failing cell: i=128, group=1, tube=1');
  await expect(page.locator('#error-context')).toContainText('Interface pressure (solver units): lower -0.002, upper 1.1');
  for (const pattern of [/0\.74/, /0\.1/, /2\.68/, /Ncrit\s*[:=]?\s*4/i,
    /64 surface intervals.*11 tubes/, /SLOR: on/i, /ISMOM\s*[:=]?\s*3/i, /Inlet 32/]) {
    await expect(page.locator('#error-context')).toContainText(pattern);
  }
  await page.locator('#copy-error-debug').click();
  await expect.poll(() => page.evaluate(() => window.errorCopies.length)).toBe(1);
  await expect(page.locator('#error-copy-status')).toContainText(/copied/i);
  const copied = await page.evaluate(() => JSON.parse(window.errorCopies[0]));
  expect(copied.schemaVersion).toBe(1);
  expect(copied.input).toEqual(submitted.request.caseData);
  expect(copied.request.id).toBe(submitted.request.id);
  expect(copied.input.geometrySource.id).toBe('rae2822-mses');
  expect(copied.error).toMatchObject({ message: rawError, code: failure.code });
  expect(copied.actualMach).toBe(.1);
  expect(copied.requestedMach).toBe(.74);
  await expect(page.locator('#error-debug-text')).toBeHidden();

  await page.evaluate(() => { window.clipboardDenied = true; });
  await page.locator('#copy-error-debug').click();
  await expect(page.locator('#error-debug-text')).toBeVisible();
  expect(JSON.parse(await page.locator('#error-debug-text').inputValue())).toEqual(copied);
  await expect(page.locator('#error-copy-status')).toContainText(/copy|clipboard/i);
  expect(await page.evaluate(() => window.errorCopies.length)).toBe(1);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.locator('#error-description').scrollIntoViewIfNeeded();
  await expect(page.locator('#copy-error-debug')).toBeVisible();
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  expect(await page.locator('#error-context').evaluate(el => el.scrollWidth <= el.clientWidth)).toBe(true);
  await page.screenshot({ path: testInfo.outputPath('error-context-phone.png'), fullPage: true });
  await testInfo.attach('error-context-phone', { path: testInfo.outputPath('error-context-phone.png'), contentType: 'image/png' });
  await testInfo.attach('copied-diagnostics', { body: JSON.stringify(copied, null, 2), contentType: 'application/json' });
  expect(errors).toEqual([]);
});

test('late old-job errors are ignored and the next NLR failure cannot inherit the prior RAE case or Mach', async ({ page }) => {
  const errors = await open(page);
  await configureRae(page);
  const old = await start(page);
  await send(page, { type: 'flow-stage', stage: 'euler', mach: .2, targetMach: .74 });
  await send(page, failure);
  await expect(page.locator('#error-description')).toHaveText(rawError);
  await page.locator('#preset').selectOption('nlr7301');
  await page.locator('#quad-mach').fill('.185');
  await page.locator('#alpha-number').fill('6');
  await page.locator('#quad-reynolds').fill('2510000');
  await page.locator('#quad-ncrit').fill('9');
  await page.locator('#grid-intervals').selectOption('16');
  await page.locator('#grid-tubes').selectOption('7');
  await page.locator('#euler-ismom').selectOption('4');
  const current = await start(page);
  expect(current.request.id).not.toBe(old.request.id);
  await expect(page.locator('#error-message')).toBeHidden();
  await send(page, { ...failure, message: 'Late RAE failure must be ignored' }, old.index);
  await expect(page.locator('#error-message')).toBeHidden();
  await expect(page.locator('#status')).toHaveText('Solving');
  await send(page, { type: 'error', message: 'Current NLR boundary-layer failure', code: 'test-nlr-domain', stage: 'coupled' });
  await expect(page.locator('#error-description')).toHaveText('Current NLR boundary-layer failure');
  await expect(page.locator('#error-context')).toContainText(/NLR\s*7301/);
  await expect(page.locator('#error-context')).not.toContainText(/RAE\s*2822/);
  await expect(page.locator('#error-context')).not.toContainText(/0\.74|0\.10(?:0)?/);
  await page.locator('#copy-error-debug').click();
  await expect.poll(() => page.evaluate(() => window.errorCopies.length)).toBe(1);
  const copied = await page.evaluate(() => JSON.parse(window.errorCopies[0]));
  expect(copied.schemaVersion).toBe(1);
  expect(copied.input).toEqual(current.request.caseData);
  expect(copied.input.geometrySource.id).toBe('nlr7301');
  expect(copied.input.eulerIsmom).toBe(4);
  expect(copied.error).toMatchObject({ message: 'Current NLR boundary-layer failure', code: 'test-nlr-domain' });
  expect(copied.request.id).toBe(current.request.id);
  expect(copied.requestedMach).toBe(.185);
  expect(copied.actualMach).not.toBe(.1);
  expect(copied.actualMach).not.toBe(.2);
  expect(errors).toEqual([]);
});
