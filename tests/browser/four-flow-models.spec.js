import { test, expect } from '@playwright/test';
import { parseCoordinates } from '../../src/geometry/parse.js';
import { getBenchmarkAirfoil } from '../../src/geometry/benchmark-airfoils.js';
import fs from 'node:fs';

test.afterEach(async ({ page }, info) => {
  const capture = await page.evaluate(() => ({ requests: window.heldFlowRequests,
    workerStubs: window.heldFlowWorkers?.length, nativeWorkersCreated: 0, forwardedRequests: 0 })).catch(() => null);
  const path = `docs/four-solver-ui/browser-case-${info.titlePath.at(-1).split(' ')[0]}${process.env.FOUR_UI_CAPTURE_SUFFIX ?? ''}.json`;
  fs.writeFileSync(path, JSON.stringify({ title: info.title, status: info.status, capture }, null, 2) + '\n');
});

async function open(page) {
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.addInitScript(() => {
    window.heldFlowRequests = []; window.heldFlowWorkers = [];
    // No native Worker is created and no request is forwarded.
    window.Worker = class {
      constructor() { window.heldFlowWorkers.push(this); }
      postMessage(data) { this.request = structuredClone(data); window.heldFlowRequests.push(this.request); }
      terminate() { this.stopped = true; }
    };
  });
  await page.goto('/');
  await page.waitForFunction(() => window.heldFlowRequests.length === 1);
  await page.locator('#stop-button').click();
  return errors;
}
async function request(page) {
  const count = await page.evaluate(() => window.heldFlowRequests.length);
  await page.locator('#solve-button').click();
  await expect.poll(() => page.evaluate(() => window.heldFlowRequests.length)).toBe(count + 1);
  const result = await page.evaluate(() => window.heldFlowRequests.at(-1));
  await page.locator('#stop-button').click();
  return result;
}
async function download(page) {
  const promise = page.waitForEvent('download');
  await page.locator('#coordinate-export-button').click();
  const result = await promise, stream = await result.createReadStream();
  let text = ''; for await (const chunk of stream) text += chunk;
  return { name: result.suggestedFilename(), text };
}

test('four explicit modes route only their visible physics and preserve independent settings', async ({ page }) => {
  const errors = await open(page);
  expect(await page.locator('#flow-model option').evaluateAll(nodes => nodes.map(n => n.value)))
    .toEqual(['inviscid', 'coupled', 'streamtube-grid', 'streamtube-bl']);
  for (const id of ['quad-viscous', 'subcritical-conditions', 'outer-mesh', 'wake-extent', 'mach'])
    await expect(page.locator('#' + id)).toHaveCount(0);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#quad-mach').fill('.28');
  await page.locator('#grid-intervals').selectOption('32');
  await page.locator('#grid-elliptic').uncheck();
  let r = await request(page);
  expect(r.caseData).toMatchObject({ flowModel: 'streamtube-grid', quadBoundaryLayers: false, mach: .28, gridIntervals: 32, gridEllipticSmoothing: false });
  for (const key of ['reynolds', 'ncrit', 'materialTrips', 'transitionMode']) expect(r.caseData[key]).toBeUndefined();
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await expect(page.locator('#quad-mach')).toHaveValue('0.2');
  await expect(page.locator('#grid-intervals')).toHaveValue('16');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await expect(page.locator('#quad-viscous-conditions')).toBeVisible();
  await expect(page.locator('#viscous-conditions')).toBeHidden();
  await page.locator('#quad-mach').fill('.25');
  await page.locator('#quad-transition').selectOption('automatic');
  await page.locator('#quad-reynolds').fill('2700000');
  await page.locator('#quad-ncrit').fill('4');
  r = await request(page);
  expect(r.caseData).toMatchObject({ flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .25,
    gridIntervals: 16, reynolds: 2700000, ncrit: 4, transitionMode: 'automatic', materialTrips: [[1, 1], [1, 1]] });
  expect(r.task).toBeUndefined();
  await page.locator('#flow-model').selectOption('coupled');
  await page.locator('#reynolds').fill('3000000');
  await page.locator('#element-list details').first().locator('summary').click();
  await page.locator('[data-key="tripUpper"]').first().fill('.3');
  r = await request(page);
  expect(r.caseData).toMatchObject({ flowModel: 'coupled', mach: 0, reynolds: 3000000, ncrit: 9, trips: [1, 1] });
  expect(r.caseData.elements[0].trips).toEqual([.3, 1]);
  for (const key of ['quadBoundaryLayers', 'materialTrips', 'gridIntervals']) expect(r.caseData[key]).toBeUndefined();
  await page.locator('#flow-model').selectOption('inviscid');
  r = await request(page);
  expect(r.caseData).toMatchObject({ flowModel: 'inviscid', mach: 0 });
  for (const key of ['quadBoundaryLayers', 'materialTrips', 'reynolds', 'gridIntervals']) expect(r.caseData[key]).toBeUndefined();
  expect(r.caseData.elements.every(e => e.trips === undefined)).toBe(true);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await expect(page.locator('#quad-mach')).toHaveValue('.28');
  await expect(page.locator('#grid-intervals')).toHaveValue('32');
  await expect(page.locator('#grid-elliptic')).not.toBeChecked();
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await expect(page.locator('#quad-mach')).toHaveValue('.25');
  await expect(page.locator('#quad-transition')).toHaveValue('automatic');
  await expect(page.locator('#quad-reynolds')).toHaveValue('2700000');
  await page.locator('#flow-model').selectOption('coupled');
  await expect(page.locator('#reynolds')).toHaveValue('3000000');
  expect(errors).toEqual([]);
});

test('NLR permits each explicit mode and a mode change rejects late results', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#preset').selectOption('nlr7301');
  let r = await request(page);
  const exact = getBenchmarkAirfoil('nlr7301');
  expect(r.caseData.elements.map(e => ({ points: e.points, trailingEdge: e.trailingEdge })))
    .toEqual(exact.elements.map(e => ({ points: e.points, trailingEdge: e.trailingEdge })));
  for (const mode of ['coupled', 'streamtube-bl']) {
    await page.locator('#flow-model').selectOption(mode);
    await expect(page.locator('#solve-button')).toBeEnabled();
    await expect(page.locator('#preset-note')).toContainText('All four solver modes are available');
  }
  await page.locator('#flow-model').selectOption('coupled');
  await expect(page.locator('#solve-button')).toBeEnabled();
  await page.locator('#flow-model').selectOption('inviscid');
  await expect(page.locator('#solve-button')).toBeEnabled();
  await page.locator('#solve-button').click();
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.evaluate(() => {
    const w = window.heldFlowWorkers.at(-1);
    w.onmessage?.({ data: { id: w.request.id, type: 'error', message: 'Late obsolete panel result' } });
  });
  await expect(page.locator('#status')).toHaveText('Inputs changed');
  await expect(page.locator('#error-message')).toBeHidden();
  await expect(page.locator('#export-button')).toBeDisabled();
  expect(errors).toEqual([]);
});

const openSurface = [{ x: 1, y: .02 }, { x: .75, y: .06 }, { x: .5, y: .08 }, { x: .25, y: .07 },
  { x: 0, y: 0 }, { x: .25, y: -.04 }, { x: .5, y: -.05 }, { x: .75, y: -.03 }, { x: 1, y: -.01 }];
const lines = points => points.map(p => `${p.x} ${p.y}`).join('\n');

test('standard open DAT imports retain both endpoints and source export; new MSES domain is explicit', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#import-file').setInputFiles({ name: 'open.dat', mimeType: 'text/plain', buffer: Buffer.from('Open finite airfoil\n' + lines(openSurface) + '\n') });
  let r = await request(page);
  expect(r.caseData.elements[0]).toMatchObject({ sourcePoints: openSurface,
    points: [...openSurface, openSurface[0]], trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 8 } });
  await page.locator('#coordinate-export > summary').click();
  let saved = JSON.parse((await download(page)).text);
  expect(saved.elements[0].sourcePoints).toEqual(openSurface);
  expect(saved.elements[0].trailingEdge).toEqual(r.caseData.elements[0].trailingEdge);
  await page.locator('#coordinate-format').selectOption('source');
  saved = parseCoordinates((await download(page)).text);
  expect(saved.elements[0].points).toEqual(openSurface);
  await page.locator('#coordinate-format').selectOption('xfoil-labeled');
  await page.locator('#coordinate-points').selectOption('wetted');
  saved = parseCoordinates((await download(page)).text);
  expect(saved.elements[0].points).toEqual(openSurface);
  await page.locator('#coordinate-format').selectOption('mses');
  await page.locator('#coordinate-export-button').click();
  await expect(page.locator('#error-message')).toContainText('all four MSES domain bounds');
  for (const [id, value] of [['x-min', '-2'], ['x-max', '3'], ['y-min', '-1'], ['y-max', '1']])
    await page.locator('#coordinate-' + id).fill(value);
  saved = parseCoordinates((await download(page)).text);
  expect(saved.domain).toEqual({ xMin: -2, xMax: 3, yMin: -1, yMax: 1 });
  expect(saved.elements[0].points).toEqual(openSurface);
  expect(errors).toEqual([]);
});

test('MSES common coordinates and legacy header survive import/export; phone has no overflow', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  const second = openSurface.map(p => ({ x: 1.1 + .3 * p.x, y: -.05 + .3 * p.y }));
  const header = [-3, 6, 3, 1.5, 1.2];
  const text = 'Measured pair\n' + header.join(' ') + '\n' + lines(openSurface) + '\n999 999\n' + lines(second) + '\n';
  await page.locator('#import-file').setInputFiles({ name: 'blade.pair', mimeType: 'text/plain', buffer: Buffer.from(text) });
  await expect(page.locator('#coordinate-import-note')).toContainText('does not override');
  await page.locator('#coordinate-export > summary').click();
  await page.locator('#coordinate-format').selectOption('source');
  const saved = parseCoordinates((await download(page)).text);
  expect(saved.name).toBe('Measured pair');
  expect(saved.header).toEqual({ kind: 'legacy-ises-mses', values: header });
  expect(saved.domain).toBeNull();
  expect(saved.elements.map(e => e.points)).toEqual([openSurface, second]);
  const r = await request(page);
  expect(r.caseData.elements.map(e => e.sourcePoints)).toEqual([openSurface, second]);
  expect(r.caseData.gridIntervals).toBe(16);
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);
  await page.screenshot({ path: `docs/four-solver-ui/browser${process.env.FOUR_UI_CAPTURE_SUFFIX ?? ''}.png`, fullPage: true });
  expect(errors).toEqual([]);
});
