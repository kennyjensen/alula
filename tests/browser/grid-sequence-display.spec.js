import fs from 'node:fs';
import { test, expect } from '@playwright/test';

// All workers are held. Saved states test presentation; no flow is solved.
async function open(page, intervals = '128') {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    window.gridWorkers = [];
    window.Worker = class {
      constructor() { window.gridWorkers.push(this); }
      postMessage(request) { this.request = structuredClone(request); }
      terminate() { this.stopped = true; }
    };
  });
  await page.goto('/');
  await page.waitForFunction(() => window.gridWorkers[0]?.request);
  await page.locator('#stop-button').click();
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await page.locator('#grid-intervals').selectOption(intervals);
  return errors;
}
async function solve(page) {
  await page.locator('#solve-button').click();
  await page.waitForFunction(() => window.gridWorkers.at(-1)?.request?.caseData?.quadBoundaryLayers);
}
async function send(page, message) {
  await page.evaluate(message => {
    const worker = window.gridWorkers.at(-1);
    worker.onmessage({ data: { ...message, id: worker.request.id } });
  }, message);
}
const saved = () => JSON.parse(fs.readFileSync('docs/current-multielement-automatic-16x9-slor-browser.json')).result;

test('grid controls count intervals and streamlines explicitly; actual surfaces use original element names', async ({ page }) => {
  const errors = await open(page);
  await expect(page.locator('#grid-intervals option:checked')).toHaveText('128 intervals · 129 points');
  await page.locator('#grid-tubes').selectOption('7');
  await expect(page.locator('#grid-tubes option:checked')).toHaveText('7 tubes · 8 streamlines');
  await page.locator('#more-mesh-settings summary').click();
  await expect(page.locator('#grid-upper-streamlines')).toHaveAttribute('placeholder', 'Auto · 8');
  await page.locator('#grid-upper-streamlines').fill('19');
  await page.locator('#grid-lower-streamlines').fill('15');
  await page.locator('#grid-gap-streamlines').fill('10');
  await solve(page);
  const input = await page.evaluate(() => window.gridWorkers.at(-1).request.caseData);
  expect(input.gridIntervals).toBe(128); expect(input.gridTubes).toBe(7);
  expect(input.gridUpperTubes).toBe(18); expect(input.gridLowerTubes).toBe(14); expect(input.gridGapTubes).toBe(9);
  const mesh = structuredClone(saved().mesh);
  mesh.initialization.tubes = [10, 13, 10];
  // Solver order differs from original element order, as for the NLR assembly.
  mesh.initialization.surfaceIntervals = [{ element: 1, intervals: 32 }, { element: 0, intervals: 16 }];
  await send(page, { type: 'mesh', mesh, stage: 'solving' });
  await expect(page.locator('#mesh-status')).toContainText('Actual streamlines, bottom → top passages [11, 14, 11]');
  await expect(page.locator('#mesh-status')).toContainText('Surface points/side: Main element 17; Flap 33');
  await page.locator('#stop-button').click();
  expect(errors).toEqual([]);
});

test('128-surface grid can request 24 base tubes per region', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#grid-tubes').selectOption('24');
  await expect(page.locator('#grid-tubes option:checked')).toHaveText('24 tubes · 25 streamlines');
  await expect(page.locator('#grid-upper-streamlines')).toHaveAttribute('placeholder', 'Auto · 25');
  await solve(page);
  const input = await page.evaluate(() => window.gridWorkers.at(-1).request.caseData);
  expect(input).toMatchObject({ gridIntervals: 128, gridTubes: 24 });
  expect(errors).toEqual([]);
});

test('intermediate iteration and Cp labels follow the current grid and clear on a level change', async ({ page }) => {
  const errors = await open(page); await solve(page);
  await send(page, { type: 'flow-stage', stage: 'coupled-grid-refinement', gridLevel: 32, requestedGridIntervals: 128 });
  await expect(page.locator('#run-meta')).toContainText('Preparing Euler/BL grid 32 → 128 intervals/side');
  await send(page, { type: 'iteration', iteration: { stage: 'coupled-grid-refinement', gridLevel: 32,
    requestedGridIntervals: 128, iteration: 3, residual: 1e-3 } });
  await expect(page.locator('#run-meta')).toContainText('Euler/BL grid 32 → 128 intervals/side · iteration 3');
  await expect(page.locator('#run-meta')).not.toContainText('Wake 1');
  await send(page, { type: 'coefficients', stage: 'coupled-grid-refinement', gridLevel: 32, requestedGridIntervals: 128,
    iteration: 3, mach: .2, targetMach: .2, coefficients: { cl: .75, cd: .02, cm: -.1 },
    pressure: { referenceChord: 1, pressureKind: 'Synthetic grid-label fixture',
      elements: [{ name: 'Main', cp: [{ x: 0, y: 0, cp: 1 }, { x: .5, y: .1, cp: -.5 }, { x: 1, y: 0, cp: 0 }] }] } });
  await expect(page.locator('#cp-status')).toContainText('Grid 32 → target 128 intervals/side');
  await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-grid-intervals', '32');
  await send(page, { type: 'flow-stage', stage: 'coupled-grid-refinement', gridLevel: 16,
    requestedGridIntervals: 128, retained: true });
  await expect(page.locator('#run-meta')).toContainText('Retained Euler/BL grid 16 → 128 intervals/side');
  await expect(page.locator('#cp-status')).toBeHidden();
  await expect(page.locator('#pressure-canvas')).not.toHaveAttribute('data-grid-intervals');
  await expect(page.locator('#cl')).toHaveText('—');
  await page.locator('#stop-button').click();
  expect(errors).toEqual([]);
});

test('a retained root with interval-only refinement metadata stays amber and cannot enable refinement', async ({ page }) => {
  const errors = await open(page); await solve(page);
  const result = structuredClone(saved());
  // Metadata-only fixture with a deliberately stale converged:true flag.
  result.gridSequence = { kind: 'coarse-to-fine', reachedTarget: false,
    actualGridIntervals: 16, requestedGridIntervals: 128, failure: { reason: 'Controlled level failure' } };
  result.automaticRefinement = result.gridSequence;
  result.stateConverged = true; result.reason = 'Controlled level failure';
  await send(page, { type: 'result', result, elapsed: 1 });
  await expect(page.locator('#status')).toHaveText('Target grid not reached');
  await expect(page.locator('#coefficient-warning')).toContainText('these values belong to the retained 16-interval grid');
  await expect(page.locator('#error-description')).toContainText('retained converged state on the 16-interval grid');
  await expect(page.locator('#diagnostics')).toContainText('Actual / requested surface intervals per side');
  await expect(page.locator('#diagnostics')).toContainText('16 / 128');
  await expect(page.locator('#run-meta')).toContainText('Grid 16 → target 128 intervals/side');
  await expect(page.locator('#cp-status')).toContainText('Grid 16 → target 128 intervals/side');
  await expect(page.locator('#cl').locator('..')).toHaveClass(/provisional/);
  await expect(page.locator('#refine-coupled-button')).toBeDisabled();
  await expect(page.locator('#grid-intervals')).toHaveValue('128');
  expect(errors).toEqual([]);
});

test('a target grid root is not marked provisional because coarser levels exist in its history', async ({ page }) => {
  const errors = await open(page, '16'); await solve(page);
  const result = structuredClone(saved());
  result.gridSequence = { kind: 'coarse-to-fine', reachedTarget: true, actualGridIntervals: 16,
    requestedGridIntervals: 16, levels: [{ intervals: 8, converged: true }] };
  result.automaticRefinement = result.gridSequence;
  await send(page, { type: 'result', result, elapsed: 1 });
  await expect(page.locator('#status')).toHaveText('Euler/BL converged · research');
  await expect(page.locator('#coefficient-warning')).not.toContainText('was not reached');
  await expect(page.locator('#cl').locator('..')).not.toHaveClass(/provisional/);
  await expect(page.locator('#cp-status')).toBeHidden();
  expect(errors).toEqual([]);
});
