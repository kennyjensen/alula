import fs from 'node:fs';
import { test, expect } from '@playwright/test';

// Presentation only: every worker is held, so no mesh or flow solve occurs.
async function open(page, { preset } = {}) {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    window.ncritWorkers = [];
    window.Worker = class {
      constructor() { window.ncritWorkers.push(this); }
      postMessage(request) { this.request = structuredClone(request); }
      terminate() { this.stopped = true; }
    };
  });
  await page.goto('/');
  await page.waitForFunction(() => window.ncritWorkers[0]?.request);
  await page.locator('#stop-button').click();
  await page.locator('#flow-model').selectOption('streamtube-bl');
  if (preset) await page.locator('#preset').selectOption(preset);
  await page.locator('#quad-transition').selectOption('automatic');
  await page.locator('#quad-ncrit').fill('9');
  await page.locator('#solve-button').click();
  await page.waitForFunction(() => window.ncritWorkers.at(-1)?.request?.caseData?.quadBoundaryLayers);
  return errors;
}
async function send(page, message) {
  await page.evaluate(message => {
    const worker = window.ncritWorkers.at(-1);
    worker.onmessage({ data: { ...message, id: worker.request.id } });
  }, message);
}
function frame(ncrit, iteration = 4) {
  return { type: 'coefficients', iteration, mach: .2, targetMach: .2, actualNcrit: ncrit, targetNcrit: 9,
    stage: 'coupled-ncrit-startup', coefficients: { cl: .75, cd: .02, cm: -.1 },
    pressure: { referenceChord: 1, pressureKind: 'Synthetic condition-label fixture',
      elements: [{ name: 'Main', cp: [{ x: 0, y: 0, cp: 1 }, { x: .5, y: .1, cp: -.5 }, { x: 1, y: 0, cp: 0 }] }] } };
}

test('live Ncrit startup labels the actual frame and clears it at a new stage or cancellation', async ({ page }) => {
  const errors = await open(page);
  await send(page, { type: 'flow-stage', stage: 'coupled-ncrit-startup', actualNcrit: 7.5, targetNcrit: 9 });
  await expect(page.locator('#run-meta')).toContainText('Ncrit 7.5 → target 9');
  await send(page, frame(7.5));
  await expect(page.locator('#cp-status')).toContainText('Ncrit 7.5 → target 9');
  await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-ncrit', '7.5');
  await expect(page.locator('#coefficient-warning')).toContainText('Ncrit 7.5 → target 9');
  await expect(page.locator('#cl')).toHaveText('0.7500');
  await expect(page.locator('#cl').locator('..')).toHaveClass(/provisional/);
  await expect(page.locator('#quad-ncrit')).toHaveValue('9');
  await send(page, { type: 'iteration', iteration: { stage: 'coupled-ncrit-startup', iteration: 5,
    actualNcrit: 7.5, targetNcrit: 9, residual: .001, ncritContinuation: { phase: 'fine' } } });
  await expect(page.locator('#run-meta')).toContainText('Fine-grid Euler/BL startup');
  await expect(page.locator('#run-meta')).toContainText('Ncrit 7.5 → target 9');
  await send(page, { type: 'flow-stage', stage: 'coupled-ncrit-startup', actualNcrit: 8, targetNcrit: 9 });
  await expect(page.locator('#cp-status')).toBeHidden();
  await expect(page.locator('#pressure-canvas')).not.toHaveAttribute('data-ncrit');
  await expect(page.locator('#cl')).toHaveText('—');
  await send(page, frame(8));
  await expect(page.locator('#cp-status')).toContainText('Ncrit 8 → target 9');
  await page.locator('#stop-button').click();
  await expect(page.locator('#cp-status')).toBeHidden();
  await expect(page.locator('#pressure-canvas')).not.toHaveAttribute('data-ncrit');
  await expect(page.locator('#cl')).toHaveText('—');
  expect(errors).toEqual([]);
});

test('a retained intermediate root remains amber and explicitly unsuccessful at the submitted Ncrit', async ({ page }) => {
  const errors = await open(page);
  // Existing display fixture, metadata changed only for this held-worker UI
  // test. This is not a numerical result or evidence for Ncrit8 physics.
  const saved = JSON.parse(fs.readFileSync('docs/current-multielement-automatic-16x9-slor-browser.json'));
  const result = structuredClone(saved.result);
  result.conditions.ncrit = 8;
  result.checkpoint.restart.options.ncrit = 8;
  result.actualNcrit = 8; result.targetNcrit = 9; result.stateConverged = true;
  result.ncritContinuation = { actualNcrit: 8, targetNcrit: 9, reachedTarget: false };
  result.reason = 'Controlled target-stage failure';
  // Intentionally retain converged:true to test the UI's condition guard.
  await send(page, { type: 'result', result, elapsed: 1 });
  await expect(page.locator('#status')).toHaveText('Target Ncrit not reached');
  await expect(page.locator('#coefficient-warning')).toContainText('these values belong to Ncrit 8');
  await expect(page.locator('#error-description')).toContainText('retained converged state at Ncrit 8');
  await expect(page.locator('#error-context')).toContainText('Requested Ncrit 9 · Actual Ncrit 8');
  await expect(page.locator('#diagnostics')).toContainText('Actual / requested Ncrit');
  await expect(page.locator('#run-meta')).toContainText('Ncrit 8 → target 9');
  await expect(page.locator('#cp-status')).toContainText('Ncrit 8 → target 9');
  await expect(page.locator('#cl').locator('..')).toHaveClass(/provisional/);
  await expect(page.locator('#quad-ncrit')).toHaveValue('9');
  expect(errors).toEqual([]);
});

test('a converged target is not marked provisional because earlier Ncrit stages exist in its history', async ({ page }) => {
  const errors = await open(page);
  const saved = JSON.parse(fs.readFileSync('docs/current-multielement-automatic-16x9-slor-browser.json'));
  const result = structuredClone(saved.result);
  result.actualNcrit = result.targetNcrit = 9;
  result.ncritContinuation = { actualNcrit: 9, targetNcrit: 9, reachedTarget: true, stateConverged: true,
    attempts: [{ actualNcrit: 4, converged: true }, { actualNcrit: 7.5, converged: true }, { actualNcrit: 8, converged: true }] };
  await send(page, { type: 'result', result, elapsed: 1 });
  await expect(page.locator('#status')).toHaveText('Euler/BL converged · research');
  await expect(page.locator('#run-meta')).toContainText('Ncrit 9');
  await expect(page.locator('#coefficient-warning')).not.toContainText('was not reached');
  await expect(page.locator('#cl').locator('..')).not.toHaveClass(/provisional/);
  await expect(page.locator('#cp-status')).toBeHidden();
  await expect(page.locator('#quad-ncrit')).toHaveValue('9');
  expect(errors).toEqual([]);
});

test('the independently qualified default NLR Ncrit9 root displays its actual loads and refined mesh', async ({ page }) => {
  const errors = await open(page, { preset: 'nlr7301' });
  // Sealed numerical receipt, unchanged. The held Worker makes this a
  // presentation regression, not a cold-start convergence qualification.
  const base = 'docs/solver-reliability/gui-defaults-ncrit4/nlr7301-refined-to9';
  const audit = JSON.parse(fs.readFileSync(`${base}/independent-strict-root-audit.json`));
  const result = JSON.parse(fs.readFileSync(`${base}/displayed-result.json`));
  expect(audit.passed).toBe(true);
  expect(result.converged).toBe(true);
  expect(result.actualNcrit).toBe(9);
  expect(result.mesh.cells).toHaveLength(5610);
  await send(page, { type: 'result', result, elapsed: 1 });
  await expect(page.locator('#status')).toHaveText('Euler/BL converged · research');
  await expect(page.locator('#error-message')).toBeHidden();
  for (const id of ['cl', 'cd', 'cm']) {
    await expect(page.locator('#' + id)).toHaveText(result[id].toFixed(id === 'cd' ? 6 : 4));
    await expect(page.locator('#' + id).locator('..')).not.toHaveClass(/provisional/);
  }
  await expect(page.locator('#mesh-status')).toContainText('5,610');
  await expect(page.locator('#quad-ncrit')).toHaveValue('9');
  await expect(page.locator('#quad-reynolds')).toHaveValue('1000000');
  await expect(page.locator('#cp-status')).toBeHidden();
  expect(errors).toEqual([]);
});
