// SPDX-License-Identifier: GPL-2.0-or-later
// Opt for one preset with MSES_GUI_DEFAULT_CASE=flap; otherwise this file checks
// every visible preset. Run with --workers=1. This performs real cold solves.
import fs from 'node:fs';
import path from 'node:path';
import { test, expect } from '@playwright/test';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';
import { guiDefaultCoupledCases, buildGuiDefaultCoupledCase, assessGuiDefaultCoupledResult } from '../../scripts/validation/gui-default-coupled-cases.js';
import { compareGuiInput } from '../../scripts/validation/gui-input-comparison.js';

const selected = process.env.MSES_GUI_DEFAULT_CASE;
if (selected && !guiDefaultCoupledCases.some(c => c.id === selected)) throw new Error(`Unknown MSES_GUI_DEFAULT_CASE: ${selected}`);
const outputRoot = process.env.MSES_GUI_DEFAULT_OUTPUT ?? 'docs/solver-reliability/gui-defaults-browser';
// Test watchdog only: Playwright uses a monotonic clock; the public worker
// retains its numerical budgets. Date elapsed time can also include host sleep
// or clock corrections, so retain both clock scopes in new receipts.
// Automatic physical-parameter continuation can take longer than one stage.
const timeoutMs = Number(process.env.MSES_GUI_DEFAULT_TIMEOUT_MS ?? 360000);
if (!Number.isInteger(timeoutMs) || timeoutMs < 1000 || timeoutMs > 900000)
  throw new Error('MSES_GUI_DEFAULT_TIMEOUT_MS must be an integer from 1000 to 900000.');
const serial = value => JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v);

for (const spec of guiDefaultCoupledCases.filter(c => !selected || c.id === selected))
test(`${spec.label}: untouched GUI coupled defaults converge through the real solver worker`, async ({ page }) => {
  test.setTimeout(timeoutMs + 90000);
  const output = path.join(outputRoot, spec.id);
  expect(fs.existsSync(output), 'Preserve earlier receipts; supply a NEW MSES_GUI_DEFAULT_OUTPUT directory.').toBe(false);
  fs.mkdirSync(output, { recursive: true });
  const sourceHashes = numericalSourceHashes(['tests/browser/gui-default-coupled-acceptance.spec.js',
    'scripts/validation/gui-default-coupled-cases.js', 'scripts/validation/gui-input-comparison.js', 'scripts/validation/solver-reliability-cases.js',
    'scripts/validation/solver-reliability-acceptance.js', 'src/worker/solver.js', 'src/ui/app.js',
    'src/ui/quad-coupled-result.js', 'src/ui/quad-coupled-ncrit.js', 'src/ui/quad-smoothing-status.js', 'index.html']);
  const expected = buildGuiDefaultCoupledCase(spec.id), started = Date.now(), startedMonotonic = performance.now(), errors = [];
  const write = (name, value) => fs.writeFileSync(path.join(output, name + '.json'), serial(value) + '\n');
  write('start', { spec, expectedInput: expected, sourceHashes, timeoutMs,
    scope: 'Actual DOM defaults and real Worker cold solve. No checkpoint injection, prepared geometry, numerical control override or substitute equation mode. Numerical acceptance does not establish aerodynamic accuracy.' });
  page.on('pageerror', e => errors.push(e.message));
  await page.exposeFunction('recordGuiDefaultProgress', event => {
    fs.appendFileSync(path.join(output, 'progress.jsonl'), serial({ seconds: (Date.now() - started) / 1000, monotonicSeconds: (performance.now() - startedMonotonic) / 1000, ...event }) + '\n');
  });
  await page.addInitScript(() => {
    const RealWorker = window.Worker;
    window.guiDefaultCapture = { request: null, terminal: null, events: [] };
    window.Worker = class extends RealWorker {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => {
          if (!this.guiDefaultRequest?.caseData?.quadBoundaryLayers) return;
          const capture = window.guiDefaultCapture;
          if (data.type === 'result' || data.type === 'error') capture.terminal = data;
          const h = data.iteration;
          capture.events.push({ type: data.type, stage: data.stage, message: data.message,
            actualNcrit: data.actualNcrit, targetNcrit: data.targetNcrit,
            ...(data.type === 'iteration' ? { iteration: { stage: h.stage, iteration: h.iteration, residual: h.residual,
              euler: h.euler, boundaryLayer: h.boundaryLayer, edgeMatching: h.edgeMatching, step: h.step, backtracks: h.backtracks,
              limiter: h.limiter, viscousLimiter: h.viscousLimiter, changes: h.changes,
              actualNcrit: h.actualNcrit, targetNcrit: h.targetNcrit } } : {}),
            ...(data.type === 'mesh' ? { mesh: { cells: data.mesh.cells.length, quality: data.mesh.quality, iteration: data.mesh.iteration } } : {}),
            ...(data.type === 'pressure' ? { pressure: { iteration: data.iteration, stage: data.stage, mach: data.mach } } : {}) });
          if (data.type === 'flow-stage' || data.type === 'iteration' && h.iteration % 5 === 0)
            window.recordGuiDefaultProgress(capture.events.at(-1)).catch(error => { capture.observationError = error.message; });
        });
      }
      postMessage(data, ...rest) {
        this.guiDefaultRequest = data;
        if (data.caseData?.quadBoundaryLayers) window.guiDefaultCapture.request = structuredClone(data);
        return super.postMessage(data, ...rest);
      }
    };
  });
  try {
    await page.goto('/');
    // The app's ordinary initial panel solve completes first. The coupled
    // run below receives no saved coupled parent and starts from scratch.
    await expect(page.locator('#status')).toContainText('Solved', { timeout: 45000 });
    const browserExpected = await page.evaluate(async id =>
      (await import('/scripts/validation/gui-default-coupled-cases.js')).buildGuiDefaultCoupledCase(id), spec.id);
    write('browser-expected', browserExpected);
    await page.locator('#flow-model').selectOption('streamtube-bl');
    await page.locator('#preset').selectOption(spec.preset);
    expect(await page.locator('#alpha-number').inputValue()).toBe('4');
    expect(await page.locator('#quad-mach').inputValue()).toBe('0.2');
    expect(await page.locator('#quad-reynolds').inputValue()).toBe('1000000');
    expect(await page.locator('#quad-ncrit').inputValue()).toBe('9');
    expect(await page.locator('#quad-transition').inputValue()).toBe('automatic');
    expect(await page.locator('#euler-ismom').inputValue()).toBe('4');
    await expect(page.locator('#grid-elliptic')).toBeChecked();
    await page.locator('#solve-button').click();
    await expect.poll(() => page.evaluate(() => window.guiDefaultCapture.request !== null), { timeout: 15000 }).toBe(true);
    const request = await page.evaluate(() => window.guiDefaultCapture.request);
    write('request', request);
    expect(request.caseData).toEqual(browserExpected);
    const inputComparison = compareGuiInput(request.caseData, expected); write('input-comparison', inputComparison);
    expect(request.parentResult).toBeUndefined(); expect(request.task).toBeUndefined();
    await expect.poll(() => page.evaluate(() => window.guiDefaultCapture.terminal !== null), { timeout: timeoutMs }).toBe(true);
    const capture = await page.evaluate(() => window.guiDefaultCapture); write('events', capture.events);
    write('terminal', capture.terminal);
    const result = capture.terminal.result;
    const acceptance = result ? assessGuiDefaultCoupledResult(request.caseData, result, result)
      : { passed: false, failures: ['worker-error'], message: capture.terminal.message };
    if (result) write('result', result);
    if (capture.observationError) errors.push(capture.observationError);
    const receipt = { passed: acceptance.passed && errors.length === 0 && changedSources(sourceHashes).length === 0,
      seconds: (Date.now() - started) / 1000, monotonicSeconds: (performance.now() - startedMonotonic) / 1000, acceptance, sourceHashes, sourceChanges: changedSources(sourceHashes), errors,
      inputHash: sha256(path.join(output, 'request.json')), resultHash: result ? sha256(path.join(output, 'result.json')) : null,
      conditions: result?.conditions, families: result?.families, inputComparison, physicalAcceptance: false };
    write('report', receipt);
    expect(capture.terminal.type).toBe('result'); expect(receipt.passed, JSON.stringify(acceptance)).toBe(true);
    await expect(page.locator('#status')).toHaveText('Euler/BL converged · research');
    await expect(page.locator('#error-message')).toBeHidden();
    for (const id of ['cl', 'cd', 'cm']) await expect(page.locator('#' + id)).not.toHaveText('—');
    console.log(JSON.stringify({ case: spec.id, passed: receipt.passed, seconds: receipt.seconds, monotonicSeconds: receipt.monotonicSeconds, families: receipt.families, output }));
  } catch (error) {
    const capture = await page.evaluate(() => window.guiDefaultCapture).catch(() => null);
    write('failure', { message: error.message, capture, errors, sourceChanges: changedSources(sourceHashes), seconds: (Date.now() - started) / 1000, monotonicSeconds: (performance.now() - startedMonotonic) / 1000 });
    throw error;
  }
});
