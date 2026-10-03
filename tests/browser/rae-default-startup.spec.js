import fs from 'node:fs';
import { test, expect } from '@playwright/test';
import { raeSourceHash } from '../../scripts/validation/rae-status.js';

// Observe the actual worker; never replace its request, solver, or result.
test('RAE 64x9 converges from the normal inviscid form with Standard startup', async ({ page, browser }) => {
  test.setTimeout(900_000);
  const sourceHash = raeSourceHash();
  await page.addInitScript(() => {
    const NativeWorker = window.Worker;
    window.raeRuns = [];
    window.Worker = class extends NativeWorker {
      postMessage(request, ...rest) {
        const run = { input: request.caseData, id: request.id, progress: [] };
        window.raeRuns.push(run);
        this.addEventListener('message', ({ data }) => {
          if (data.id !== run.id) return;
          if (data.type === 'iteration') run.progress.push(data.iteration);
          if (data.type === 'error') run.error = data.message;
          if (data.type === 'result') {
            const r = data.result;
            run.result = { status: r.status, mach: r.mach, alpha: r.alpha,
              residual: r.diagnostics?.equationResidual, iterations: r.diagnostics?.iterations,
              reason: r.diagnostics?.reason, settings: r.solverSettings,
              cells: r.mesh?.cells?.length, gridValid: r.mesh?.quality?.valid,
              upwind: r.flow?.solverInput?.upwind, elapsedMs: data.elapsed };
          }
        });
        return super.postMessage(request, ...rest);
      }
    };
  });
  await page.goto('/');
  await page.locator('#stop-button').click();
  await page.locator('#preset').selectOption('rae2822-mses');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#grid-intervals').selectOption('64');
  await page.locator('#grid-tubes').selectOption('9');
  await page.locator('#quad-mach').fill('0.74');
  await page.locator('#alpha-number').fill('2.68');
  await expect(page.locator('#euler-startup')).toHaveValue('standard');
  await page.locator('#solve-button').click();
  await page.waitForFunction(() => {
    const run = window.raeRuns.at(-1);
    return run?.input?.gridIntervals === 64 && (run.result || run.error);
  }, null, { timeout: 850_000 });
  const run = await page.evaluate(() => window.raeRuns.at(-1));
  expect(raeSourceHash()).toBe(sourceHash);
  const receipt = { ...run, sourceHash, route: 'browser-form', browser: browser.version(), recordedAt: new Date().toISOString() };
  if (process.env.RAE_APP_RECEIPT) fs.writeFileSync(process.env.RAE_APP_RECEIPT, JSON.stringify(receipt, null, 2) + '\n');
  console.log(JSON.stringify(run.result ?? run.error));
  expect(run.error).toBeUndefined();
  expect(run.input).toMatchObject({ eulerStartup: 'standard', gridIntervals: 64, gridTubes: 9,
    mach: .74, alpha: 2.68, quadBoundaryLayers: false, eulerIsmom: 4 });
  expect(run.result.status).toBe('research-converged');
  expect(run.result.settings).toMatchObject({ requestedEulerStartup: 'standard', eulerStartup: 'harmonic-shock', maxIterations: 80 });
  expect(run.result.residual).toBeLessThanOrEqual(1e-10);
  expect(run.result).toMatchObject({ mach: .74, alpha: 2.68, gridValid: true, cells: 9120 });
  expect(run.result.upwind).toMatchObject({ mucon: 1, mcrit: .99 });
  await expect(page.locator('#status')).toHaveText('Euler converged · research');
});
