import fs from 'node:fs';
import { test, expect } from '@playwright/test';
import { raeSourceHash } from '../../scripts/validation/rae-status.js';

// Exercise automatic defaults or explicit mesh overrides through the real
// form. Only exact default receipts can certify docs/rae-status.md.
const tubeCounts = (process.env.RAE_FINE_TUBES ?? '7,9,11,24').split(',').map(Number);
const spacing = process.env.RAE_FINE_SPACING ?? 'source';
const requestedEnds = process.env.RAE_FINE_ENDS ?? '64';
const configurations = tubeCounts.map(tubes => ({ tubes, ends: requestedEnds }));
if (spacing === 'automatic' && requestedEnds === 'auto' && tubeCounts.includes(24))
  configurations.push({ tubes: 24, ends: '32' });
for (const { tubes, ends } of configurations) {
  // The denser 24-tube configuration also uses the public stagnation-cell shape
  // control. This is explicit test input, not a hidden solver override.
  const aspectRatio = process.env.RAE_FINE_ASPECT === undefined ? (tubes === 24 && spacing !== 'automatic' ? 2.5 : undefined) : Number(process.env.RAE_FINE_ASPECT);
  test(`RAE 128x${tubes} converges with ${spacing} spacing${aspectRatio === undefined ? '' : ` and stagnation aspect ${aspectRatio}`} / ${ends} end intervals through the app`, async ({ page, browser }) => {
    test.setTimeout(900_000);
    page.setDefaultTimeout(30_000);
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
                startupAttempts: r.diagnostics?.startupAttempts, gridSelection: r.mesh?.initialization?.gridSelection,
                upwind: r.flow?.solverInput?.upwind, elapsedMs: data.elapsed };
            }
          });
          return super.postMessage(request, ...rest);
        }
      };
    });
    await page.goto('/');
    if (await page.locator('#stop-button').isEnabled()) await page.locator('#stop-button').click();
    await page.locator('#preset').selectOption('rae2822-mses');
    await page.locator('#flow-model').selectOption('streamtube-grid');
    await page.locator('#grid-intervals').selectOption('128');
    await page.locator('#grid-tubes').selectOption(String(tubes));
    await page.locator('#grid-inlet').selectOption(ends);
    await page.locator('#grid-outlet').selectOption(ends);
    if (spacing === 'automatic') {
      await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
      await expect(page.locator('#grid-match-aspect')).not.toBeChecked();
    } else {
      await page.locator('#more-mesh-settings > summary').click();
      await page.locator('#grid-surface-spacing').selectOption(spacing);
    }
    if (aspectRatio !== undefined) {
      await page.locator('#grid-match-aspect').check();
      await page.locator('#grid-aspect-ratio').fill(String(aspectRatio));
    }
    if (spacing === 'curvature') {
      await page.locator('#grid-le-ratio').fill('0.8');
      await page.locator('#grid-te-ratio').fill('0.4');
    }
    await page.locator('#quad-mach').fill('0.74');
    await page.locator('#alpha-number').fill('2.68');
    await expect(page.locator('#euler-startup')).toHaveValue('standard');
    console.log('Fine-grid controls set; starting real worker.');
    await page.locator('#solve-button').click();
    await expect.poll(async () => {
      const progress = await page.evaluate(() => {
        const run = window.raeRuns.at(-1);
        return { grid: run?.input?.gridIntervals, last: run?.progress?.at(-1), done: Boolean(run?.result || run?.error) };
      });
      console.log(JSON.stringify({ grid: progress.grid, iteration: progress.last?.iteration, residual: progress.last?.residual, done: progress.done }));
      return progress.grid === 128 && progress.done;
    }, { timeout: 850_000, intervals: [15_000] }).toBe(true);
    const run = await page.evaluate(() => window.raeRuns.at(-1));
    expect(raeSourceHash()).toBe(sourceHash);
    const receipt = { ...run, sourceHash, route: 'browser-form', browser: browser.version(), recordedAt: new Date().toISOString() };
    if (process.env.RAE_FINE_RECEIPT) {
      const path = configurations.length === 1 ? process.env.RAE_FINE_RECEIPT
        : `${process.env.RAE_FINE_RECEIPT}-${tubes}${ends === requestedEnds ? '' : `-ends${ends}`}.json`;
      fs.writeFileSync(path, JSON.stringify(receipt, null, 2) + '\n');
    }
    console.log(JSON.stringify(run.result ?? run.error));
    expect(run.error).toBeUndefined();
    expect(run.input).toMatchObject({ eulerStartup: 'standard', gridIntervals: 128, gridTubes: tubes,
      ...(ends === 'auto' ? {} : { gridInletIntervals: Number(ends), gridOutletIntervals: Number(ends) }), gridSurfaceSpacing: spacing,
      ...(aspectRatio === undefined ? {} : { gridStagnationAspectRatio: aspectRatio }),
      ...(spacing === 'curvature' ? { gridCurvatureSpacing: { exponent: .5, leadingSpacingRatio: .8, trailingSpacingRatio: .4 } } : {}),
      mach: .74, alpha: 2.68, quadBoundaryLayers: false, eulerIsmom: 4 });
    if (ends === 'auto') {
      expect(run.input.gridInletIntervals).toBeUndefined();
      expect(run.input.gridOutletIntervals).toBeUndefined();
    }
    if (spacing === 'automatic') {
      expect(run.result.gridSelection).toMatchObject({ balanced: true, surfaceSpacing: 'source', normalSpacing: 'automatic' });
      if (tubes === 24 && ends === '32') expect(run.result.startupAttempts).toHaveLength(2);
    }
    expect(run.result.status).toBe('research-converged');
    expect(run.input.gridStagnationAspectRatio).toBe(aspectRatio);
    expect(run.result.settings).toMatchObject({ requestedEulerStartup: 'standard', eulerStartup: 'harmonic-shock', maxIterations: 80 });
    expect(run.result.residual).toBeLessThanOrEqual(1e-10);
    expect(run.result).toMatchObject({ mach: .74, alpha: 2.68, gridValid: true });
    expect(run.result.upwind).toMatchObject({ mucon: 1, mcrit: .99 });
    await expect(page.locator('#status')).toHaveText('Euler converged · research');
  });
}
