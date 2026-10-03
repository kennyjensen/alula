import { test, expect } from '@playwright/test';
import { raeSourceHash } from '../../scripts/validation/rae-status.js';

test('NLR 7301 64x24 converges viscous at the reported subsonic conditions', async ({ page }, testInfo) => {
  test.setTimeout(2_400_000);
  const sourceHash = raeSourceHash();
  await page.addInitScript(() => {
    const NativeWorker = window.Worker;
    window.nlrRuns = [];
    window.Worker = class extends NativeWorker {
      postMessage(request, ...rest) {
        const run = { input: request.caseData, id: request.id, stages: [], history: [], started: performance.now() };
        window.nlrRuns.push(run);
        this.addEventListener('message', ({ data }) => {
          if (data.id !== run.id) return;
          if (data.type === 'error') run.error = data;
          if (data.type === 'flow-stage') run.stages.push({ ...data, elapsed: performance.now() - run.started });
          if (data.type === 'iteration') {
            const h = data.iteration;
            run.history.push({ ...h, elapsed: performance.now() - run.started });
            console.log('NLR iteration', JSON.stringify({ stage: h.stage, iteration: h.iteration, residual: h.residual }));
          }
          if (data.type === 'result') {
            const r = data.result;
            run.result = { status: r.status, mach: r.mach, alpha: r.alpha, residual: r.diagnostics?.equationResidual,
              elapsed: data.elapsed, cl: r.cl, cd: r.cd, cm: r.cm,
              iterations: r.diagnostics?.iterations,
              reason: r.diagnostics?.reason ?? r.initialization?.attempts?.at(-1)?.reason,
              gridValid: r.mesh?.quality?.valid, cellCount: r.mesh?.cells?.length,
              initialization: r.initialization };
          }
        });
        return super.postMessage(request, ...rest);
      }
    };
  });
  page.on('console', message => { if (message.text().startsWith('NLR iteration')) console.log(message.text()); });
  await page.goto('/'); await page.locator('#stop-button').click();
  await page.selectOption('#preset', 'nlr7301');
  await page.selectOption('#flow-model', 'streamtube-bl');
  await page.fill('#quad-mach', '0.185'); await page.fill('#alpha-number', '6');
  await page.fill('#quad-reynolds', '2510000'); await page.fill('#quad-ncrit', '9');
  await page.selectOption('#quad-transition', 'automatic');
  await page.selectOption('#grid-intervals', '64'); await page.selectOption('#grid-tubes', '24');
  await page.selectOption('#grid-inlet', '32'); await page.selectOption('#grid-outlet', '32');
  await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await page.click('#solve-button');
  await page.waitForFunction(() => {
    const r = window.nlrRuns.at(-1);
    return r?.input?.gridTubes === 24 && (r.result || r.error);
  }, null, { timeout: 2_350_000 });
  const run = await page.evaluate(() => window.nlrRuns.at(-1));
  run.sourceHash = sourceHash;
  run.route = 'browser-form';
  run.recordedAt = new Date().toISOString();
  expect(raeSourceHash()).toBe(sourceHash);
  const receipt = testInfo.outputPath('nlr-run.json');
  await (await import('node:fs/promises')).writeFile(receipt, JSON.stringify(run));
  await testInfo.attach('nlr-run', { path: receipt, contentType: 'application/json' });
  console.log(JSON.stringify(run.result ?? run.error));
  expect(run.input).toMatchObject({ mach: .185, alpha: 6, reynolds: 2510000, ncrit: 9,
    gridIntervals: 64, gridTubes: 24, gridInletIntervals: 32, gridOutletIntervals: 32, quadBoundaryLayers: true });
  expect(run.input.elements.map(e => e.trailingEdge.lowerIndex)).toEqual([424, 216]);
  expect(run.error).toBeUndefined();
  expect(run.result).toMatchObject({ status: 'research-coupled-equations-converged', mach: .185, alpha: 6, gridValid: true });
  expect(run.result.residual).toBeLessThanOrEqual(1e-10);
  expect(run.result.cellCount).toBe(19824);
  expect(run.result.iterations).toBeLessThanOrEqual(40);
  expect(run.stages.some(s => s.couplingStartup?.method === 'subsonic-pressure-handoff')).toBe(true);
  expect(run.stages.some(s => s.boundaryLayerInitialization?.panelProfile?.accepted === true)).toBe(true);
});
