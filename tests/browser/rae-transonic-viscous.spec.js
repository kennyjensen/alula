import { test, expect } from '@playwright/test';
import { buildReliabilityCase } from '../../scripts/validation/solver-reliability-cases.js';
import { raeSourceHash } from '../../scripts/validation/rae-status.js';

for (const [intervals, tubes, ends, reynolds, ncrit, cells, mach = .74] of [[64, 11, 32, 2710000, 4, 7504], [16, 7, 16, 1000000, 9, 1340],
  [16, 7, 16, 2710000, 4, 1340], [16, 9, 16, 1000000, 9, 1920], [8, 7, 16, 1000000, 9, 1340],
  [32, 7, 16, 1000000, 9, 2020], [32, 9, 16, 1000000, 9, 3216], [16, 7, 16, 1000000, 9, 1600, .76]])
test(`RAE 2822 ${intervals}x${tubes} Re${reynolds} N${ncrit} M${mach} converges viscous at the reported transonic conditions`, async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  const sourceHash = raeSourceHash();
  await page.addInitScript(() => {
    const NativeWorker = window.Worker;
    window.raeRuns = [];
    window.Worker = class extends NativeWorker {
      postMessage(request, ...rest) {
        const run = { input: request.caseData, id: request.id, stages: [], history: [], started: performance.now() };
        window.raeRuns.push(run);
        this.addEventListener('message', ({ data }) => {
          if (data.id !== run.id) return;
          if (data.type === 'error') run.error = data;
          if (data.type === 'flow-stage') run.stages.push({ ...data, elapsed: performance.now() - run.started });
          if (data.type === 'iteration') {
            const h = data.iteration;
            run.history.push({ ...h, elapsed: performance.now() - run.started });
            console.log('RAE iteration', JSON.stringify({ stage: h.stage, iteration: h.iteration, residual: h.residual }));
          }
          if (data.type === 'result') {
            const r = data.result;
            run.result = { status: r.status, mach: r.mach, alpha: r.alpha, residual: r.diagnostics?.equationResidual,
              elapsed: data.elapsed, cl: r.cl, cd: r.cd, cm: r.cm,
              iterations: r.diagnostics?.iterations,
              reason: r.diagnostics?.reason ?? r.initialization?.attempts?.at(-1)?.reason,
              gridValid: r.mesh?.quality?.valid, cellCount: r.mesh?.cells?.length,
              initialization: r.initialization, diagnostics: r.diagnostics, solverSettings: r.solverSettings,
              machRecovery: r.machRecovery, resolutionRecovery: r.resolutionRecovery, quality: r.mesh?.quality, seedRecovery: r.seedRecovery, transitionResolutionRecovery: r.transitionResolutionRecovery };
          }
        });
        return super.postMessage(request, ...rest);
      }
    };
  });
  page.on('console', message => { if (message.text().startsWith('RAE iteration')) console.log(message.text()); });
  await page.goto('/'); await page.locator('#stop-button').click();
  await page.selectOption('#preset', 'rae2822-mses');
  await page.selectOption('#flow-model', 'streamtube-bl');
  await page.fill('#quad-mach', String(mach)); await page.fill('#alpha-number', '2.68');
  await page.fill('#quad-reynolds', String(reynolds)); await page.fill('#quad-ncrit', String(ncrit));
  await page.selectOption('#quad-transition', 'automatic');
  // These regression cases use natural transition. The RAE preset now sets
  // benchmark trip limits, so explicitly restore this suite's intended inputs.
  await page.fill('#quad-trip-upper', '1');
  await page.fill('#quad-trip-lower', '1');
  await page.selectOption('#grid-intervals', String(intervals)); await page.selectOption('#grid-tubes', String(tubes));
  await page.selectOption('#grid-inlet', String(ends)); await page.selectOption('#grid-outlet', String(ends));
  await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await page.click('#solve-button');
  await page.waitForFunction(tubes => {
    const r = window.raeRuns.at(-1);
    return r?.input?.gridTubes === tubes && (r.result || r.error);
  }, tubes, { timeout: 550_000 });
  const run = await page.evaluate(() => window.raeRuns.at(-1));
  run.sourceHash = sourceHash;
  run.route = 'browser-form';
  run.recordedAt = new Date().toISOString();
  expect(raeSourceHash()).toBe(sourceHash);
  const receipt = testInfo.outputPath(`rae${intervals}x${tubes}-run.json`);
  await (await import('node:fs/promises')).writeFile(receipt, JSON.stringify(run));
  await testInfo.attach(`rae${intervals}x${tubes}-run`, { path: receipt, contentType: 'application/json' });
  console.log(JSON.stringify(run.result ? { ...run.result, initialization: undefined } : run.error));
  expect(run.input).toMatchObject({ mach, alpha: 2.68, reynolds, ncrit,
    gridIntervals: intervals, gridTubes: tubes, gridInletIntervals: ends, gridOutletIntervals: ends, quadBoundaryLayers: true });
  expect(run.input.elements[0].points).toHaveLength(129);
  const { caseData } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-bl',
    changes: { mach, alpha: 2.68, reynolds, ncrit, gridIntervals: intervals, gridTubes: tubes, gridInletIntervals: ends, gridOutletIntervals: ends } });
  const { elements: expectedElements, ...expectedControls } = caseData;
  const { elements: actualElements, ...actualControls } = run.input;
  expect(actualControls).toEqual(expectedControls);
  expect(actualElements).toHaveLength(expectedElements.length);
  expect(actualElements).toEqual(expectedElements);
  expect(run.error).toBeUndefined();
  expect(run.result).toMatchObject({ mach, alpha: 2.68, gridValid: true });
  expect(['research-coupled-equations-converged', 'research-coupled-changes-converged']).toContain(run.result.status);
  if (run.result.status === 'research-coupled-equations-converged') {
    expect(run.result.residual).toBeLessThanOrEqual(1e-10);
  } else {
    const convergence = run.result.diagnostics.convergence;
    expect(convergence).toMatchObject({ method: 'mses-changes', converged: true, stable: true,
      rmsTolerance: 1e-4, maximumTolerance: 1e-3 });
    // Small damped steps alone must never certify a stalled solve. Require
    // the full Newton corrections and equation families to pass as well.
    for (const key of ['accepted', 'newton', 'residual']) {
      const families = Object.values(convergence[key]);
      expect(families.length).toBeGreaterThan(0);
      for (const family of families) {
        expect(Number.isInteger(family.count)).toBe(true);
        expect(family.count).toBeGreaterThan(0);
        expect(family.rms).toBeGreaterThanOrEqual(0);
        expect(family.rms).toBeLessThanOrEqual(1e-4);
        expect(family.maximum).toBeGreaterThanOrEqual(0);
        expect(family.maximum).toBeLessThanOrEqual(1e-3);
      }
    }
    expect(run.result.residual).toBeLessThanOrEqual(1e-3);
  }
  expect(run.result.diagnostics.upwind).toMatchObject({ mucon: 1, mcrit: .99 });
  expect(run.result.diagnostics.hybrid.ismom).toBe(4);
  expect(run.result.quality.minCornerSine).toBeGreaterThan(.05);
  expect(run.result.diagnostics.maxMach).toBeGreaterThan(1);
  expect(run.result.cellCount).toBe(cells);
  const transitionRestart = intervals === 16 && (tubes === 9 || mach === .76);
  if (transitionRestart) {
    const recovery = run.result.transitionResolutionRecovery;
    const budget = run.result.solverSettings.maxIterations;
    expect(budget).toBe(40);
    expect(recovery.accepted).toBe(true);
    expect(recovery.plan.normalFactor).toBe(1);
    // Total work includes startup extensions and the discarded seed retry.
    // Bound each phase instead of applying the old single-startup total cap.
    expect(recovery.originalIterations).toBeLessThanOrEqual(3 * budget);
    expect(recovery.precedingIterations).toBeLessThanOrEqual(recovery.originalIterations + budget);
    expect(recovery.iterations).toBeLessThanOrEqual(budget);
    expect(run.result.iterations).toBe(recovery.precedingIterations + recovery.iterations);
  } else expect(run.result.iterations).toBeLessThanOrEqual(160);
  const resolutionRestart = intervals === 8 || intervals === 32 && tubes === 7;
  if (resolutionRestart) {
    expect(run.result.resolutionRecovery).toMatchObject({ accepted: true, requestedGridIntervals: intervals,
      finalNominalGridIntervals: Math.max(16, intervals), operatingConditionsChanged: false, equationsChanged: false });
  }
  if (intervals >= 32) {
    expect(run.result.initialization.coarseInitialization).toMatchObject({ accepted: true,
      requestedGrid: { intervals, tubes }, equationsChanged: false,
      mapping: { nominalGridIntervals: intervals, initialGuessOnly: true } });
  }
  if (intervals < 32 || resolutionRestart) {
    const seed = run.result.initialization.boundaryLayer.panelEdgeGuess;
    expect(seed?.accepted === true || run.result.seedRecovery?.accepted === true).toBe(true);
  }

});
