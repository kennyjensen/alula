import { test, expect } from '@playwright/test';
import { buildReliabilityCase } from '../../scripts/validation/solver-reliability-cases.js';
import { raeSourceHash } from '../../scripts/validation/rae-status.js';

for (const [gridIntervals, gridTubes, alpha, cells, budget] of [[32, 9, 4, 2352, 40], [16, 7, 8, 1340, 120]])
test(`NACA 0012 ${gridIntervals}x${gridTubes} alpha ${alpha} converges viscous at Mach 0.2`, async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  const sourceHash = raeSourceHash();
  await page.addInitScript(() => {
    const NativeWorker = window.Worker;
    window.nacaRuns = [];
    window.Worker = class extends NativeWorker {
      postMessage(request, ...rest) {
        const run = { input: request.caseData, id: request.id, stages: [], history: [], started: performance.now() };
        window.nacaRuns.push(run);
        this.addEventListener('message', ({ data }) => {
          if (data.id !== run.id) return;
          if (data.type === 'error') run.error = data;
          if (data.type === 'flow-stage') run.stages.push({ ...data, elapsed: performance.now() - run.started });
          if (data.type === 'iteration') {
            const h = data.iteration;
            run.history.push({ ...h, elapsed: performance.now() - run.started });
            console.log('NACA32 iteration', JSON.stringify({ stage: h.stage, iteration: h.iteration, residual: h.residual }));
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
  page.on('console', message => { if (message.text().startsWith('NACA32 iteration')) console.log(message.text()); });
  await page.goto('/'); await page.locator('#stop-button').click();
  await page.selectOption('#preset', 'single');
  await page.selectOption('#flow-model', 'streamtube-bl');
  await page.fill('#quad-mach', '0.2'); await page.fill('#alpha-number', String(alpha));
  await page.fill('#quad-reynolds', '1000000'); await page.fill('#quad-ncrit', '9');
  await page.selectOption('#quad-transition', 'automatic');
  await page.selectOption('#grid-intervals', String(gridIntervals)); await page.selectOption('#grid-tubes', String(gridTubes));
  await page.selectOption('#grid-inlet', '16'); await page.selectOption('#grid-outlet', '16');
  await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await page.click('#solve-button');
  await page.waitForFunction(tubes => {
    const r = window.nacaRuns.at(-1);
    return r?.input?.gridTubes === tubes && (r.result || r.error);
  }, gridTubes, { timeout: 550_000 });
  const run = await page.evaluate(() => window.nacaRuns.at(-1));
  run.sourceHash = sourceHash;
  run.route = 'browser-form';
  run.recordedAt = new Date().toISOString();
  expect(raeSourceHash()).toBe(sourceHash);
  const receipt = testInfo.outputPath(`naca${gridIntervals}-alpha${alpha}-run.json`);
  await (await import('node:fs/promises')).writeFile(receipt, JSON.stringify(run));
  await testInfo.attach('naca32-run', { path: receipt, contentType: 'application/json' });
  console.log(JSON.stringify(run.result ? { ...run.result, initialization: undefined } : run.error));
  expect(run.input).toMatchObject({ mach: .2, alpha, reynolds: 1000000, ncrit: 9,
    gridIntervals, gridTubes, gridInletIntervals: 16, gridOutletIntervals: 16, quadBoundaryLayers: true });
  expect(run.input.elements.map(e => e.trailingEdge.lowerIndex)).toEqual([160]);
  const { caseData } = buildReliabilityCase({ preset: 'single', mode: 'streamtube-bl',
    changes: { mach: .2, alpha, gridIntervals, gridTubes, gridInletIntervals: 16, gridOutletIntervals: 16 } });
  const { elements: expectedElements, ...expectedControls } = caseData;
  const { elements: actualElements, ...actualControls } = run.input;
  expect(actualControls).toEqual(expectedControls);
  expect(actualElements).toHaveLength(expectedElements.length);
  for (let e = 0; e < actualElements.length; e++) {
    expect(actualElements[e].trailingEdge).toEqual(expectedElements[e].trailingEdge);
    for (const key of ['points', 'sourcePoints']) {
      expect(actualElements[e][key]).toHaveLength(expectedElements[e][key].length);
      actualElements[e][key].forEach((p, i) => {
        // JS transcendental rounding can differ between Node and Chromium.
        expect(p.x).toBeCloseTo(expectedElements[e][key][i].x, 14);
        expect(p.y).toBeCloseTo(expectedElements[e][key][i].y, 14);
      });
    }
  }
  expect(run.error).toBeUndefined();
  expect(run.result).toMatchObject({ status: 'research-coupled-equations-converged', mach: .2, alpha, gridValid: true });
  expect(run.result.residual).toBeLessThanOrEqual(1e-10);
  expect(run.result.cellCount).toBe(cells);
  expect(run.result.iterations).toBeLessThanOrEqual(budget);
  const seed = run.result.initialization.boundaryLayer.panelEdgeGuess;
  expect(seed.accepted).toBe(true);
  if (alpha === 4) {
    expect(seed.profile.accepted).toBe(false);
    expect(seed.profile.reason).toContain('lost resolved transition');
  }
});
