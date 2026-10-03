import { test, expect } from '@playwright/test';
import { buildReliabilityCase } from '../../scripts/validation/solver-reliability-cases.js';
import { naca4Standard, transform } from '../../src/geometry/airfoil.js';
import { prepareAirfoilElement } from '../../src/geometry/airfoil-element.js';
import { raeSourceHash } from '../../scripts/validation/rae-status.js';

test('three-element custom assembly converges viscous after grid recovery', async ({ page }, testInfo) => {
  test.setTimeout(900_000);
  const sourceHash = raeSourceHash();
  await page.addInitScript(() => {
    const NativeWorker = window.Worker;
    window.customRuns = [];
    window.Worker = class extends NativeWorker {
      postMessage(request, ...rest) {
        const run = { input: request.caseData, id: request.id, stages: [], history: [], started: performance.now() };
        window.customRuns.push(run);
        this.addEventListener('message', ({ data }) => {
          if (data.id !== run.id) return;
          if (data.type === 'error') run.error = data;
          if (data.type === 'flow-stage') run.stages.push({ ...data, elapsed: performance.now() - run.started });
          if (data.type === 'iteration') {
            const h = data.iteration;
            run.history.push({ ...h, elapsed: performance.now() - run.started });
            console.log('CUSTOM32 iteration', JSON.stringify({ stage: h.stage, iteration: h.iteration, residual: h.residual }));
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
  page.on('console', message => { if (message.text().startsWith('CUSTOM32 iteration')) console.log(message.text()); });
  await page.goto('/'); await page.locator('#stop-button').click();
  await page.selectOption('#preset', 'three');
  await page.locator('.element-card').first().locator('summary').click();
  await page.locator('[data-element="0"][data-key="deflection"]').fill('-12');
  await page.locator('[data-element="0"][data-key="deflection"]').press('Tab');
  await page.selectOption('#flow-model', 'streamtube-bl');
  await page.fill('#quad-mach', '0.2'); await page.fill('#alpha-number', '2');
  await page.fill('#quad-reynolds', '1000000'); await page.fill('#quad-ncrit', '4');
  await page.selectOption('#quad-transition', 'automatic');
  await page.selectOption('#grid-intervals', '32'); await page.selectOption('#grid-tubes', '9');
  await page.selectOption('#grid-inlet', '16'); await page.selectOption('#grid-outlet', '16');
  await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await page.click('#solve-button');
  await page.waitForFunction(() => {
    const r = window.customRuns.at(-1);
    return r?.input?.gridTubes === 9 && (r.result || r.error);
  }, null, { timeout: 850_000 });
  const run = await page.evaluate(() => window.customRuns.at(-1));
  run.sourceHash = sourceHash;
  run.route = 'browser-form';
  run.recordedAt = new Date().toISOString();
  expect(raeSourceHash()).toBe(sourceHash);
  const receipt = testInfo.outputPath('custom32-run.json');
  await (await import('node:fs/promises')).writeFile(receipt, JSON.stringify(run));
  await testInfo.attach('custom32-run', { path: receipt, contentType: 'application/json' });
  console.log(JSON.stringify(run.result ? { ...run.result, initialization: undefined } : run.error));
  expect(run.input).toMatchObject({ mach: .2, alpha: 2, reynolds: 1000000, ncrit: 4,
    gridIntervals: 32, gridTubes: 9, gridInletIntervals: 16, gridOutletIntervals: 16, quadBoundaryLayers: true });
  expect(run.input.elements.map(e => e.trailingEdge.lowerIndex)).toEqual([160, 160, 160]);
  const { caseData } = buildReliabilityCase({ preset: 'three', mode: 'streamtube-bl',
    changes: { mach: .2, alpha: 2, ncrit: 4, gridIntervals: 32, gridTubes: 9, gridInletIntervals: 16, gridOutletIntervals: 16 } });
  caseData.elements[0] = prepareAirfoilElement({ name: 'Slat', points: transform(naca4Standard('0012', 160),
    { chord: .2, x: -.22, y: .05, angle: 12 }) });
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
  expect(run.result).toMatchObject({ status: 'research-coupled-equations-converged', mach: .2, alpha: 2, gridValid: true });
  expect(run.result.residual).toBeLessThanOrEqual(1e-10);
  expect(run.result.cellCount).toBe(12474);
  expect(run.result.iterations).toBeLessThanOrEqual(40);
  expect(Math.max(...run.history.filter(h => h.stage === 'euler').map(h => h.iteration))).toBeLessThanOrEqual(40);
  expect(run.history.some(h => h.startupStrategy === 'passage-grid-recovery')).toBe(true);
  const seed = run.result.initialization.boundaryLayer.panelEdgeGuess;
  expect(seed.accepted).toBe(true);
  expect(seed.profile.accepted).toBe(false);
  expect(seed.profile.reason).toContain('lost resolved transition');
});
