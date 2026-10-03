import { test, expect } from '@playwright/test';
import fs from 'node:fs/promises';
import { buildReliabilityCase } from '../../scripts/validation/solver-reliability-cases.js';
import { assessSolverResult } from '../../scripts/validation/solver-reliability-acceptance.js';
import { raeSourceHash } from '../../scripts/validation/rae-status.js';

test('flap default viscous form converges with automatic inlet and outlet', async ({ page }, testInfo) => {
  const preset = 'flap';
  test.setTimeout(600_000);
  const sourceHash = raeSourceHash();
  await page.addInitScript(() => {
    const NativeWorker = window.Worker;
    window.presetRuns = [];
    window.Worker = class extends NativeWorker {
      postMessage(request, ...rest) {
        const run = { id: request.id, input: request.caseData, history: [] };
        window.presetRuns.push(run);
        this.addEventListener('message', ({ data }) => {
          if (data.id !== run.id) return;
          if (data.type === 'error') run.error = data;
          if (data.type === 'iteration') {
            const h = data.iteration;
            run.history.push({ stage: h.stage, iteration: h.iteration, residual: h.residual });
            if (h.iteration % 10 === 0) console.log('preset iteration', JSON.stringify(run.history.at(-1)));
          }
          if (data.type === 'result') {
            const r = data.result;
            run.elapsed = data.elapsed;
            run.result = { converged: r.converged, status: r.status, mach: r.mach, alpha: r.alpha,
              cl: r.cl, cd: r.cd, cm: r.cm, families: r.families, diagnostics: r.diagnostics,
              convergence: r.convergence, solverSettings: r.solverSettings,
              conditions: r.conditions, boundaryLayer: r.boundaryLayer,
              mesh: { quality: r.mesh?.quality }, cells: r.mesh?.cells?.length,
              materialTrips: r.materialTrips, referenceReynolds: r.referenceReynolds,
              initialization: r.initialization };
          }
        });
        return super.postMessage(request, ...rest);
      }
    };
  });
  page.on('console', m => { if (m.text().startsWith('preset iteration')) console.log(m.text()); });
  await page.goto('/'); await page.locator('#stop-button').click();
  await page.selectOption('#preset', preset);
  await page.selectOption('#flow-model', 'streamtube-bl');
  await page.fill('#quad-mach', '0.2'); await page.fill('#alpha-number', '4');
  await page.fill('#quad-reynolds', '1000000'); await page.fill('#quad-ncrit', '9');
  await page.selectOption('#quad-transition', 'automatic');
  await page.selectOption('#grid-intervals', '16'); await page.selectOption('#grid-tubes', '7');
  await page.selectOption('#grid-inlet', 'auto'); await page.selectOption('#grid-outlet', 'auto');
  await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
  await expect(page.locator('#grid-elliptic')).toBeChecked();
  await page.click('#solve-button');
  await page.waitForFunction(() => {
    const r = window.presetRuns.at(-1);
    return r?.input?.quadBoundaryLayers && (r.result || r.error);
  }, undefined, { timeout: 550_000 });
  const run = await page.evaluate(() => window.presetRuns.at(-1));
  run.sourceHash = sourceHash; run.finalSourceHash = raeSourceHash();
  run.route = 'browser-form'; run.recordedAt = new Date().toISOString();
  run.acceptance = assessSolverResult(run.input, run.result);
  const receipt = testInfo.outputPath(`${preset}-default-run.json`);
  await fs.writeFile(receipt, JSON.stringify(run));
  await testInfo.attach('default-run', { path: receipt, contentType: 'application/json' });
  expect(run.finalSourceHash).toBe(sourceHash);
  const { caseData } = buildReliabilityCase({ preset, mode: 'streamtube-bl' });
  const { elements: actual, ...controls } = run.input;
  const { elements: expected, ...expectedControls } = caseData;
  expect(controls).toEqual(expectedControls);
  expect(actual).toHaveLength(expected.length);
  actual.forEach((e, k) => {
    expect(e.trailingEdge).toEqual(expected[k].trailingEdge);
    for (const key of ['points', 'sourcePoints']) {
      expect(e[key]).toHaveLength(expected[k][key].length);
      e[key].forEach((p, i) => {
        expect(p.x).toBeCloseTo(expected[k][key][i].x, 14);
        expect(p.y).toBeCloseTo(expected[k][key][i].y, 14);
      });
    }
  });
  expect(run.error).toBeUndefined();
  expect(run.acceptance.passed, JSON.stringify(run.acceptance)).toBe(true);
  expect(run.result.referenceReynolds).toBe(1e6);
  expect(run.result.conditions.ncrit).toBe(9);
  expect(run.result.materialTrips).toEqual(actual.map(() => [1, 1]));
  expect(run.result.cells).toBe(3531);
  expect(run.result.diagnostics.convergence.converged).toBe(true);
});
