import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';

test('public GUI solves the retained unsmoothed automatic-transition foil and displays actual surface status', async ({ page }) => {
  test.setTimeout(180000);
  const baselinePath = 'docs/current-quad-automatic-mrchue-public.json', baseline = JSON.parse(fs.readFileSync(baselinePath));
  assert.deepEqual(changedSources(baseline.sourceHashes), []);
  assert.equal(baseline.result.converged, true);
  const sourceHashes = numericalSourceHashes(['tests/browser/quad-automatic-default.spec.js', 'index.html',
    'src/ui/app.js', 'src/ui/plots.js', 'src/ui/quad-coupled-result.js', 'src/ui/quad-coupled-coefficients.js',
    'src/ui/quad-mesh-progress.js', 'src/worker/solver.js']);
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args); this.addEventListener('message', ({ data }) => {
          if (data.type === 'result') window.lastResult = data.result;
          if (data.type === 'mesh') window.meshMessages = (window.meshMessages ?? 0) + 1;
        });
      }
    };
  });
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved', { timeout: 45000 });
  await page.locator('#flow-model').selectOption('streamtube-bl');
  // This archived numerical baseline predates the new SLOR-on GUI default.
  await page.locator('#grid-elliptic').uncheck();
  await page.locator('#quad-transition').selectOption('automatic');
  for (const id of ['quad-trip-upper', 'quad-trip-lower']) await expect(page.locator('#' + id)).toHaveValue('1');
  await expect(page.locator('#quad-trip-note')).toContainText('allow natural transition');
  // Mode changes retain each mode's global trip values.
  await page.locator('#quad-transition').selectOption('fixed-trip');
  await expect(page.locator('#quad-trip-upper')).toHaveValue('0.05');
  await page.locator('#quad-transition').selectOption('automatic');
  await expect(page.locator('#quad-trip-upper')).toHaveValue('1');
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler/BL converged · research', { timeout: 120000 });
  const result = await page.evaluate(() => {
    const r = window.lastResult;
    return { converged: r.converged, families: r.families, conditions: r.conditions, cl: r.cl, cd: r.cd, cm: r.cm,
      coefficientStatus: r.coefficientStatus, physicalAcceptance: r.physicalAcceptance,
      surfaceStatus: r.boundaryLayer.surfaces.map(s => ({ element: s.element, side: s.side, forced: s.forced,
        kind: s.transitionKind, fraction: s.transitionFraction, point: s.transitionPoint,
        terminalRegime: s.stations.at(-1).regime, terminalAmplification: s.stations.at(-1).amplification })),
      surfaces: r.boundaryLayer.surfaces.length, wakes: r.boundaryLayer.wakes.length,
      restart: r.restart, meshMessages: window.meshMessages, initialization: { factor: r.initialization.boundaryLayer.thicknessFactor,
        eulerIterations: r.initialization.euler.iterations, coupledIterations: r.history.length - 1 } };
  });
  expect(result.conditions.transitionMode).toBe('automatic'); expect(result.physicalAcceptance).toBe(false);
  expect(result.surfaces).toBe(4); expect(result.wakes).toBe(2); expect(result.meshMessages).toBeGreaterThan(20);
  expect(Math.max(...Object.values(result.families))).toBeLessThan(1e-10); expect(result.initialization.factor).toBe(1);
  expect(result.coefficientStatus).toBe('research-unvalidated'); expect(result.cd).toBeGreaterThan(0);
  for (const key of ['cl', 'cd', 'cm']) {
    expect(Number.isFinite(result[key])).toBe(true);
    await expect(page.locator('#' + key)).toHaveText(result[key].toFixed(key === 'cd' ? 6 : 4));
  }
  expect(result.surfaceStatus.filter(s => s.kind === 'natural').length).toBe(2);
  expect(result.surfaceStatus.filter(s => s.kind === 'laminar-to-te').length).toBe(2);
  for (const s of result.surfaceStatus) {
    expect(s.forced).toBe(false);
    if (s.kind === 'natural') expect(s.fraction).toBeLessThan(.2);
    else { expect(Math.abs(s.fraction - 1)).toBeLessThan(1e-12); expect(s.terminalRegime).toBe('laminar'); expect(s.terminalAmplification).toBeLessThan(9); }
  }
  await expect(page.locator('#diagnostics')).toContainText('Natural'); await expect(page.locator('#diagnostics')).toContainText('Laminar to TE');
  await expect(page.locator('#coefficient-warning')).toContainText('Unvalidated');
  let maximumStateDifference = 0, maximumNormalizedStateDifference = 0, maximumNodeDifference = 0, maximumSurfaceNodeDifference = 0,
    largestStateDifference, largestNodeDifference;
  const old = [...baseline.restart.initialEuler.x, ...baseline.restart.initialBL], current = [...result.restart.initialEuler.x, ...result.restart.initialBL];
  expect(current.length).toBe(old.length); current.forEach((v, i) => {
    const difference = Math.abs(v - old[i]);
    if (difference > maximumStateDifference) { maximumStateDifference = difference; largestStateDifference = { column: i, browser: v, node: old[i] }; }
    maximumNormalizedStateDifference = Math.max(maximumNormalizedStateDifference, difference / Math.max(1, Math.abs(v), Math.abs(old[i])));
  });
  result.restart.initialEuler.nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    const q = baseline.restart.initialEuler.nodes[g][i][j], error = Math.hypot(p.x - q.x, p.y - q.y);
    if (error > maximumNodeDifference) { maximumNodeDifference = error; largestNodeDifference = { g, i, j, browser: p, node: q }; }
    baseline.restart.input.bodies.forEach((body, b) => {
      if (i >= body.leadingIndex && i <= body.trailingIndex && (g === b && j === row.length - 1 || g === b + 1 && j === 0))
        maximumSurfaceNodeDifference = Math.max(maximumSurfaceNodeDifference, error);
    });
  })));
  const nodes = baseline.restart.initialEuler.nodes.flat(2);
  const domainScale = Math.max(baseline.input.referenceChord,
    Math.max(...nodes.map(p => p.x)) - Math.min(...nodes.map(p => p.x)),
    Math.max(...nodes.map(p => p.y)) - Math.min(...nodes.map(p => p.y)));
  // Internal thicknesses are scaled by sqrt(Re) and may exceed unity. Use
  // the same normalized comparison for every primitive; retain raw errors.
  // The earlier absolute-only failure is saved in the before-browser log.
  const downloading = page.waitForEvent('download'); await page.locator('#export-button').click();
  const download = await downloading, exported = JSON.parse(fs.readFileSync(await download.path()));
  expect(exported.input.transitionMode).toBe('automatic');
  expect(exported.result.boundaryLayer.surfaces.map(s => s.transitionKind)).toEqual(result.surfaceStatus.map(s => s.kind));
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  await page.screenshot({ path: 'docs/current-quad-automatic-mrchue-gui.png', fullPage: true });
  expect(errors).toEqual([]); expect(changedSources(sourceHashes)).toEqual([]);
  const restart = result.restart;
  delete result.restart;
  const comparison = { maximumStateDifference, maximumNormalizedStateDifference, largestStateDifference, maximumNodeDifference,
    domainScale, maximumNormalizedNodeDifference: maximumNodeDifference / domainScale, maximumSurfaceNodeDifference, largestNodeDifference,
    originalAbsoluteNodeLimitPassed: maximumNodeDifference < 1e-9 };
  fs.writeFileSync('docs/current-quad-automatic-mrchue-gui.json', JSON.stringify({ date: new Date().toISOString(), sourceHashes,
    source: { path: baselinePath, sha256: sha256(baselinePath) }, physicalAcceptance: false,
    scope: 'Fresh public GUI/worker automatic default solve, actual per-surface transition/laminar-to-TE display, coefficients, export, mode controls and phone layout. Complete state/geometry agree with the audited Node cold case; independent physical accuracy remains open.',
    ...comparison, result, restart }, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v, 2) + '\n');
  console.log(JSON.stringify({ result, comparison }));
  expect(maximumNormalizedStateDifference).toBeLessThan(1e-9);
  // Free-grid coordinates scale with the complete farfield domain. Keep a
  // separate chord-scale comparison on both displaced walls of every foil.
  expect(maximumNodeDifference / domainScale).toBeLessThan(1e-9);
  expect(maximumSurfaceNodeDifference / baseline.input.referenceChord).toBeLessThan(1e-9);
});
