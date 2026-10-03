import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('default GUI solves four BLs and two wakes with resolved automatic transition roots', async ({ page }) => {
  test.setTimeout(180000);
  const sourceHashes = numericalSourceHashes(['tests/browser/transition-root-precision.spec.js', 'index.html', 'src/ui/app.js',
    'src/ui/quad-coupled-result.js', 'src/worker/solver.js']);
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args); this.addEventListener('message', ({ data }) => {
          if (data.type === 'result') window.precisionResult = data.result;
          if (data.type === 'mesh') window.precisionMeshes = (window.precisionMeshes ?? 0) + 1;
        });
      }
    };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved', { timeout: 45000 });
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await page.locator('#quad-transition').selectOption('automatic');
  await expect(page.locator('#quad-trip-upper')).toHaveValue('1'); await expect(page.locator('#quad-trip-lower')).toHaveValue('1');
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler/BL converged · research', { timeout: 120000 });
  const captured = await page.evaluate(() => {
    const r = window.precisionResult;
    return JSON.parse(JSON.stringify({ restart: r.restart, result: { converged: r.converged, families: r.families,
      conditions: r.conditions, cl: r.cl, cd: r.cd, cm: r.cm, coefficientStatus: r.coefficientStatus,
      physicalAcceptance: r.physicalAcceptance, surfaces: r.boundaryLayer.surfaces.length, wakes: r.boundaryLayer.wakes.length,
      transitions: r.boundaryLayer.surfaces.map(s => ({ element: s.element, side: s.side, kind: s.transitionKind,
        fraction: s.transitionFraction, point: s.transitionPoint })), history: r.history,
      meshMessages: window.precisionMeshes, initialization: r.initialization } },
    (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
  });
  const r = captured.result;
  expect(r.conditions.transitionMode).toBe('automatic'); expect(r.physicalAcceptance).toBe(false);
  expect(r.surfaces).toBe(4); expect(r.wakes).toBe(2); expect(r.meshMessages).toBeGreaterThan(20);
  expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-10);
  expect(r.transitions.filter(t => t.kind === 'natural')).toHaveLength(2);
  expect(r.transitions.filter(t => t.kind === 'laminar-to-te')).toHaveLength(2);
  expect(r.coefficientStatus).toBe('research-unvalidated');
  for (const key of ['cl', 'cd', 'cm']) {
    expect(Number.isFinite(r[key])).toBe(true);
    await expect(page.locator('#' + key)).toHaveText(r[key].toFixed(key === 'cd' ? 6 : 4));
  }
  const downloading = page.waitForEvent('download'); await page.locator('#export-button').click();
  const download = await downloading, exported = JSON.parse(fs.readFileSync(await download.path()));
  expect(exported.input.transitionMode).toBe('automatic');
  expect(exported.result.boundaryLayer.surfaces.map(s => s.transitionKind)).toEqual(r.transitions.map(s => s.kind));
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
  await page.screenshot({ path: 'docs/current-transition-root-precision-gui.png', fullPage: true });
  expect(errors).toEqual([]); expect(changedSources(sourceHashes)).toEqual([]);
  fs.writeFileSync('docs/current-transition-root-precision-gui.json', JSON.stringify({ date: new Date().toISOString(), sourceHashes,
    physicalAcceptance: false, inProgress: false, scope: 'Fresh default public GUI solve with resolved automatic roots: four surface BLs, two wakes, live mesh, coefficients, export and phone layout. Complete browser restart retained for independent numerical audit.',
    ...captured }) + '\n');
  console.log(JSON.stringify({ families: r.families, iterations: r.history.length - 1, transitions: r.transitions,
    surfaces: r.surfaces, wakes: r.wakes, meshMessages: r.meshMessages }));
});
