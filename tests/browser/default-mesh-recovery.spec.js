// SPDX-License-Identifier: GPL-2.0-or-later
// Actual browser/worker regression for the default three-element mesh.
import fs from 'node:fs';
import path from 'node:path';
import { test, expect } from '@playwright/test';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('three-element defaults produce a convex smoothed mesh through the real worker', async ({ page }) => {
  test.setTimeout(210000);
  const output = process.env.MSES_DEFAULT_MESH_OUTPUT ?? 'docs/solver-reliability/three-default-browser';
  expect(fs.existsSync(output), 'Supply a new MSES_DEFAULT_MESH_OUTPUT directory.').toBe(false);
  fs.mkdirSync(output, { recursive: true });
  const serial = x => JSON.stringify(x, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v);
  const write = (name, data) => fs.writeFileSync(path.join(output, `${name}.json`), serial(data) + '\n');
  const sourceHashes = numericalSourceHashes(['src/ui/app.js', 'src/worker/solver.js', 'index.html',
    'tests/browser/default-mesh-recovery.spec.js']);
  const errors = [], start = Date.now();
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    const Original = window.Worker;
    window.meshRecovery = { requests: [], events: [], meshes: [] };
    window.Worker = class extends Original {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => {
          if (!this.meshRequest) return;
          if (data.type === 'mesh') {
            window.meshRecovery.meshes.push(data.mesh);
            window.meshRecovery.events.push({ type: data.type, stage: data.stage, cells: data.mesh.cells.length });
          } else window.meshRecovery.events.push(data);
        });
      }
      postMessage(data, ...args) {
        this.meshRequest = data.meshOnly === true;
        if (this.meshRequest) window.meshRecovery.requests.push(structuredClone(data));
        return super.postMessage(data, ...args);
      }
    };
  });
  try {
    await page.goto('/');
    await expect(page.locator('#status')).toContainText('Solved', { timeout: 45000 });
    await page.locator('#flow-model').selectOption('streamtube-grid');
    await page.locator('#preset').selectOption('three');
    await expect(page.locator('#grid-elliptic')).toBeChecked();
    await expect(page.locator('#grid-intervals')).toHaveValue('16');
    await expect(page.locator('#grid-tubes')).toHaveValue('7');
    await expect(page.locator('#alpha-number')).toHaveValue('4');
    await page.locator('#mesh-button').click();
    await page.waitForFunction(() => window.meshRecovery.events.some(e => e.type === 'mesh-ready' || e.type === 'error'),
      undefined, { timeout: 180000 });
    const capture = await page.evaluate(() => window.meshRecovery);
    write('capture', capture);
    const last = capture.meshes.at(-1), smoothing = last?.initialization?.gridSmoothing;
    const terminal = capture.events.find(e => e.type === 'mesh-ready' || e.type === 'error');
    expect(terminal.type, terminal.message).toBe('mesh-ready');
    expect(capture.requests).toHaveLength(1);
    expect(capture.events.some(e => e.type === 'iteration' || e.type === 'result')).toBe(false);
    expect(smoothing.converged).toBe(true);
    expect(smoothing.harmonicPassages).toContain(3);
    expect(smoothing.regions[3].harmonicFallback).toMatchObject({ converged: true,
      fixedBoundaries: true, retainedAdditionalPoissonSpacingControl: false });
    expect(last.quality.valid).toBe(true);
    // Independent corner checks on the connectivity sent to the browser.
    let minimumCornerSine = Infinity;
    for (const cell of last.cells) {
      expect(cell).toHaveLength(4);
      const p = cell.map(id => last.vertices[id]);
      for (let j = 0; j < 4; j++) {
        const a = p[j], b = p[(j + 1) % 4], c = p[(j + 2) % 4];
        const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
        minimumCornerSine = Math.min(minimumCornerSine, (ux * vy - uy * vx) / (Math.hypot(ux, uy) * Math.hypot(vx, vy)));
      }
    }
    expect(minimumCornerSine).toBeGreaterThan(0);
    await expect(page.locator('#status')).toHaveText('Mesh ready');
    await expect(page.locator('#mesh-status')).toContainText('harmonic SLOR used');
    await page.locator('.geometry-panel').screenshot({ path: path.join(output, 'mesh.png') });
    const changes = changedSources(sourceHashes);
    write('report', { passed: errors.length === 0 && changes.length === 0, seconds: (Date.now() - start) / 1000,
      cells: last.cells.length, minimumCornerSine, smoothing, errors, sourceHashes, sourceChanges: changes,
      scope: 'Cold real-browser mesh-only default. No Euler/BL convergence or physical accuracy claim.' });
    expect(errors).toEqual([]); expect(changes).toEqual([]);
  } catch (error) {
    write('failure', { message: error.message, errors, seconds: (Date.now() - start) / 1000,
      capture: await page.evaluate(() => window.meshRecovery).catch(() => null),
      sourceHashes, sourceChanges: changedSources(sourceHashes) });
    throw error;
  }
});
