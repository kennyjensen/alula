import fs from 'node:fs';
import path from 'node:path';
import { test, expect } from '@playwright/test';
import { raeGrids, raeSourceHash, readRaeReceipts, verifiedRaeReceipt, renderRaeStatus } from '../../scripts/validation/rae-status.js';
import { buildReliabilityCase } from '../../scripts/validation/solver-reliability-cases.js';

// An opt-in inventory, including failures. Exercise the actual form and Worker;
// never inject a solver request or alter the numerical stopping criteria.
const output = process.env.RAE_STATUS_AUDIT;
const existing = readRaeReceipts(), sourceHash = raeSourceHash();
for (const viscous of [false, true]) for (const [n, t] of raeGrids) {
  const id = `rae${n}x${t}-${viscous ? 'viscous' : 'inviscid'}`;
  test(`${id} records the cold browser outcome`, async ({ page, browser }) => {
    test.skip(!output, 'Set RAE_STATUS_AUDIT to an evidence directory. Run with --workers=1.');
    test.skip(existing.some(r => r.input?.gridIntervals === n && r.input?.gridTubes === t
      && r.input?.quadBoundaryLayers === viscous && verifiedRaeReceipt(r, sourceHash)), 'Current browser success already recorded.');
    test.setTimeout(1_860_000);
    page.setDefaultTimeout(30_000);
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
            if (data.type === 'flow-stage') run.stage = data.stage;
            if (data.type === 'error') { run.error = data.message; run.failure = data; }
            if (data.type === 'result') {
              const r = data.result;
              run.result = { status: r.status, mach: r.mach, alpha: r.alpha,
                residual: r.diagnostics?.equationResidual, iterations: r.diagnostics?.iterations,
                reason: r.reason ?? r.diagnostics?.reason, settings: r.solverSettings,
                cells: r.mesh?.cells?.length, gridValid: r.mesh?.quality?.valid,
                upwind: r.solverInput?.upwind ?? r.flow?.solverInput?.upwind,
                diagnostics: r.diagnostics, lastRejectedStep: r.lastRejectedStep,
                elapsedMs: data.elapsed };
            }
          });
          return super.postMessage(request, ...rest);
        }
      };
    });
    await page.goto('/');
    if (await page.locator('#stop-button').isEnabled()) await page.locator('#stop-button').click();
    await page.locator('#preset').selectOption('rae2822-mses');
    await page.locator('#flow-model').selectOption(viscous ? 'streamtube-bl' : 'streamtube-grid');
    await page.locator('#grid-intervals').selectOption(String(n));
    await page.locator('#grid-tubes').selectOption(String(t));
    await page.locator('#quad-mach').fill('0.74');
    await page.locator('#alpha-number').fill('2.68');
    if (viscous) {
      await page.locator('#quad-reynolds').fill('2700000');
      await page.locator('#quad-ncrit').fill('4');
    }
    await expect(page.locator('#grid-surface-spacing')).toHaveValue('automatic');
    await expect(page.locator('#grid-inlet')).toHaveValue('auto');
    await expect(page.locator('#grid-outlet')).toHaveValue('auto');
    await page.locator('#solve-button').click();
    const started = Date.now();
    let timedOut = false;
    try {
      await expect.poll(async () => {
        const p = await page.evaluate(() => {
          const r = window.raeRuns.at(-1);
          return { stage: r?.stage, last: r?.progress.at(-1), done: !!(r?.result || r?.error) };
        });
        console.log(JSON.stringify({ id, seconds: Math.round((Date.now() - started) / 1000),
          stage: p.stage, iteration: p.last?.iteration, residual: p.last?.residual, done: p.done }));
        return p.done;
      }, { timeout: 1_800_000, intervals: [30_000] }).toBe(true);
    } catch (error) {
      if (Date.now() - started < 1_790_000) throw error;
      timedOut = true;
    }
    const run = await page.evaluate(() => window.raeRuns.at(-1));
    const receipt = { ...run, sourceHash, route: 'browser-form', browser: browser.version(),
      recordedAt: new Date().toISOString(), ...(timedOut ? { error: 'Browser audit stopped after 30 minutes; convergence not established.', timedOut: true } : {}) };
    fs.mkdirSync(output, { recursive: true });
    fs.writeFileSync(path.join(output, `${id}.json`), JSON.stringify(receipt, null, 2) + '\n');
    expect(raeSourceHash()).toBe(sourceHash);
    const { caseData } = buildReliabilityCase({ preset: 'rae2822-mses', mode: viscous ? 'streamtube-bl' : 'streamtube-grid',
      changes: { mach: .74, alpha: 2.68, gridIntervals: n, gridTubes: t,
        ...(viscous ? { reynolds: 2.7e6, ncrit: 4 } : { eulerStartup: 'standard' }) } });
    expect(run.input).toEqual(caseData);
    const { progress, ...summary } = receipt;
    summary.lastProgress = progress.at(-1);
    summary.auditOutcome = verifiedRaeReceipt(receipt, sourceHash) ? 'converged' : timedOut ? 'timeout' : 'unconverged';
    const receipts = readRaeReceipts().filter(r => r.input?.gridIntervals !== n || r.input?.gridTubes !== t
      || r.input?.quadBoundaryLayers !== viscous);
    receipts.push(summary);
    fs.writeFileSync('docs/rae-status-results.json', JSON.stringify(receipts, null, 2) + '\n');
    fs.writeFileSync('docs/rae-status.md', renderRaeStatus(receipts, sourceHash));
    console.log(JSON.stringify({ id, outcome: summary.auditOutcome, result: run.result, error: receipt.error }));
    if (timedOut) await page.locator('#stop-button').click();
  });
}
