// SPDX-License-Identifier: GPL-2.0-or-later
// Diagnostic observation of the real cold Worker, stopping only after its
// first startup. Does not assert GUI convergence or change any Newton update.
import fs from 'node:fs';
import path from 'node:path';
import { test, expect } from '@playwright/test';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('capture the real browser first coupled startup', async ({ page }) => {
  test.setTimeout(240000);
  const output = process.env.MSES_STARTUP_CAPTURE_OUTPUT;
  expect(output && !fs.existsSync(output)).toBeTruthy(); fs.mkdirSync(output, { recursive: true });
  const write = (name, v) => fs.writeFileSync(path.join(output, name + '.json'), JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x) + '\n');
  const hashes = numericalSourceHashes(['tests/browser/coupled-startup-capture.spec.js']);
  const original = fs.readFileSync('src/worker/solver.js', 'utf8');
  const anchor = 'const raw=solveCoupledStreamtubeAssembly(caseData,{';
  expect(original.split(anchor).length).toBe(2);
  const modified = original.replace(anchor, anchor + `
          onIterationCheckpoint:(checkpoint,details)=>self.postMessage({id,type:'diagnostic-checkpoint',checkpoint,details}),`)
    .replace("onStage:stage=>{progress.stage(stage);self.postMessage({id,type:'flow-stage',...stage});}",
      "onStage:stage=>{if(stage.startupAttempt===2)throw new Error('DIAGNOSTIC_FIRST_START_CAPTURED');progress.stage(stage);self.postMessage({id,type:'flow-stage',...stage});}");
  expect(modified).toContain('DIAGNOSTIC_FIRST_START_CAPTURED');
  fs.writeFileSync(path.join(output, 'worker-original.js'), original);
  fs.writeFileSync(path.join(output, 'worker-observed.js'), modified);
  write('start', { sourceHashes: hashes, scope: 'Exact GUI request and real Worker with one checkpoint observer and an explicit stop before startup2. No changed numerical controls. This diagnostic is not a cold-convergence certificate.' });
  await page.route('**/src/worker/solver.js', route => route.fulfill({ status: 200, contentType: 'application/javascript', body: modified }));
  await page.exposeFunction('recordStartupCapture', d => {
    if (d.type === 'diagnostic-checkpoint') {
      if (!fs.existsSync(path.join(output, 'initial-checkpoint.json'))) write('initial-checkpoint', { checkpoint: d.checkpoint, details: d.details });
      write('latest-checkpoint', { checkpoint: d.checkpoint, details: d.details });
    }
    else fs.appendFileSync(path.join(output, 'events.jsonl'), JSON.stringify(d) + '\n');
  });
  await page.addInitScript(() => {
    const OriginalWorker = Worker;
    window.captureFirst = { terminal: null, latest: null, request: null };
    window.Worker = class extends OriginalWorker {
      constructor(...args) { super(...args); this.addEventListener('message', ({ data }) => {
        if (!this.captureRequest?.caseData?.quadBoundaryLayers) return;
        if (data.type === 'diagnostic-checkpoint') window.captureFirst.latest = data;
        if (data.type === 'result' || data.type === 'error') window.captureFirst.terminal = data;
        if (['diagnostic-checkpoint', 'iteration', 'flow-stage', 'error'].includes(data.type)) window.recordStartupCapture(data);
      }); }
      postMessage(data, ...rest) { this.captureRequest = data;
        if (data.caseData?.quadBoundaryLayers) window.captureFirst.request = data;
        return super.postMessage(data, ...rest); }
    };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved', { timeout: 45000 });
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await page.locator('#preset').selectOption('rae2822-mses'); await page.locator('#solve-button').click();
  await expect.poll(() => page.evaluate(() => window.captureFirst.terminal !== null), { timeout: 180000 }).toBe(true);
  const saved = await page.evaluate(() => window.captureFirst); write('request', saved.request); write('terminal', saved.terminal);
  expect(saved.latest?.checkpoint.restart.options.ncrit).toBe(9);
  write('first-start-checkpoint', { checkpoint: saved.latest.checkpoint, details: saved.latest.details });
  const sourceChanges = changedSources(hashes); expect(sourceChanges).toEqual([]);
  write('report', { captured: true, sourceChanges, families: saved.latest.checkpoint.families,
    terminalMessage: saved.terminal.message, details: { stage: saved.latest.details.stage,
      startupAttempt: saved.latest.details.startupAttempt, iterationOffset: saved.latest.details.iterationOffset,
      final: saved.latest.details.history.at(-1) } });
});
