// SPDX-License-Identifier: GPL-2.0-or-later
// Actual public cold default solve, with no worker/input/result substitution.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { chromium } from 'playwright';
import { directChannelConservation } from '../tests/oracles/streamtube.js';
import { numericalSourceHashes, changedSources } from './validation/provenance.js';
const started = performance.now();
const output = process.argv[2] ?? 'docs/default-euler-projected-browser.json';
const report = { date: new Date().toISOString(), physicalAcceptance: false,
  sourceHashes: numericalSourceHashes(['scripts/check-default-euler-browser.js', 'index.html', 'src/ui/app.js', 'tests/oracles/streamtube.js']),
  scope: 'Actual public default two-element cold quad Euler solve in Chromium/WASM at 390x844, smoothing off and ordinary iteration limit. No BL, mocked fields, prescribed restart or modified iteration controls.', failures: [] };
const browser = await chromium.launch({ headless: true });
try {
  const page = await browser.newPage({ viewport: { width: 390, height: 844 } });
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.addInitScript(() => {
    window.eulerCheck = { frames: [] };
    const OriginalWorker = window.Worker;
    window.Worker = class extends OriginalWorker {
      postMessage(data) {
        if (data.caseData?.flowModel === 'streamtube-grid') {
          window.eulerCheck.input = data.caseData;
          const receive = this.onmessage;
          this.onmessage = event => {
            receive(event); const q = event.data;
            if (q.type === 'mesh' && q.mesh.iteration) window.eulerCheck.frames.push({ iteration: q.mesh.iteration,
              quality: q.mesh.quality, label: document.getElementById('mesh-update').textContent,
              busy: !document.getElementById('stop-button').hidden });
            if (q.type === 'result') { window.eulerCheck.result = q.result; window.eulerCheck.done = true; }
            if (q.type === 'error') { window.eulerCheck.error = q.message; window.eulerCheck.done = true; }
          };
        }
        return super.postMessage(data);
      }
    };
  });
  await page.goto('http://localhost:5173/');
  await page.waitForFunction(() => document.getElementById('status').textContent.includes('Solved'), null, { timeout: 30000 });
  assert.equal(await page.locator('#preset').inputValue(), 'flap');
  assert.equal(await page.locator('#grid-intervals').inputValue(), '16');
  assert.equal(await page.locator('#grid-tubes').inputValue(), '7');
  assert.equal(await page.locator('#grid-elliptic').isChecked(), false);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#solve-button').click();
  await page.waitForFunction(() => window.eulerCheck.done, null, { timeout: 120000 });
  Object.assign(report, await page.evaluate(() => window.eulerCheck));
  report.ui = await page.evaluate(() => ({ status: document.getElementById('status').textContent,
    warning: document.getElementById('coefficient-warning').textContent,
    warningVisible: !document.getElementById('coefficient-warning').hidden,
    coefficientColor: getComputedStyle(document.getElementById('cl')).color,
    grid: document.getElementById('mesh-update').textContent }));
  assert.deepEqual(errors, []); assert.equal(report.error, undefined);
  const r = report.result;
  assert.equal(r.model, 'research-streamtube-euler'); assert.equal(r.boundaryLayer, null);
  assert.equal(r.flow.projectedSteps, true); assert.equal(r.solverSettings.projectedSteps, true);
  assert.equal(r.solverSettings.maxIterations, 40); assert.equal(r.solverInput.bodies.length, 2);
  assert.equal(report.frames.length, r.flow.history.length - 1);
  for (const [i, f] of report.frames.entries()) {
    assert.equal(f.iteration.iteration, i + 1); assert.ok(f.busy);
    assert.ok(f.label.includes(`Grid update ${i + 1}`));
  }
  report.restart = { input: r.solverInput, initialEuler: { x: r.flow.x, nodes: r.flow.nodes } };
  report.conservation = r.flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) }, r.conditions.gamma));
  if (!r.flow.converged) {
    report.failures.push(`Euler remains unconverged: ${r.flow.reason}`);
    assert.equal(r.status, 'unconverged'); assert.ok(report.ui.warning.includes('Unconverged'));
    assert.equal(report.ui.coefficientColor, 'rgb(239, 189, 104)');
  }
  if (!r.mesh.quality.valid) report.failures.push('Final mesh fails the published quality floor.');
  for (const c of report.conservation) for (const key of ['maxLocal', 'total', 'internalCancellation'])
    if (c[key].some(v => Math.abs(v) > 2e-9)) report.failures.push(`Independent ${key} balance is unfinished.`);
  report.browserChecksPassed = true;
} catch (error) { report.failures.push(error.message); report.error = error.stack; report.browserChecksPassed = false; }
finally { await browser.close(); }
report.seconds = (performance.now() - started) / 1000; report.changedSources = changedSources(report.sourceHashes);
if (report.changedSources.length) report.failures.push('Numerical/UI source changed during the browser solve.');
fs.writeFileSync(output, JSON.stringify(report,
  (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v, 2) + '\n');
console.log(JSON.stringify({ seconds: report.seconds, browserChecksPassed: report.browserChecksPassed,
  diagnostics: report.result?.flow?.diagnostics, ui: report.ui, quality: report.result?.mesh.quality, failures: report.failures, error: report.error }));
if (report.failures.length) process.exitCode = 1;
