import { test, expect } from '@playwright/test';
import { readFileSync, writeFileSync } from 'node:fs';
import { streamtubeConvergenceDiagnostics } from '../../src/euler/streamtube-result.js';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('the explicit unsmoothed Mach 0.4 startup reports sonic capacity and retains the valid mesh before Newton', async ({page}) => {
  test.setTimeout(45000);
  const sourceHashes=numericalSourceHashes(['tests/browser/quad-failure-status.spec.js','index.html','src/ui/app.js','src/worker/solver.js']);
  const errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.addInitScript(()=>{
    window.sonicCheck={iterations:[]};const Original=window.Worker;
    window.Worker=class extends Original {
      postMessage(data){
        if(data.caseData?.flowModel==='streamtube-grid'){
          window.sonicCheck.input=data.caseData;
          this.addEventListener('message',({data:q})=>{
            if(q.type==='mesh')window.sonicCheck.mesh=q.mesh;
            if(q.type==='iteration')window.sonicCheck.iterations.push(q.iteration);
            if(q.type==='error')window.sonicCheck.failure=q;
            if(q.type==='result')window.sonicCheck.result=q.result;
          });
        }
        return super.postMessage(data);
      }
    };
  });
  await page.goto('/');await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#grid-elliptic').uncheck();
  await page.locator('#quad-mach').fill('.4');await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler initialization blocked',{timeout:30000});
  await expect(page.locator('#error-message')).toContainText('65.8% excess');
  await expect(page.locator('#error-message')).toContainText('Euler iterations have not started');
  await expect(page.locator('#error-message')).toContainText('shocks are not supported');
  const r=await page.evaluate(()=>window.sonicCheck);
  expect(r.input.mach).toBe(.4);expect(r.input.elements).toHaveLength(2);
  expect(r.input.gridIntervals).toBe(16);expect(r.input.gridTubes).toBe(7);expect(r.input.gridEllipticSmoothing).toBe(false);
  expect(r.failure.code).toBe('streamtube-sonic-capacity');
  expect(r.failure.diagnostics.section).toEqual({i:74,group:1,tube:0});
  expect(r.failure.diagnostics.capacityRatio).toBeCloseTo(1.6575579734226285,10);
  expect(r.iterations).toEqual([]);expect(r.result).toBeUndefined();
  expect(r.mesh.cells).toHaveLength(4587);expect(r.mesh.quality.valid).toBe(true);expect(r.mesh.flow).toBeUndefined();
  await expect(page.locator('#cl')).toHaveText('—');await expect(page.locator('#flow-legend')).toBeHidden();
  await expect(page.locator('#show-mesh')).toBeChecked();await expect(page.locator('#stop-button')).toBeHidden();
  await expect(page.locator('#mesh-status')).toContainText('4,587 cells');
  expect(errors).toEqual([]);expect(changedSources(sourceHashes)).toEqual([]);
  writeFileSync('/tmp/mses-mach04-gui-diagnostic.json',JSON.stringify({date:new Date().toISOString(),physicalAcceptance:false,sourceHashes,...r})+'\n');
});

// Presentation regression only: replay real, saved endpoints after applying
// the production diagnostic formatter. Numerical evidence is in the original
// cold-run reports and independent frozen audit, not in this replay.
test('quad status distinguishes stalled Newton from converged equations on a rejected grid', async ({ page }) => {
  const errors = []; page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    const Original = window.Worker;
    window.Worker = class extends Original {
      postMessage(data) {
        if (data.caseData?.flowModel === 'streamtube-grid' && window.replayQuad) {
          queueMicrotask(() => this.dispatchEvent(new MessageEvent('message', {
            data: { type: 'result', id: data.id, result: window.replayQuad, elapsed: 1 },
          })));
          return;
        }
        return super.postMessage(data);
      }
    };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  for (const [name, equationsClosed, smooth, rejectedCells] of [
    ['quad-32-11-unsmoothed', false, false, 1],
    ['quad-16-7-smoothed', true, true, 1],
    ['quad-32-11-smoothed', false, true, 2],
  ]) {
    const r = JSON.parse(readFileSync(`docs/${name}.json`)).result;
    Object.assign(r.diagnostics, streamtubeConvergenceDiagnostics(r.flow, r.mesh.quality));
    expect(r.diagnostics.residualConverged).toBe(equationsClosed);
    expect(r.diagnostics.gridValid).toBe(false);
    expect(r.diagnostics.reason).toContain(`${rejectedCells} non-convex or degenerate cell`);
    if (!equationsClosed) expect(r.diagnostics.reason).toContain(r.flow.lastRejectedStep.message);
    await page.evaluate(result => { window.replayQuad = result; }, r);
    await page.locator('#grid-elliptic').setChecked(smooth);
    await page.locator('#solve-button').click();
    await expect(page.locator('#status')).toHaveText(equationsClosed ? 'Euler converged · grid rejected' : 'Euler unconverged');
    await expect(page.locator('#error-message')).toContainText(equationsClosed ? 'Euler equations converged, but the final grid was rejected' : r.flow.lastRejectedStep.message);
    await expect(page.locator('#mesh-status')).toContainText(smooth ? 'SLOR smoothed' : 'SLOR not run');
    await expect(page.locator('#diagnostics')).toContainText(smooth ? 'completed' : 'not run');
    await expect(page.locator('#flow-iteration')).toContainText(equationsClosed ? 'equations converged · grid rejected · provisional' : 'unconverged · provisional');
    await expect(page.locator('#coefficient-warning')).toBeVisible();
    await expect(page.locator('#cl')).toHaveCSS('color', 'rgb(239, 189, 104)');
    await expect(page.locator('#export-button')).toBeEnabled();
  }
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  expect(errors).toEqual([]);
});
