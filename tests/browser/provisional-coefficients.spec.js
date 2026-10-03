import { test, expect } from '@playwright/test';

async function open(page) {
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => {
          if (data.type === 'result') window.lastCoefficientResult = data.result;
        });
      }
      postMessage(data) {
        // Bound the real numerical worker via its iteration control. No
        // replacement field, mesh, equation or response is used in this test.
        if (window.coefficientTestLimit !== undefined && data.caseData)
          data = { ...data, caseData: { ...data.caseData, maxIterations: window.coefficientTestLimit } };
        if (window.watchLiveMesh && data.caseData) {
          const receive = this.onmessage;
          this.onmessage = event => {
            receive(event);
            // Observe after the actual UI handler, not from an earlier
            // EventTarget listener's microtask checkpoint.
            const message = event.data;
            if (message.type !== 'mesh') return;
            window.liveMeshFrames.push({ mesh: message.mesh,
              image: document.getElementById('geometry-canvas').toDataURL(),
              label: document.getElementById('mesh-status').textContent,
              progress: document.getElementById('mesh-update').textContent,
              flowProgress: document.getElementById('flow-iteration').textContent,
              flowLegendVisible: !document.getElementById('flow-legend').hidden,
              speedChange: document.getElementById('flow-change-readout').textContent,
              busy: !document.getElementById('stop-button').hidden });
            if (message.mesh.iteration?.iteration === 2) document.getElementById('stop-button').click();
          };
        }
        return super.postMessage(data);
      }
    };
  });
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  await page.locator('#preset').selectOption('single'); await expect(page.locator('#status')).toContainText('Solved');
  return errors;
}

for (const [preset, smooth] of [['single', false], ['flap', false], ['flap', true]]) test(`quad ${preset}${smooth ? ' smoothed' : ''} draws every accepted grid while busy and Stop retains the second iterate`, async ({ page }) => {
  test.setTimeout(90000);
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#preset').selectOption(preset);
  await page.locator('#grid-intervals').selectOption(preset === 'single' ? '8' : '16');
  await page.locator('#alpha-number').fill(preset === 'single' ? '2' : '4');
  await page.locator('#grid-elliptic').setChecked(smooth);
  await page.evaluate(() => { window.watchLiveMesh = true; window.liveMeshFrames = []; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Stopped', { timeout: smooth ? 75000 : 35000 });
  const frames = await page.evaluate(() => window.liveMeshFrames), moving = frames.filter(f => f.mesh.iteration?.iteration > 0);
  expect(moving.map(f => f.mesh.iteration.iteration)).toEqual([1, 2]);
  const initial = frames.filter(f => !f.mesh.iteration?.iteration).at(-1);
  expect(initial.mesh.flow.iteration).toBe(0);
  expect(initial.flowProgress).toContain('Initial flow estimate');
  expect(initial.flowLegendVisible).toBe(true);
  expect(initial.busy).toBe(true);
  expect(initial.mesh.initialization.gridSmoothing?.converged ?? false).toBe(smooth);
  let previous = initial;
  for (const f of moving) {
    expect(f.busy).toBe(true); expect(f.label).toContain(`Euler iteration ${f.mesh.iteration.iteration}`);
    expect(f.mesh.vertices).not.toEqual(previous.mesh.vertices); expect(f.image).not.toBe(previous.image);
    const movement = Math.max(...f.mesh.vertices.map((p, i) => Math.hypot(p.x - previous.mesh.vertices[i].x, p.y - previous.mesh.vertices[i].y)));
    expect(f.mesh.iteration.maximumNodeMovement).toBe(movement);
    const fromStart = Math.max(...f.mesh.vertices.map((p, i) => Math.hypot(p.x - initial.mesh.vertices[i].x, p.y - initial.mesh.vertices[i].y)));
    expect(f.mesh.iteration.maximumNodeMovementFromInitial).toBe(fromStart);
    expect(f.progress).toContain(`Grid update ${f.mesh.iteration.iteration}`);
    expect(f.progress).toContain(`${(movement / f.mesh.iteration.movementReferenceChord).toExponential(2)} c_ref this step`);
    expect(f.progress).toContain(`${(fromStart / f.mesh.iteration.movementReferenceChord).toExponential(2)} c_ref from start`);
    expect(f.label).toContain(`max move ${(movement / f.mesh.iteration.movementReferenceChord).toExponential(1)} c_ref`);
    expect(f.mesh.initialization.flowSolved).toBe(false);
    expect(f.flowLegendVisible).toBe(true);
    expect(f.flowProgress).toContain(`Euler iteration ${f.mesh.iteration.iteration}`);
    expect(f.flowProgress).toContain('solving · provisional');
    expect(f.mesh.flow.residual).toBe(f.mesh.iteration.residual);
    const maximumSpeedChange = Math.max(...f.mesh.flow.lines.flatMap((line, k) => line.speedRatios.map((q, i) => Math.abs(q - previous.mesh.flow.lines[k].speedRatios[i]))));
    expect(f.mesh.flow.maximumSpeedChangeFromPrevious).toBe(maximumSpeedChange);
    expect(f.speedChange).toContain(`${(maximumSpeedChange * 100).toFixed(3)}% of U∞`);
    expect(f.mesh.flow.lines.flatMap(line => line.speedChanges)).toEqual(f.mesh.flow.lines.flatMap((line, k) =>
      line.speedRatios.map((q, i) => q - initial.mesh.flow.lines[k].speedRatios[i])));
    expect(f.mesh.flow.lines.map(l => l.speedRatios)).not.toEqual(previous.mesh.flow.lines.map(l => l.speedRatios));
    expect(f.mesh.initialization.gridSpacing).toEqual(initial.mesh.initialization.gridSpacing);
    previous = f;
  }
  await expect(page.locator('#mesh-status')).toContainText('Euler iteration 2');
  await expect(page.locator('#mesh-status')).toContainText('stopped');
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).toBe(moving.at(-1).image);
  await expect(page.locator('#mesh-update')).toHaveAttribute('data-iteration', '2');
  await expect(page.locator('#flow-iteration')).toContainText('stopped · provisional');
  await expect(page.locator('#flow-iteration')).toHaveCSS('color', 'rgb(239, 189, 104)');
  await page.locator('#show-mesh').uncheck();
  await expect(page.locator('#flow-legend')).toBeVisible();
  const flowOnly = await page.locator('#geometry-canvas').evaluate(c => c.toDataURL());
  await page.locator('#flow-color-mode').selectOption('change');
  await expect(page.locator('#flow-view-note')).toContainText('corresponding tube sections');
  await expect(page.locator('#flow-speed-ticks')).toContainText('≤−0.10');
  const changesOnly = await page.locator('#geometry-canvas').evaluate(c => c.toDataURL());
  expect(changesOnly).not.toBe(flowOnly);
  await page.locator('#flow-color-mode').selectOption('speed');
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).toBe(flowOnly);
  await page.locator('#show-streamlines').uncheck();
  await expect(page.locator('#flow-legend')).toBeHidden();
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).not.toBe(flowOnly);
  await page.locator('#show-streamlines').check();
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).toBe(flowOnly);
  await page.locator('#show-mesh').check();
  await page.locator('#compare-mesh').check();
  await expect(page.locator('#mesh-update')).toContainText('gold dashed: starting grid; blue: current grid');
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).not.toBe(moving.at(-1).image);
  await page.locator('#compare-mesh').uncheck();
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).toBe(moving.at(-1).image);
  if (smooth) {
    await page.locator('#compare-mesh').check();
    await page.setViewportSize({ width: 390, height: 844 });
    expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
    await page.locator('.geometry-panel').screenshot({ path: '/tmp/mses-live-grid-comparison-mobile.png' });
  }
  await expect(page.locator('#solve-button')).toBeEnabled();
  console.log(JSON.stringify({ preset, smooth, cells: initial.mesh.cells.length,
    updates: moving.map(f => ({ iteration: f.mesh.iteration.iteration, residual: f.mesh.iteration.residual,
      maximumNodeMovement: f.mesh.iteration.maximumNodeMovement,
      maximumNodeMovementFromInitial: Math.max(...f.mesh.vertices.map((p, i) => Math.hypot(p.x - initial.mesh.vertices[i].x, p.y - initial.mesh.vertices[i].y))) })) }));
  expect(errors).toEqual([]);
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await expect(page.locator('#flow-legend')).toBeVisible();
  await page.locator('#show-mesh').uncheck();
  await page.locator('.geometry-panel').screenshot({ path: `/tmp/mses-live-flow-${preset}-${smooth ? 'smooth' : 'raw'}-mobile.png` });
  await page.locator('#flow-color-mode').selectOption('change');
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.locator('.geometry-panel').screenshot({ path: `/tmp/mses-live-flow-changes-${preset}-mobile.png` });
  await page.locator('#alpha-number').fill('3');
  await expect(page.locator('#flow-legend')).toBeHidden();
});

async function checkProvisional(page) {
  await expect(page.locator('#coefficient-warning')).toBeVisible();
  await expect(page.locator('#coefficient-warning')).toContainText('Unconverged');
  await expect(page.locator('#error-message')).toContainText('did not converge');
  const result = await page.evaluate(() => window.lastCoefficientResult);
  expect(result.status).toBe('unconverged'); expect(result.coefficientStatus).toBe('unconverged');
  for (const k of ['cl', 'cm', 'cd']) {
    expect(Number.isFinite(result[k])).toBe(true);
    await expect(page.locator(`#${k}`)).toHaveText(result[k].toFixed(k === 'cd' ? 6 : 4));
    await expect(page.locator(`#${k}`)).toHaveCSS('color', 'rgb(239, 189, 104)');
    await expect(page.locator(`#${k}`)).toHaveAttribute('title', /Unconverged/);
  }
  const download = page.waitForEvent('download'); await page.locator('#export-button').click();
  const stream = await (await download).createReadStream(); let text = ''; for await (const part of stream) text += part;
  const exported = JSON.parse(text).result;
  expect(exported.status).toBe('unconverged'); expect(exported.coefficientStatus).toBe('unconverged');
  expect(exported.cl).toBe(result.cl); expect(exported.cd).toBe(result.cd);
  expect(exported.warnings.join(' ')).toContain('provisional');
  return result;
}

test('quad last-iterate pressure coefficients display in amber, export their status and reset on input/mesh changes', async ({ page }) => {
  test.setTimeout(90000);
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await page.locator('#grid-intervals').selectOption('8'); await page.locator('#alpha-number').fill('2');
  await page.evaluate(() => { window.coefficientTestLimit = 1; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler unconverged');
  const r = await checkProvisional(page);
  await expect(page.locator('#flow-iteration')).toContainText('unconverged · provisional');
  expect(r.mesh.flow.iteration).toBe(r.diagnostics.iterations);
  expect(r.mesh.flow.residual).toBe(r.diagnostics.equationResidual);
  const recorded = r.flow.history.at(-1).iteration;
  expect(recorded).toBe(1);
  expect(r.mesh.flow.maximumSpeedChangeFromPrevious).toBeGreaterThan(0);
  await expect(page.locator('#flow-change-readout')).not.toContainText('0.000%');
  expect(r.coefficients.dragKind).toBe('pressure'); expect(r.coefficients.physicalValidation).toBe(false);
  await expect(page.locator('#drag-label')).toHaveText('Pressure drag');
  await expect(page.locator('#drag-description')).toContainText('no viscous drag');
  await expect(page.locator('.metrics')).toBeVisible();
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.screenshot({ path: '/tmp/mses-provisional-coefficients-mobile.png', fullPage: false });
  await page.locator('#mesh-button').click();
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#flow-legend')).toBeHidden();
  await expect(page.locator('#cl')).toHaveText('—'); await expect(page.locator('#coefficient-warning')).toBeHidden();
  await expect(page.locator('#export-button')).toBeDisabled();
  await page.locator('#alpha-number').fill('0');
  await page.evaluate(() => { delete window.coefficientTestLimit; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Euler converged · research', { timeout: 35000 });
  await expect(page.locator('#flow-iteration')).toContainText('equations converged · research, unvalidated');
  await expect(page.locator('#coefficient-warning')).toBeHidden();
  await expect(page.locator('.metric.provisional')).toHaveCount(0);
  expect(await page.evaluate(() => window.lastCoefficientResult.coefficientStatus)).toBe('research-unvalidated');
  await page.locator('#reference').fill('0'); await page.locator('#solve-button').click();
  await expect(page.locator('#error-message')).toContainText('positive reference chord');
  await expect(page.locator('#cl')).toHaveText('—'); await expect(page.locator('#export-button')).toBeDisabled();
  expect(errors).toEqual([]);
});

test('coupled BL worker exposes unfinished loads through the same warning/color/export UI', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#flow-model').selectOption('coupled');
  await page.evaluate(() => { window.coefficientTestLimit = 1; });
  await page.locator('#solve-button').click();
  await expect(page.locator('#status')).toHaveText('Flow unconverged');
  await checkProvisional(page);
  await expect(page.locator('#drag-label')).toHaveText('Drag coefficient');
  await expect(page.locator('#drag-description')).toContainText('Squire–Young');
  await page.locator('#alpha-number').fill('3');
  await expect(page.locator('#coefficient-warning')).toBeHidden();
  await expect(page.locator('.metric.provisional')).toHaveCount(0);
  await expect(page.locator('#cl')).toHaveText('—'); await expect(page.locator('#export-button')).toBeDisabled();
  expect(errors).toEqual([]);
});
