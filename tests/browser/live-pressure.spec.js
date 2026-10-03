import { test, expect } from '@playwright/test';

// Presentation tests only. Every Worker is held, including the page's initial
// request; no native Worker, mesh preparation or flow solve is started.
async function open(page) {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    window.heldPressureWorkers = [];
    window.heldPressureRequests = [];
    window.Worker = class {
      constructor() { window.heldPressureWorkers.push(this); }
      postMessage(request) {
        this.request = structuredClone(request);
        window.heldPressureRequests.push(this.request);
      }
      terminate() { this.stopped = true; }
    };
  });
  await page.goto('/');
  await page.waitForFunction(() => window.heldPressureRequests.length === 1);
  await page.locator('#stop-button').click();
  await page.locator('#flow-model').selectOption('streamtube-bl');
  return errors;
}

async function start(page) {
  const before = await page.evaluate(() => window.heldPressureRequests.length);
  await page.locator('#solve-button').click();
  await expect.poll(() => page.evaluate(() => window.heldPressureRequests.length)).toBe(before + 1);
  await expect(page.locator('#stop-button')).toBeVisible();
  return page.evaluate(() => ({
    index: window.heldPressureWorkers.length - 1,
    request: window.heldPressureRequests.at(-1),
  }));
}

async function send(page, message, workerIndex = -1) {
  await page.evaluate(({ message, workerIndex }) => {
    const worker = window.heldPressureWorkers.at(workerIndex);
    worker.onmessage({ data: { id: worker.request.id, ...message } });
  }, { message, workerIndex });
}

function frame(iteration, cp, mach = .2, targetMach = mach) {
  return {
    type: 'coefficients', iteration, mach, targetMach,
    coefficients: { cl: iteration / 10, cd: iteration / 1000, cm: -iteration / 100 },
    pressure: {
      referenceChord: 2, pressureKind: 'Synthetic live pressure',
      elements: [
        { name: 'Live main', cp: [0, .4, .8].map(x => ({ x, y: .1, cp })) },
        { name: 'Live flap', cp: [1, 1.3, 1.6].map(x => ({ x, y: -.1, cp })) },
      ],
    },
  };
}

const canvasImage = page => page.locator('#pressure-canvas').evaluate(canvas => canvas.toDataURL());

for (const mode of ['streamtube-bl', 'streamtube-grid']) {
  test(`${mode} displays provisional coefficients throughout Euler steps`, async ({ page }) => {
    const errors = await open(page);
    await page.locator('#flow-model').selectOption(mode);
    await start(page);
    await send(page, { type: 'flow-stage', stage: 'euler', mach: .47, targetMach: .74 });
    for (const iteration of [0, 1]) {
      await send(page, { type: 'iteration', iteration: { iteration, residual: .1 / (iteration + 1) } });
      await send(page, { type: 'coefficients', kind: 'euler-pressure', iteration, mach: .47, targetMach: .74,
        coefficients: { cl: .5 + iteration / 10, cd: .02, cm: -.1 } });
      await expect(page.locator('#overlay-cl')).toHaveText((.5 + iteration / 10).toFixed(4));
      await expect(page.locator('#overlay-cd')).toHaveText('0.020000');
      await expect(page.locator('#overlay-cm')).toHaveText('-0.1000');
      await expect(page.locator('#overlay-name')).toHaveText('Main element + flap');
      await expect(page.locator('#overlay-mach')).toHaveText('0.470');
      await expect(page.locator('#overlay-alpha')).toHaveText('4.00');
      await expect(page.locator('#overlay-cl')).toHaveText((.5 + iteration / 10).toFixed(4));
      await expect(page.locator('#overlay-cm')).toHaveText('-0.1000');
      await expect(page.locator('#overlay-cd')).toHaveText('0.020000');
      await expect(page.locator('#overlay-ld')).toHaveText(((.5 + iteration / 10) / .02).toFixed(1));
      await expect(page.locator('#overlay-re')).toHaveText(mode === 'streamtube-bl' ? '1.00e+6' : '—');
      await expect(page.locator('#overlay-ncr')).toHaveText(mode === 'streamtube-bl' ? '9.0' : '—');
      await expect(page.locator('#airfoil-overlay')).toHaveClass(/provisional/);
      await expect(page.locator('#airfoil-overlay')).not.toContainText('Provisional');

      await expect(page.locator('#coefficient-warning')).toContainText('CD excludes viscous drag');
      await expect(page.locator('#coefficient-warning')).toContainText('target 0.740');
      await expect(page.locator('#stop-button')).toBeVisible();
    }
    if (mode === 'streamtube-bl') {
      await send(page, frame(2, -.25));
      await expect(page.locator('#coefficient-warning')).not.toContainText('excludes viscous drag');
    }
    await send(page, { type: 'coefficients', iteration: 3, mach: .3, actualAlpha: 2.5, actualNcrit: 4,
      coefficients: { cl: .7, cm: -.2, cd: null } });
    await expect(page.locator('#overlay-ld')).toHaveText('—');
    await expect(page.locator('#overlay-cd')).toHaveText('—');
    await expect(page.locator('#overlay-alpha')).toHaveText('2.50');
    await expect(page.locator('#overlay-ncr')).toHaveText(mode === 'streamtube-bl' ? '4.0' : '—');
    await page.setViewportSize({ width: 390, height: 844 });
    await page.evaluate(() => scrollTo(0, 0));
    await expect(page.locator('.metrics')).toHaveCount(0);
    expect(await page.locator('.pressure-panel').evaluate(node => node.getBoundingClientRect().bottom <= innerHeight + 2)).toBe(true);
    expect(await page.locator('#airfoil-overlay').evaluate(node => {
      const overlay = node.getBoundingClientRect(), plot = node.closest('.canvas-wrap').getBoundingClientRect();
      return overlay.left >= plot.left && overlay.right <= plot.right && overlay.bottom <= plot.bottom;
    })).toBe(true);
    await page.locator('#settings-tab').click();
    await page.locator('#stop-button').click();
    await expect(page.locator('#overlay-cl')).toHaveText('—');
    expect(errors).toEqual([]);
  });
}

async function hoverPressure(page) {
  // The UI chooses the nearest current plotted point, so this checks its
  // actual pointer handler without reproducing the plot's pixel transform.
  await page.locator('#pressure-canvas').dispatchEvent('pointermove', { clientX: 200, clientY: 200 });
}

async function expectCleared(page) {
  await expect(page.locator('#cp-status')).toBeHidden();
  await expect(page.locator('#pressure-canvas')).not.toHaveAttribute('data-iteration');
  await expect(page.locator('#pressure-canvas')).not.toHaveAttribute('data-mach');
  await expect(page.locator('#cp-tooltip')).toBeHidden();
  await expect(page.locator('#cp-hover')).toHaveText('Hover to inspect');
  await expect(page.locator('.pressure-panel')).not.toHaveClass(/provisional/);
  await hoverPressure(page);
  await expect(page.locator('#cp-tooltip')).toBeHidden();
  // Inspect the two element colors, rather than comparing every axis/text
  // pixel across input changes which can also change the panel's layout.
  const curvePixels = await page.locator('#pressure-canvas').evaluate(canvas => {
    const pixels = canvas.getContext('2d').getImageData(0, 0, canvas.width, canvas.height).data;
    let count = 0;
    for (let i = 0; i < pixels.length; i += 4) {
      const main = pixels[i] === 137 && pixels[i + 1] === 229 && pixels[i + 2] === 208;
      const flap = pixels[i] === 214 && pixels[i + 1] === 186 && pixels[i + 2] === 125;
      if ((main || flap) && pixels[i + 3] > 0) count++;
    }
    return count;
  });
  expect(curvePixels).toBe(0);
}

for (const transition of ['fixed-trip', 'automatic']) {
  test(`${transition} cold quad solve draws successive live coefficient pressure and uses current hover data`, async ({ page }) => {
    const errors = await open(page);
    await page.locator('#quad-transition').selectOption(transition);
    const { request } = await start(page);
    expect(request.caseData).toMatchObject({ mach: .2, quadBoundaryLayers: true });
    expect(request.task).toBeUndefined(); // Ordinary cold route, not Mach continuation/refinement.
    const blank = await canvasImage(page);
    await send(page, frame(2, -.25));
    await expect(page.locator('#cp-status')).toContainText('Live');
    await expect(page.locator('#cp-status')).toContainText('iteration 2');
    await expect(page.locator('#cp-status')).toContainText('Mach 0.200');
    await expect(page.locator('#cp-status')).toContainText('provisional');
    await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-iteration', '2');
    await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-mach', '0.2');
    await expect(page.locator('.pressure-panel')).toHaveClass(/provisional/);
    await expect(page.locator('#overlay-cl')).toHaveText('0.2000');
    const first = await canvasImage(page);
    expect(first).not.toBe(blank);
    await hoverPressure(page);
    await expect(page.locator('#cp-tooltip')).toContainText('Live');
    await expect(page.locator('#cp-tooltip')).toContainText('Cp -0.250');

    await send(page, frame(3, .625));
    await expect(page.locator('#cp-tooltip')).toBeHidden();
    await expect(page.locator('#cp-hover')).toHaveText('Hover to inspect');
    await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-iteration', '3');
    expect(await canvasImage(page)).not.toBe(first);
    await hoverPressure(page);
    await expect(page.locator('#cp-tooltip')).toContainText('Cp 0.625');
    await expect(page.locator('#cp-tooltip')).not.toContainText('-0.250');
    await expect(page.locator('#cp-hover')).toHaveText('Cp 0.625');
    await expect(page.locator('#overlay-cl')).toHaveText('0.3000');
    await expect(page.locator('#status')).toHaveText('Solving');
    await expect(page.locator('#stop-button')).toBeVisible();
    await expect(page.locator('#export-button')).toBeDisabled();
    expect(errors).toEqual([]);
  });
}

test('cold pressure-only frames, Mach stages and unavailable messages clear the correct live data', async ({ page }) => {
  const errors = await open(page);
  await start(page);
  const blank = await canvasImage(page);
  const cold = frame(1, -.25);
  delete cold.coefficients;
  await send(page, { ...cold, type: 'pressure' });
  await expect(page.locator('#cp-status')).toContainText('iteration 1');
  await expect(page.locator('#cp-status')).toContainText('Mach 0.200');
  await expect(page.locator('#overlay-cl')).toHaveText('—');
  expect(await canvasImage(page)).not.toBe(blank);
  await send(page, { type: 'pressure-unavailable', mach: .2, message: 'Synthetic unavailable pressure' });
  await expectCleared(page);
  expect(await canvasImage(page)).toBe(blank);

  await send(page, frame(7, -.25, .47, .74));
  await expect(page.locator('#cp-status')).toContainText('Mach 0.470');
  await expect(page.locator('#cp-status')).toContainText('target 0.740');
  await send(page, { type: 'flow-stage', stage: 'coupled-mach', mach: .5375, targetMach: .74 });
  await expectCleared(page);
  expect(await canvasImage(page)).toBe(blank);
  await send(page, frame(1, .625, .5375, .74));
  await expect(page.locator('#cp-status')).toContainText('iteration 1');
  await expect(page.locator('#cp-status')).toContainText('Mach 0.537');
  await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-mach', '0.5375');
  await expect(page.locator('#cp-status')).toContainText('target 0.740');
  await hoverPressure(page);
  await expect(page.locator('#cp-tooltip')).toContainText('Cp 0.625');
  await send(page, { type: 'coefficients-unavailable', mach: .5375, message: 'Synthetic missing wall pressure' });
  await expectCleared(page);
  expect(await canvasImage(page)).toBe(blank);
  await expect(page.locator('#overlay-cl')).toHaveText('—');
  await expect(page.locator('#coefficient-warning')).toContainText('unavailable');
  await expect(page.locator('#stop-button')).toBeVisible();
  expect(errors).toEqual([]);
});

test('input changes, new jobs and Stop discard live Cp and ignore obsolete worker messages', async ({ page }) => {
  const errors = await open(page);
  const old = await start(page);
  await send(page, frame(4, -.25));
  await hoverPressure(page);
  await page.locator('#quad-reynolds').fill('1100000');
  await expect(page.locator('#status')).toHaveText('Inputs changed');
  await expectCleared(page);
  const current = await start(page);
  await expectCleared(page);
  await send(page, frame(5, .625));
  const currentImage = await canvasImage(page);
  await send(page, frame(99, -.4, .74, .74), old.index);
  await send(page, { type: 'pressure-unavailable', mach: .2 }, old.index);
  await send(page, { type: 'flow-stage', stage: 'coupled-mach', mach: .74, targetMach: .74 }, old.index);
  // Also reject a mismatched ID delivered on the current Worker object.
  await send(page, { ...frame(100, -.4), id: old.request.id }, current.index);
  await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-iteration', '5');
  await expect(page.locator('#overlay-cl')).toHaveText('0.5000');
  expect(await canvasImage(page)).toBe(currentImage);
  await page.locator('#stop-button').click();
  await expect(page.locator('#status')).toHaveText('Stopped');
  await expectCleared(page);
  await expect(page.locator('#overlay-cl')).toHaveText('—');
  await expect(page.locator('#export-button')).toBeDisabled();
  await send(page, frame(6, -.25), current.index);
  await expectCleared(page);
  await start(page);
  await expectCleared(page);
  expect(errors).toEqual([]);
});

test('a saved actual final pressure replaces the live provisional curve and a new run clears it', async ({ page }) => {
  const errors = await open(page);
  await page.locator('#preset').selectOption('rae2822-mses');
  await page.locator('#quad-mach').fill('.2');
  // The saved physical display root uses Ncrit4. Match the submitted case
  // rather than presenting that root as the form's default Ncrit9 result.
  await page.locator('#quad-ncrit').fill('4');
  await start(page);
  await send(page, frame(12, .625));
  await hoverPressure(page);
  const liveImage = await canvasImage(page);
  // This file already contains the actual Worker display result. Replay it
  // directly: no numerical reconstruction or presentation-adapter recompute.
  await page.evaluate(async () => {
    const response = await fetch('/docs/rae2822/real-browser-mses-coarse-cold32/state.json');
    if (!response.ok) throw new Error('Missing saved coupled display fixture.');
    const saved = await response.json();
    const worker = window.heldPressureWorkers.at(-1);
    worker.onmessage({ data: { id: worker.request.id, type: 'result', result: saved.result, elapsed: saved.elapsed } });
    window.finalPressureFixture = { cl: saved.result.cl, mach: saved.result.mach, converged: saved.result.converged };
  });
  const final = await page.evaluate(() => window.finalPressureFixture);
  expect(final.converged).toBe(true);
  await expect(page.locator('#status')).toHaveText('Euler/BL converged · research');
  await expect(page.locator('#stop-button')).toBeHidden();
  await expect(page.locator('#cp-status')).toBeHidden();
  await expect(page.locator('.pressure-panel')).not.toHaveClass(/provisional/);
  await expect(page.locator('#pressure-canvas')).not.toHaveAttribute('data-iteration');
  await expect(page.locator('#pressure-canvas')).toHaveAttribute('data-mach', String(final.mach));
  await expect(page.locator('#cp-tooltip')).toBeHidden();
  expect(await canvasImage(page)).not.toBe(liveImage);
  await hoverPressure(page);
  await expect(page.locator('#cp-tooltip')).toContainText('RAE 2822');
  await expect(page.locator('#cp-tooltip')).not.toContainText('Live');
  await expect(page.locator('#overlay-cl')).toHaveText(final.cl.toFixed(4));
  await expect(page.locator('#export-button')).toBeEnabled();
  await start(page);
  await expectCleared(page);
  await expect(page.locator('#export-button')).toBeDisabled();
  expect(errors).toEqual([]);
});


test('partial live loads retain provisional lift and moment when wake drag is unavailable', async ({ page }) => {
  const errors = await open(page);
  await start(page);
  const partial = frame(2, -.25, .47, .74);
  partial.coefficients.cd = null;
  await send(page, partial);
  await expect(page.locator('#overlay-cl')).toHaveText('0.2000');
  await expect(page.locator('#overlay-cm')).toHaveText('-0.0200');
  await expect(page.locator('#overlay-cd')).toHaveText('—');
  await expect(page.locator('#coefficient-warning')).toContainText('CD is unavailable');
  await expect(page.locator('#coefficient-warning')).toContainText('target 0.740');
  await send(page, frame(3, .625, .47, .74));
  await expect(page.locator('#overlay-cd')).toHaveText('0.003000');
  await expect(page.locator('#coefficient-warning')).not.toContainText('CD is unavailable');
  expect(errors).toEqual([]);
});

test('Mach contours display the current Euler field and remain zoomable', async ({ page }) => {
  const errors = await open(page);
  await start(page);
  await page.locator('#show-mach-contours').check();
  await expect(page.locator('#mach-contour-status')).toContainText('after an Euler flow update');
  const mesh = await page.evaluate(async () => {
    const { createStreamtubeBodySystem } = await import('/src/euler/streamtube-body.js');
    const { intrinsicBodyFixture } = await import('/tests/fixtures/intrinsic-body.js');
    const { streamtubeMeshSnapshot } = await import('/src/euler/streamtube-mesh-preview.js');
    const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 3 }));
    const flow = system.evaluate(system.initial);
    return streamtubeMeshSnapshot({system, nodes: flow.nodes, diagnostics: {}, flow, iteration: { iteration: 1 }});
  });
  await send(page, { type: 'mesh', mesh, stage: 'solving' });
  await page.locator('#show-streamlines').uncheck();
  await page.locator('#show-mesh').uncheck();
  await expect(page.locator('#mach-contour-status')).toContainText('ΔM = 0.05');
  await expect(page.locator('#mach-contour-status')).toContainText('unconverged');
  const before = await page.locator('#geometry-canvas').evaluate(c => c.toDataURL());
  await page.locator('#show-mach-contours').uncheck();
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).not.toBe(before);
  await page.locator('#show-mach-contours').check();
  await page.locator('#geometry-canvas').hover();
  await page.locator('.geometry-panel').screenshot({ path: '/tmp/mses-mach-contours.png' });
  await page.mouse.wheel(0, -200);
  await page.waitForTimeout(100);
  expect(await page.locator('#geometry-canvas').evaluate(c => c.toDataURL())).not.toBe(before);
  expect(errors).toEqual([]);
});
