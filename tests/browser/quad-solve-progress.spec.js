import { test, expect } from '@playwright/test';

// Synthetic Worker messages exercise the real page without running a solver.
test('method status keeps log recovery visible, records attempts and resets for a new request', async ({ page }) => {
  const errors = []; page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    window.Worker = class {
      postMessage(request) {
        window.emitSolveProgress = message => this.onmessage?.({ data: { ...message, id: request.id } });
      }
      terminate() {}
    };
  });
  await page.goto('/');
  await page.locator('#flow-model').selectOption('streamtube-bl');
  await page.locator('#solve-button').click();
  const emit = message => page.evaluate(message => window.emitSolveProgress(message), message);
  const stage = fields => emit({ type: 'flow-stage', ...fields });
  const iteration = fields => emit({ type: 'iteration', iteration: {
    stage: 'coupled', residual: .5, euler: .1, boundaryLayer: .5, edgeMatching: .2, ...fields,
  } });
  const current = page.locator('#solve-method-current');
  await expect(page.locator('#solution-status-title')).toHaveText('Solution Status');
  await expect(page.locator('#solution-status #export-button')).toHaveCount(1);
  await expect(page.locator('.workspace-heading')).toHaveCount(0);
  for (const id of ['status', 'run-meta', 'solve-methods', 'coefficient-warning', 'error-message',
    'streamline-status', 'mesh-status', 'mesh-update', 'flow-iteration', 'flow-change-readout', 'cp-status']) {
    await expect(page.locator(`#solution-status #${id}`)).toHaveCount(1);
    await expect(page.locator(`.geometry-panel #${id}`)).toHaveCount(0);
  }
  expect(await page.evaluate(() => {
    const pressure = document.querySelector('.pressure-panel').getBoundingClientRect();
    const status = document.querySelector('#solution-status').getBoundingClientRect();
    return status.top >= pressure.bottom;
  })).toBe(true);

  await stage({ stage: 'boundary-layer-initialization', startupAttempt: 1,
    boundaryLayerInitialization: { method: 'MRCHUE surfaces and ISET-style wake guess', wakeInitialization: 'iset-linear-shape', thicknessFactor: 1 },
    shearCoordinate: 'linear', hkFloorLinearization: 'exact' });
  await expect(current).toContainText('ISET wake guess');
  await stage({ stage: 'coupled', startupAttempt: 1, shearCoordinate: 'linear', hkFloorLinearization: 'exact' });
  await expect(current).toContainText('Ordinary Euler/BL startup');
  await expect(current).toContainText('linear shear');
  await emit({ type: 'mesh', stage: 'initial', mesh: {
    topology: 'intrinsic-quadrilateral-streamtubes',
    vertices: [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 1, y: 1 }, { x: 0, y: 1 }], cells: [[0, 1, 2, 3]],
    quality: { valid: true }, initialization: { tubes: [1] },
    iteration: { stage: 'coupled', iteration: 0, residual: .5 },
  } });
  await expect(page.locator('#mesh-update')).toHaveText('Starting Euler/BL grid · before the first Newton update');
  await page.locator('#show-mesh').uncheck();
  await page.locator('#show-streamlines').uncheck();
  await expect(page.locator('#solution-status #mesh-update')).toBeVisible();
  await expect(page.locator('#solution-status #mesh-status')).toBeVisible();

  await stage({ stage: 'boundary-layer-initialization', startupAttempt: 2, actualNcrit: 4, targetNcrit: 9,
    shearCoordinate: 'linear', hkFloorLinearization: 'native' });
  await stage({ stage: 'coupled', startupAttempt: 2, actualNcrit: 4, targetNcrit: 9 });
  await iteration({ iteration: 40, startupAttempt: 2 });
  await expect(current).toContainText('Lower-Ncrit recovery');
  await expect(current).toContainText('native Hk derivative');
  await stage({ stage: 'coupled', startupAttempt: 2, shearRecovery: { shearCoordinate: 'logarithmic' } });
  await iteration({ iteration: 41, startupAttempt: 2 });
  await expect(current).toContainText('Logarithmic shear recovery');
  await expect(page.locator('#run-meta')).toContainText('iteration 41');
  await page.locator('#solve-method-summary').click();
  const history = page.locator('#solve-method-history');
  await expect(history).toContainText('Ordinary Euler/BL startup');
  await expect(history).toContainText('ISET wake guess');
  await expect(history).toContainText('Lower-Ncrit recovery');
  await expect(history).toContainText('Logarithmic shear recovery');
  await stage({ stage: 'transition-refinement', startupAttempt: 2 });
  await expect(current).toContainText('Transition refinement');
  await expect(current).toContainText('logarithmic shear');
  await stage({ stage: 'coupled-ncrit-startup', startupAttempt: 2, actualNcrit: 5, targetNcrit: 9 });
  await expect(current).toContainText('Ncrit continuation');
  await expect(current).toContainText('Ncrit 5 → 9');
  await stage({ stage: 'coupled-mach', startupAttempt: 1, mach: .25, targetMach: .3 });
  await expect(current).toContainText('Mach continuation');
  await expect(current).toContainText('Mach 0.250 → 0.300');
  await stage({ stage: 'coupled-mach', mach: .72, targetMach: .74, gridLevel: 32,
    gridStrategy: 'target-mach-before-refinement' });
  await iteration({ stage: 'coupled-mach', iteration: 2, mach: .72, targetMach: .74,
    dissipation: { method: 'mses', mcrit: .81, targetMcrit: .99 } });
  await expect(current).toContainText('coarse-to-fine continuation');
  await expect(current).toContainText('MCRIT 0.810 → 0.990');
  await iteration({ stage: 'coupled-mach', iteration: 6, mach: .72, targetMach: .74,
    progress: { action: 'damped-recovery', cause: 'two-state-cycle', frozenMcrit: .81, stateChange: .001 },
    residualContext: { mcrit: .81, equationsChanged: false },
    prescribedResidual: { euler: .3, boundaryLayer: .6, edgeMatching: .2 } });
  await expect(current).toContainText('damped Newton recovery');
  await expect(current).toContainText('two state cycle');
  await expect(page.locator('#solution-status #run-meta')).toContainText('iteration 6');
  await stage({ stage: 'coupled-operating-point', mach: .72, targetMach: .74,
    actualAlpha: 2.2, targetAlpha: 2.68, stepMethod: 'combined' });
  await expect(current).toContainText('Mach–alpha continuation');
  await expect(current).toContainText('α 2.20° → 2.68°');
  await expect(current).toContainText('combined');
  await stage({ stage: 'coupled-grid-refinement', mach: .72, targetMach: .74,
    actualAlpha: 2.2, targetAlpha: 2.68, gridLevel: 64, gridStrategy: 'refine-retained-operating-point' });
  await expect(current).toContainText('refine retained operating point');
  await expect(current).toContainText('grid 64');
  await expect(current).not.toContainText('combined');
  await stage({ stage: 'coupled-grid-refinement', gridLevel: 64, gridStrategy: 'shock-triggered-refinement' });
  await expect(current).toContainText('shock-triggered refinement');
  await page.setViewportSize({ width: 390, height: 844 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.locator('#settings-tab').click();
  await page.locator('#stop-button').click();
  await expect(current).toContainText('Last method:');
  await expect(history).toContainText('Logarithmic shear recovery');
  await page.locator('#solve-button').click();
  await expect(page.locator('#solve-methods')).toBeHidden();
  await stage({ stage: 'euler', startupAttempt: 0 });
  await expect(current).toContainText('Euler precursor');
  await expect(page.locator('#solve-method-history li')).toHaveCount(1);
  await expect(history).not.toContainText('Logarithmic');
  await expect(current).not.toContainText('MCRIT');
  await expect(current).not.toContainText('damped Newton recovery');
  expect(errors).toEqual([]);
});
