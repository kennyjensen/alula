import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { test, expect } from '@playwright/test';
import { measureSurfaceGridSpacing, streamtubeSurfaceNodes } from '../../src/geometry/tests/surface-grid-spacing.js';

// Reconstruct logical rows directly from exported connectivity so the fixed
// boundary and spacing checks do not depend on the solver's layout helpers.
function meshRegions(mesh) {
  const { streamwiseSegments: nx, tubes } = mesh.initialization;
  assert.ok(Number.isInteger(nx) && nx >= 2);
  let offset = 0;
  return tubes.map(nt => {
    const indices = Array.from({ length: nx + 1 }, () => Array(nt + 1).fill(null));
    for (let i = 0; i < nx; i++) for (let j = 0; j < nt; j++) {
      const cell = mesh.cells[offset + i * nt + j];
      assert.equal(cell.length, 4);
      [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]].forEach(([a, b], k) => {
        if (indices[a][b] !== null) assert.equal(indices[a][b], cell[k]);
        indices[a][b] = cell[k];
      });
    }
    offset += nx * nt;
    return { nx, nt, nodes: indices.map(row => row.map(id => mesh.vertices[id])) };
  });
}

const repository = new URL('../../', import.meta.url);
const recordedSources = ['index.html', 'src/ui/app.js', 'src/euler/streamtube-result.js',
  'src/euler/streamtube-body-initializer.js', 'src/euler/streamtube-mesh-preview.js',
  'src/geometry/facing-outline-spacing.js', 'src/euler/streamtube-potential-guides.js',
  'src/euler/streamtube-elliptic-initializer.js', 'src/geometry/elliptic-streamtube-grid.js',
  'src/geometry/paired-boundary-slor.js', 'src/geometry/orthogonal-boundary-control.js',
  'src/geometry/boundary-stretch-control.js', 'src/geometry/passage-count-planner.js',
  'src/geometry/outer-streamtube-spacing.js', 'tests/browser/iset-slor-preview.spec.js'];
function sourceHashes() {
  return Object.fromEntries(recordedSources.map(path => [path,
    createHash('sha256').update(readFileSync(new URL(path, repository))).digest('hex')]));
}

test('accepted two-element GUI mesh uses paired SLOR with the same initial grid as the unchecked build', async ({ page }) => {
  test.setTimeout(60000);
  const errors = [], started = performance.now(), sources = sourceHashes();
  page.on('pageerror', error => errors.push(error.message));
  await page.addInitScript(() => {
    const OriginalWorker = window.Worker;
    window.isetEvents = []; window.isetRequests = [];
    window.Worker = class extends OriginalWorker {
      constructor(...args) {
        super(...args);
        this.addEventListener('message', ({ data }) => window.isetEvents.push(data));
      }
      postMessage(data, ...args) { window.isetRequests.push(data); super.postMessage(data, ...args); }
    };
  });
  await page.goto('/');
  await expect(page.locator('#status')).toContainText('Solved');
  await expect(page.locator('#preset')).toHaveValue('flap');
  await page.locator('#flow-model').selectOption('streamtube-grid');
  await expect(page.locator('#grid-crosslines')).toHaveCount(0);
  await expect(page.locator('#grid-smoothing-method')).toHaveCount(0);
  await expect(page.locator('#grid-elliptic')).toBeVisible();
  await page.locator('#mesh-button').click();
  await page.waitForFunction(() => window.isetEvents.some(e => e.type === 'mesh-ready' || e.type === 'error'),
    undefined, { timeout: 45000 });
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#mesh-status')).toContainText('SLOR not run');
  const unchecked = await page.evaluate(() => window.isetEvents.find(e => e.type === 'mesh').mesh);
  await page.locator('#grid-elliptic').check();
  await page.evaluate(() => { window.isetEvents = []; window.isetRequests = []; });
  await page.locator('#mesh-button').click();
  let terminalWaitError = null;
  try {
    await page.waitForFunction(() => window.isetEvents.some(e => e.type === 'mesh-ready' || e.type === 'error'),
      undefined, { timeout: 45000 });
  } catch (error) { terminalWaitError = error.message; }
  const { events, requests, status } = await page.evaluate(() => ({
    events: window.isetEvents, requests: window.isetRequests,
    status: { main: document.querySelector('#status')?.textContent,
      mesh: document.querySelector('#mesh-status')?.textContent },
  }));
  const snapshots = events.filter(e => e.type === 'mesh');
  const finalSourceHashes = sourceHashes();
  // Save the real outcome before terminal success and geometric assertions,
  // including rejected candidates or a timeout's latest available mesh.
  writeFileSync('/tmp/mses-iset-default-browser.json', JSON.stringify({
    case: 'actual potential-block GUI default with fixed boundaries',
    input: requests[0]?.caseData, requests, status, terminalWaitError, pageErrors: errors,
    seconds: (performance.now() - started) / 1000,
    initialMesh: snapshots[0]?.mesh, finalMesh: snapshots.at(-1)?.mesh,
    smoothing: snapshots.at(-1)?.mesh.initialization?.gridSmoothing,
    events: events.map(({ mesh, ...event }) => mesh ? { ...event,
      meshSummary: { vertices: mesh.vertices.length, cells: mesh.cells.length,
        quality: mesh.quality, gridSmoothing: mesh.initialization?.gridSmoothing } } : event),
    sourceHashes: sources, finalSourceHashes,
    changedSources: recordedSources.filter(path => sources[path] !== finalSourceHashes[path]),
    scope: 'Public default mesh regression and spacing screen; no physical airfoil acceptance.',
  }, null, 2) + '\n');
  await page.locator('.geometry-panel').screenshot({ path: '/tmp/mses-iset-slor-preview.png' });
  expect(terminalWaitError).toBeNull();
  await expect(page.locator('#status')).toHaveText('Mesh ready');
  await expect(page.locator('#mesh-status')).toContainText('SLOR smoothed');
  await expect(page.locator('#mesh-status')).not.toContainText('SLOR not run');
  expect(requests).toHaveLength(1);
  expect(requests[0]).toMatchObject({ meshOnly: true, caseData: { gridCrosslinePlacement: 'potential', gridEllipticSmoothing: true } });
  expect(events[0]).toMatchObject({ type: 'mesh', stage: 'initial' });
  expect(events.some(e => e.type === 'iteration' || e.type === 'error' || e.type === 'result')).toBe(false);
  expect(events.filter(e => e.type === 'mesh-ready')).toHaveLength(1);
  expect(events.some(e => e.stage === 'smoothing')).toBe(true);
  const initial = snapshots[0].mesh, mesh = snapshots.at(-1).mesh;
  expect(initial).toEqual(unchecked);
  expect(mesh.quality.valid).toBe(true);
  expect(mesh.initialization.flowSolved).toBe(false);
  const boundary = mesh.initialization.potentialCrosslines;
  expect(boundary.outerSpacing.every(side => side.indexInterpolation === 'pchip')).toBe(true);
  expect(boundary.outerSpacing.every(side => side.endSpacing === 'vinokur' && side.exactModernMsetX === false)).toBe(true);
  expect(mesh.cells).toEqual(initial.cells);
  expect(mesh.cells).toHaveLength(2919);
  expect(initial.initialization.cutStationSpacing).toBe('physical-x');
  expect(initial.initialization.gridSmoothing).toBeUndefined();
  expect(initial.initialization.gridRepair).toBeUndefined();
  expect(mesh.initialization.massFlows).toEqual(initial.initialization.massFlows);
  expect(mesh.initialization.gridSpacing).toEqual(initial.initialization.gridSpacing);
  const smoothing = mesh.initialization.gridSmoothing;
  expect(smoothing.boundaryControl).toBe('wall-angle');
  expect(smoothing.converged).toBe(true);
  expect(smoothing.retainedOriginal).not.toBe(true);
  expect(smoothing.regions).toHaveLength(3);
  for (const region of smoothing.regions) expect(region.coordinateEquations).toMatchObject({
    discretization: 'giles-1985', metricUpdate: 'each line', stationSpacing: 'prescribed nonuniform xi',
    lineLinearization: 'full-metrics', lineGrouping: 'boundary-pairs', lineSearch: 'armijo',
    boundaryConditions: { lower: 'fixed', upper: 'fixed' } });
  expect(smoothing.stationCoordinate).toMatchObject({ sharedAcrossPassages: true, exactModernMsetGauge: false,
    source: 'primary upper body/stagnation outline' });
  const xi = smoothing.stationCoordinate.xi;
  expect(xi).toHaveLength(mesh.initialization.streamwiseSegments + 1);
  expect(xi[0]).toBe(0); expect(xi.at(-1)).toBe(1);
  expect(xi.every((x, i) => !i || x > xi[i - 1])).toBe(true);
  expect(mesh.cells.every(cell => cell.length === 4)).toBe(true);

  const before = meshRegions(initial), after = meshRegions(mesh);
  const fixture = JSON.parse(readFileSync(new URL('tests/fixtures/default-grid-surface-spacing.json', repository)));
  const surfaceSpacing = fixture.surfaces.map(surface => measureSurfaceGridSpacing({ ...surface,
    points: streamtubeSurfaceNodes(mesh, surface) }));
  expect(surfaceSpacing.every(s => s.passed)).toBe(true);
  after.forEach((region, g) => region.nodes.forEach((row, i) => row.forEach((point, j) => {
    if (i === 0 || i === region.nx || j === 0 || j === region.nt)
      assert.deepEqual(point, before[g].nodes[i][j], `Fixed boundary moved in region ${g}, station ${i}, row ${j}`);
  })));
  const outerSpacing = [];
  for (const [stage, regions] of [['initial', before], ['final', after]]) {
    for (const [g, j, side] of [[0, 0, 'lower'], [regions.length - 1, regions.at(-1).nt, 'upper']]) {
      const row = regions[g].nodes.map(points => points[j]);
      const lengths = row.slice(1).map((p, i) => Math.hypot(p.x - row[i].x, p.y - row[i].y));
      assert.ok(lengths.every(length => length > 0 && Number.isFinite(length)));
      let maximumAdjacentRatio = 1, peakStation = null;
      for (let i = 1; i < lengths.length; i++) {
        const ratio = Math.max(lengths[i] / lengths[i - 1], lengths[i - 1] / lengths[i]);
        if (ratio > maximumAdjacentRatio) { maximumAdjacentRatio = ratio; peakStation = i; }
      }
      outerSpacing.push({ stage, side, maximumAdjacentRatio, peakStation });
      assert.ok(maximumAdjacentRatio <= 1.5,
        `${stage} ${side} farfield spacing ratio ${maximumAdjacentRatio} exceeds 1.5 at station ${peakStation}`);
    }
  }
  for (const cell of mesh.cells) for (let k = 0; k < 4; k++) {
    const [a, b, c] = [0, 1, 2].map(d => mesh.vertices[cell[(k + d) % 4]]);
    assert.ok((b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x) > 0,
      'Final quadrilateral has a nonpositive corner Jacobian.');
  }
  await expect(page.locator('#mesh-status')).toContainText('SLOR smoothed');
  await expect(page.locator('#cl')).toHaveText('—');
  expect(errors).toEqual([]);
  console.log(JSON.stringify({ case: 'actual potential-block GUI default with fixed boundaries',
    outerSpacing, surfaceSpacing: surfaceSpacing.map(s => s.lengthRatio),
    cells: mesh.cells.length, seconds: (performance.now() - started) / 1000,
    scope: 'Public default mesh regression and spacing screen; no physical airfoil acceptance.' }));
});
