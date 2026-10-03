import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { resolveStreamtubeGridControls } from '../src/euler/streamtube-result.js';
import { repeatedGridStepLimit } from '../src/euler/streamtube-iteration-progress.js';
const input = { gridIntervals: 128, gridTubes: 24, mach: .74, alpha: 2.68, eulerIsmom: 4,
  quadBoundaryLayers: false, gridSurfaceSpacing: 'automatic', gridEllipticSmoothing: true };
test('automatic fine-grid preparation preserves operating point and explicit end counts', () => {
  const before = structuredClone(input);
  assert.deepEqual(resolveStreamtubeGridControls(input), { surfaceSpacing: 'source', normalSpacing: 'automatic',
    inletIntervals: 64, outletIntervals: 64, balanced: true });
  assert.deepEqual(resolveStreamtubeGridControls({ ...input, gridInletIntervals: 32, gridOutletIntervals: 48,
    gridStagnationAspectRatio: 3 }), { surfaceSpacing: 'source', normalSpacing: 'stagnation',
    inletIntervals: 32, outletIntervals: 48, balanced: true });
  assert.deepEqual(input, before);
});
test('explicit distributions and coarser or coupled grids retain their preparation', () => {
  for (const gridSurfaceSpacing of ['supplied', 'curvature', 'source'])
    assert.deepEqual(resolveStreamtubeGridControls({ ...input, gridSurfaceSpacing }), {
      surfaceSpacing: gridSurfaceSpacing, normalSpacing: 'supplied', inletIntervals: undefined, outletIntervals: undefined, balanced: false });
  for (const override of [{ gridIntervals: 64 }, { quadBoundaryLayers: true }, { mach: .2 }])
    assert.deepEqual(resolveStreamtubeGridControls({ ...input, ...override }), {
      surfaceSpacing: 'supplied', normalSpacing: 'supplied', inletIntervals: undefined, outletIntervals: undefined, balanced: false });
});
const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/rae128-grid-stall-history.json', import.meta.url))).history;
test('captured RAE cell collapse is detected before further microscopic steps', () => {
  const before = structuredClone(fixture), result = repeatedGridStepLimit(fixture, 1e-10);
  assert.deepEqual(result.cell, { group: 1, i: 95, tube: 1, corner: 0 });
  assert.deepEqual(result.iterations, [9, 10, 11]);
  assert.deepEqual(fixture, before);
});
test('grid monitor does not stop improving, converged, changing-law or unrelated limited steps', () => {
  assert.equal(repeatedGridStepLimit(fixture.slice(1), 1e-10), null);
  for (const edit of [h => { h.at(-1).residual *= .9; }, h => { h.at(-1).residual = 1e-11; },
    h => { h.at(-1).step = .1; }, h => { h.at(-1).dissipation.mcrit = .99; },
    h => { h.at(-1).rejections = []; }, h => { h.at(-1).rejections[0].diagnostics.cell.i++; }]) {
    const h = structuredClone(fixture); edit(h); assert.equal(repeatedGridStepLimit(h, 1e-10), null);
  }
});

test('transonic viscous precursors share the harmonic seed without changing grid controls', async () => {
  const { resolveStreamtubeStartup } = await import('../src/euler/streamtube-result.js');
  assert.equal(resolveStreamtubeStartup({ ...input, gridIntervals: 64, quadBoundaryLayers: true }), 'harmonic-shock');
  for (const gridIntervals of [8,16,32])
    assert.equal(resolveStreamtubeStartup({ ...input, gridIntervals, quadBoundaryLayers: true }), 'harmonic-shock');
  assert.equal(resolveStreamtubeStartup({ ...input, gridIntervals: 32, quadBoundaryLayers: false }), 'standard');
  for (const override of [{ mach: .2 }, { gridEllipticSmoothing: false }, { eulerIsmom: 2 }])
    assert.equal(resolveStreamtubeStartup({ ...input, quadBoundaryLayers: true, ...override }), 'standard');
});
