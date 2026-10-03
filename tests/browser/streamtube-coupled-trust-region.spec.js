import { test, expect } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directChannelConservation } from '../oracles/streamtube.js';

async function solve(page, data) {
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  return page.evaluate(data => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createCoupledStreamtubeBody,solveCoupledStreamtubeBody}=await import(${JSON.stringify(module)});
      const s=createCoupledStreamtubeBody(data.input,data.options),meshes=[];
      const r=solveCoupledStreamtubeBody(s,{...data.controls,stepMethod:'dogleg',
        onMesh:m=>meshes.push({iteration:m.iteration,quality:m.quality,flowSolved:m.initialization.flowSolved})});
      self.postMessage({result:r,meshes});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(data);
  }), data);
}

test('browser coupled trust region closes four surfaces and two wakes with conservation and live iteration meshes', async ({ page }) => {
  const { result: r, meshes } = await solve(page, {
    input: intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    controls: { maxIterations: 16, tolerance: 1e-10, initialTrustRadius: .01 } });
  expect(r.converged).toBe(true); expect(r.stepMethod).toBe('dogleg');
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(meshes.map(m => m.iteration.iteration)).toEqual(r.history.slice(1).map(h => h.iteration));
  expect(meshes.every(m => m.quality.valid && !m.flowSolved)).toBe(true);
  for (let g = 0; g < r.flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.flow.nodes[g], sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
  console.log(JSON.stringify({ case: 'controlled two-element', unknowns: r.x.length, iterations: r.history.length - 1, families: r.families }));
});

test('browser WASM crosses the actual 7079-unknown default coupled material-trip event without claiming convergence', async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL('../fixtures/default-coupled-trip-event.json', import.meta.url)));
  const { result: r, meshes } = await solve(page, {
    input: fixture.input, options: { ...fixture.options, initialEuler: fixture.initialEuler, initialBL: fixture.initialBL },
    controls: { maxIterations: 1, tolerance: 1e-8, initialTrustRadius: 1 } });
  expect(r.x).toHaveLength(7079); expect(r.history).toHaveLength(2); expect(meshes).toHaveLength(1);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.boundaryLayer.stations).toHaveLength(295);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(r.transitionEvents).toBe(true); expect(r.history[1].activeChange).toBe(true);
  expect(r.history[1].meritComparable).toBe(false);
  expect(r.history[1].actualReduction).toBeNull(); expect(r.history[1].reductionRatio).toBeNull();
  expect(r.history[1].transitionChanges).toHaveLength(1);
  expect(r.history[1].transitionChanges[0]).toMatchObject({ body: 0, side: 'upper', from: 4, to: 5 });
  expect(r.boundaryLayer.surfaces[0].transition).toBe(5);
  expect(r.history[1].residual).toBeLessThan(r.history[0].residual);
  expect(r.converged).toBe(false); expect(r.status).toBe('unconverged');
  expect(r.mesh.quality.valid).toBe(true); expect(r.mesh.initialization.flowSolved).toBe(false);
  console.log(JSON.stringify({ case: 'retained default two-element', unknowns: r.x.length, history: r.history, quality: r.mesh.quality }));
});

test('browser WASM takes a projected coupled step at the actual default corner limit', async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL('../fixtures/default-coupled-geometry-limit.json', import.meta.url)));
  const { result: r, meshes } = await solve(page, {
    input: fixture.input, options: { ...fixture.options, initialEuler: fixture.initialEuler, initialBL: fixture.initialBL },
    controls: { maxIterations: 1, initialTrustRadius: 1 } });
  expect(r.x).toHaveLength(7079); expect(r.history).toHaveLength(2); expect(meshes).toHaveLength(1);
  // Model ranking selects Cauchy here; retain the full nonlinear reduction
  // and geometry requirements previously checked for projected Newton.
  expect(r.history[1].stepKind).toBe('projected-gradient'); expect(r.history[1].projection.converged).toBe(true);
  expect(r.history[1].actualReduction).toBeGreaterThan(80); expect(Math.abs(r.history[1].reductionRatio - 1)).toBeLessThan(.01);
  expect(r.history[1].residual).toBeLessThan(r.history[0].residual);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(meshes[0].quality.valid).toBe(true); expect(meshes[0].quality.minCornerSine).toBeGreaterThan(0);
  expect(meshes[0].flowSolved).toBe(false); expect(r.converged).toBe(false);
  console.log(JSON.stringify({ case: 'retained default corner constraint', step: r.history[1], quality: r.mesh.quality }));
});

for (const radius of [.25, 1]) test(`browser WASM corrects nonlinear coupled corner curvature at radius ${radius} and reruns material-trip preparation`, async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL('../fixtures/default-coupled-curvature-limit.json', import.meta.url)));
    const { result: r, meshes } = await solve(page, {
      input: fixture.input, options: { ...fixture.options, initialEuler: fixture.initialEuler, initialBL: fixture.initialBL },
      controls: { maxIterations: 1, initialTrustRadius: radius } });
    expect(r.x).toHaveLength(7079); expect(r.history).toHaveLength(2); expect(meshes).toHaveLength(1);
    const h = r.history[1]; expect(h.stepKind).toBe('projected-gradient-soc'); expect(h.correction.corrected).toBe(true);
    expect(h.trialRadius).toBe(radius); expect(h.correction.correctionNorm).toBeLessThan(1e-8 * radius);
    expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
    expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
    expect(meshes[0].quality.valid).toBe(true); expect(r.mesh.quality.valid).toBe(true);
    expect(r.mesh.quality.minCornerSine).toBeGreaterThan(1e-12); expect(r.converged).toBe(false);
    if (radius === .25) { expect(h.actualReduction).toBeGreaterThan(80); expect(h.reductionRatio).toBeGreaterThan(.99); }
    else {
      expect(h.activeChange).toBe(true); expect(h.actualReduction).toBeNull(); expect(h.reductionRatio).toBeNull();
      expect(h.transitionChanges.some(c => c.body === 0 && c.side === 'lower' && c.from === 6 && c.to === 5)).toBe(true);
    }
    console.log(JSON.stringify({ case: 'default nonlinear corner correction', radius, step: h, quality: r.mesh.quality }));
});

test('browser restores the saved coupled grid without requiring a valid undisplaced auxiliary grid', async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL('../fixtures/default-coupled-displaced-restart.json', import.meta.url)));
  const { result: r, meshes } = await solve(page, {
    input: fixture.input, options: { ...fixture.options, initialEuler: fixture.initialEuler, initialBL: fixture.initialBL },
    controls: { maxIterations: 0 } });
  expect(r.x).toHaveLength(7079); expect(r.history).toHaveLength(1); expect(meshes).toHaveLength(0);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.mesh.quality.valid).toBe(true); expect(r.converged).toBe(false);
  for (const [key, expected] of Object.entries(fixture.expectedFamilies)) expect(Math.abs(r.families[key] - expected)).toBeLessThan(1e-8);
});

test('browser WASM uses the corrected native wake derivative in a full default coupled Newton step', async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL('../fixtures/default-coupled-displaced-restart.json', import.meta.url)));
  const { result: r, meshes } = await solve(page, {
    input: fixture.input, options: { ...fixture.options, initialEuler: fixture.initialEuler, initialBL: fixture.initialBL },
    controls: { maxIterations: 1, initialTrustRadius: 2 } });
  expect(r.x).toHaveLength(7079); expect(r.history).toHaveLength(2); expect(meshes).toHaveLength(1);
  const step = r.history[1]; expect(step.stepKind).toBe('newton-ray'); expect(step.trialRadius).toBe(2);
  expect(step.actualReduction).toBeGreaterThan(13); expect(Math.abs(step.reductionRatio - 1)).toBeLessThan(.01);
  expect(step.residual).toBeLessThan(r.history[0].residual);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(meshes[0].quality.valid).toBe(true); expect(r.mesh.quality.valid).toBe(true); expect(r.converged).toBe(false);
  console.log(JSON.stringify({ case: 'corrected DILW derivative', step, quality: r.mesh.quality }));
});

test('browser WASM advances the retained default state on its active Hk branch without invalid derivative probes', async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL('../fixtures/default-coupled-active-limit.json', import.meta.url)));
  const { result: r, meshes } = await solve(page, {
    input: fixture.input, options: { ...fixture.options, initialEuler: fixture.initialEuler, initialBL: fixture.initialBL },
    controls: { maxIterations: 1, initialTrustRadius: 2 } });
  expect(r.x).toHaveLength(7079); expect(r.history).toHaveLength(2); expect(meshes).toHaveLength(1);
  expect(r.reason).toBe('iteration limit'); expect(r.history[1].residual).toBeLessThan(r.history[0].residual);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.mesh.quality.valid).toBe(true); expect(meshes[0].quality.valid).toBe(true); expect(r.converged).toBe(false);
  console.log(JSON.stringify({ case: 'exact active Hk branch', step: r.history[1], quality: r.mesh.quality }));
});

test('browser WASM projects the default coupled Newton direction at its compressible shape bound', async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL('../fixtures/default-coupled-shape-bound.json', import.meta.url)));
  const { result: r, meshes } = await solve(page, {
    input: fixture.input, options: { ...fixture.options, initialEuler: fixture.initialEuler, initialBL: fixture.initialBL },
    controls: { maxIterations: 1, initialTrustRadius: 16 } });
  expect(r.x).toHaveLength(7079); expect(r.history).toHaveLength(2); expect(meshes).toHaveLength(1);
  const h = r.history[1]; expect(h.stepKind).toBe('projected-newton'); expect(h.trialRadius).toBe(16);
  expect(h.projection.converged).toBe(true); expect(h.activeChange).toBeUndefined();
  expect(h.actualReduction).toBeGreaterThan(8.8); expect(Math.abs(h.reductionRatio - 1)).toBeLessThan(.02);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(meshes[0].quality.valid).toBe(true); expect(r.mesh.quality.minCornerSine).toBeGreaterThan(.009);
  expect(meshes[0].flowSolved).toBe(false); expect(r.converged).toBe(false); expect(r.reason).toBe('iteration limit');
  console.log(JSON.stringify({ case: 'projected Newton at default shape bound', step: h, quality: r.mesh.quality }));
});
