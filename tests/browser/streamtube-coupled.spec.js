import { test, expect } from '@playwright/test';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directChannelConservation } from '../oracles/streamtube.js';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';
import { readFileSync } from 'node:fs';

test('simultaneous two-element Euler and all BL/wake equations converge in a browser WASM worker', async ({ page }) => {
  await page.goto('/');
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const r = await page.evaluate(input => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createCoupledStreamtubeBody,solveCoupledStreamtubeBody}=await import(${JSON.stringify(module)});
      const system=createCoupledStreamtubeBody(data);
      self.postMessage(solveCoupledStreamtubeBody(system,{maxIterations:8,tolerance:1e-10}));
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(r.status).toBe('research-coupled-equations-converged'); expect(r.converged).toBe(true);
  for (const error of Object.values(r.families)) expect(error).toBeLessThan(1e-10);
  expect(r.history.length).toBeGreaterThan(1); expect(r.linearDiagnostics.solves).toBe(r.history.length - 1);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.boundaryLayer.stations).toHaveLength(31);
  expect(r.mesh.quality.valid).toBe(true); expect(r.mesh.initialization.flowSolved).toBe(true);
  expect(r.mesh.cells.every(c => c.length === 4)).toBe(true);
  expect(r.cl).toBeNull(); expect(r.cd).toBeNull(); expect(r.cm).toBeNull();
  console.log(JSON.stringify({ model: r.status, unknowns: r.x.length, iterations: r.history.length - 1,
    families: r.families, linear: r.linearDiagnostics, eulerDiagnostics: r.flow.diagnostics }));
  for (let g = 0; g < r.flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.flow.nodes[g], sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
});

test('the default coupled Newton trip transfer runs in a browser WASM worker', async ({ page }) => {
  const f = JSON.parse(readFileSync(new URL('../fixtures/default-newton-trip-crossing.json', import.meta.url)));
  await page.goto('/');
  const r = await page.evaluate(f => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const source = `self.onmessage=async({data:f})=>{try{
      const {createCoupledStreamtubeBody,solveCoupledStreamtubeBody}=await import(${JSON.stringify(module)});
      const system=createCoupledStreamtubeBody(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:Float64Array.from(f.initialBL)});
      const r=solveCoupledStreamtubeBody(system,{maxIterations:1});
      self.postMessage({history:r.history,families:r.families,unknowns:r.x.length,quality:r.mesh.quality,nodes:r.flow.nodes,
        surfaces:r.boundaryLayer.surfaces.length,wakes:r.boundaryLayer.wakes.length,linear:r.linearDiagnostics,
        pendingEvents:system.bl.activeTargets(r.x.subarray(0,system.ne)).filter(t=>t.from!==t.to)});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(f);
  }), f);
  expect(r.unknowns).toBe(8195); expect(r.history).toHaveLength(2);
  expect(r.history[1].step).toBeGreaterThanOrEqual(.0625); expect(r.history[1].activeChange).toBe(true);
  expect(r.history[1].meritComparable).toBe(false); expect(r.pendingEvents).toEqual([]);
  expect(r.surfaces).toBe(4); expect(r.wakes).toBe(2);
  expect(r.linear.solves).toBe(1); expect(r.linear.maxRelativeResidual).toBeLessThan(1e-10);
  expect(r.quality.valid).toBe(true); expect(directStreamtubeVolumeGeometry(r.nodes).valid).toBe(true);
  expect(r.families.boundaryLayer).toBeLessThan(f.families.boundaryLayer);
  console.log(JSON.stringify({ unknowns: r.unknowns, step: r.history[1].step, families: r.families }));
});
