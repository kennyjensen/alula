import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';

for (const [mode, checkpoint] of [['isentropic', 'default-coupled-mass-interpolation.json'],
  ['momentum', 'default-positive-coupled-momentum-finish.json']])
test(`the positive default multielement coupled ${mode} root recovers a BL perturbation in a WASM worker`, async ({ page }) => {
  await page.goto('/');
  const f = JSON.parse(fs.readFileSync(new URL(`../../docs/${checkpoint}`, import.meta.url))).restart;
  expect(f.input.streamwiseMode).toBe(mode);
  const r = await page.evaluate(f => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const oracle = new URL('/tests/oracles/streamtube.js', location.href).href;
    const source = `self.onmessage=async({data:f})=>{try{
      const {createCoupledStreamtubeBody,solveCoupledStreamtubeBody}=await import(${JSON.stringify(module)});
      const {directChannelConservation}=await import(${JSON.stringify(oracle)});
      const bl=Float64Array.from(f.initialBL);
      for(let k=3;k<bl.length;k+=4)bl[k]*=1.0001;
      const system=createCoupledStreamtubeBody(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:bl});
      const snapshots=[];
      const r=solveCoupledStreamtubeBody(system,{maxIterations:4,tolerance:1e-8,
        onMesh:m=>snapshots.push({iteration:m.iteration,vertices:m.vertices})});
      const conservation=r.flow.nodes.map((nodes,g)=>{
        const c=directChannelConservation({nodes,sections:r.flow.sections.map(row=>row[g]),cells:r.flow.cells.map(row=>row[g])});
        return {total:c.total,maxLocal:c.maxLocal,internalCancellation:c.internalCancellation};
      });
      self.postMessage({converged:r.converged,reason:r.reason,families:r.families,history:r.history,quality:r.mesh.quality,
        nodes:r.flow.nodes,vertices:r.mesh.vertices,surfaces:r.boundaryLayer.surfaces.length,wakes:r.boundaryLayer.wakes.length,
        unknowns:r.x.length,linear:r.linearDiagnostics,snapshots,conservation});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(f);
  }), f);
  expect(r.unknowns).toBe(19859); expect(r.surfaces).toBe(4); expect(r.wakes).toBe(2);
  expect(r.converged, r.reason).toBe(true); expect(r.quality.valid).toBe(true);
  expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-8);
  expect(r.history[0].residual).toBeGreaterThan(1e-6); expect(r.history.length).toBeGreaterThan(1);
  expect(r.snapshots.length).toBe(r.history.length - 1);
  expect(r.snapshots.at(-1).vertices).toEqual(r.vertices);
  expect(directStreamtubeVolumeGeometry(r.nodes).valid).toBe(true);
  for (const c of r.conservation) {
    for (const k of mode === 'momentum' ? [0, 1, 2, 3] : [0, 3]) {
      expect(Math.abs(c.total[k])).toBeLessThan(1e-7);
      expect(Math.abs(c.maxLocal[k])).toBeLessThan(1e-7);
    }
    expect(Math.max(...c.internalCancellation.map(Math.abs))).toBeLessThan(1e-7);
  }
  console.log(JSON.stringify({ mode, unknowns: r.unknowns, updates: r.history.length - 1, families: r.families, quality: r.quality, linear: r.linear, conservation: r.conservation }));
});
