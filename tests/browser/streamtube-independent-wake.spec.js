import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directChannelConservation } from '../oracles/streamtube.js';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';

test('independent wake banks solve simultaneously in a browser worker and convert the actual default restart', async ({ page }) => {
  await page.goto('/');
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), wakeGeometry: 'independent-banks' };
  const fixture = JSON.parse(fs.readFileSync(new URL('../fixtures/default-coupled-ises-wake-fold.json', import.meta.url)));
  const data = await page.evaluate(({ input, fixture }) => new Promise((resolve, reject) => {
    const ises = new URL('/src/euler/streamtube-coupled-ises.js', location.href).href;
    const restart = new URL('/src/euler/tests/streamtube-independent-wake-restart.js', location.href).href;
    const source = `self.onmessage=async({data:{input,fixture:f}})=>{try{
      const {solveCoupledStreamtubeIses}=await import(${JSON.stringify(ises)});
      const {independentWakeBankRestart}=await import(${JSON.stringify(restart)});
      const snapshots=[];
      const result=solveCoupledStreamtubeIses(input,{edgeMatching:'section-velocity',maxIterations:12,tolerance:1e-10,
        stepAcceptance:'admissible',onMesh:m=>snapshots.push({nodes:m.nodes,iteration:m.iteration.iteration})});
      const c=independentWakeBankRestart(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:Float64Array.from(f.initialBL)});
      self.postMessage({result,snapshots,conversion:c.diagnostics,defaultStrict:c.system.admissible(c.system.initial),
        defaultSampled:c.system.admissible(c.system.initial,{requireConvex:false})});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({ input, fixture });
  }), { input, fixture });
  const r = data.result;
  expect(r.converged).toBe(true); expect(r.x.length).toBe(362); expect(r.mesh.quality.valid).toBe(true);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.solverInput.wakeGeometry).toBe('independent-banks');
  expect(data.snapshots.map(m => m.iteration)).toEqual(r.history.map(h => h.iteration));
  expect(data.snapshots.at(-1).nodes).toEqual(r.flow.nodes);
  expect(directStreamtubeVolumeGeometry(r.flow.nodes).valid).toBe(true);
  for (let g = 0; g < r.flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.flow.nodes[g], sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
  expect(data.conversion.newUnknowns).toBe(8278); expect(data.conversion.maxGeometryChange).toBeLessThan(1e-13);
  expect(data.conversion.gapResidual).toBeLessThan(1e-13); expect(data.defaultStrict).toBe(false); expect(data.defaultSampled).toBe(true);
  console.log(JSON.stringify({ families: r.families, updates: r.history.length - 1, conversion: data.conversion }));
});

test('the actual default coupled equations close from the corrected seed in a browser worker with explicit final-grid rejection', async ({ page }) => {
  await page.goto('/');
  const saved = JSON.parse(fs.readFileSync(new URL('../../docs/default-wake-initialization.json', import.meta.url)));
  const f = saved.cases.find(c => c.wakeGeometry === 'independent-banks').restart;
  const data = await page.evaluate(f => new Promise((resolve, reject) => {
    const ises = new URL('/src/euler/streamtube-coupled-ises.js', location.href).href;
    const source = `self.onmessage=async({data:f})=>{try{
      const {solveCoupledStreamtubeIses}=await import(${JSON.stringify(ises)});
      const snapshots=[];
      const result=solveCoupledStreamtubeIses(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:Float64Array.from(f.initialBL),
        maxIterations:12,tolerance:1e-8,stepAcceptance:'admissible',iterationGeometry:'ises-sampled',
        onMesh:m=>snapshots.push({nodes:m.nodes,iteration:m.iteration.iteration})});
      self.postMessage({result,snapshots});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(f);
  }), f);
  const r = data.result;
  expect(r.x.length).toBe(8278); expect(Math.max(...r.residual.map(Math.abs))).toBeLessThan(1e-8);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.converged).toBe(false); expect(r.mesh.quality.invalidCells).toEqual([463]); expect(r.mesh.initialization.flowSolved).toBe(false);
  expect(data.snapshots.map(m => m.iteration)).toEqual(r.history.map(h => h.iteration));
  expect(data.snapshots.at(-1).nodes).toEqual(r.flow.nodes); expect(data.snapshots.at(-1).nodes).not.toEqual(data.snapshots[0].nodes);
  expect(directStreamtubeVolumeGeometry(r.flow.nodes).valid).toBe(true);
  console.log(JSON.stringify({ defaultCoupled: { families: r.families, updates: r.history.length - 1, quality: r.mesh.quality, linear: r.linearDiagnostics } }));
});
