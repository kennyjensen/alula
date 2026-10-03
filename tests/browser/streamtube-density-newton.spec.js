import { test, expect } from '@playwright/test';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directChannelConservation } from '../oracles/streamtube.js';

for (const motion of [
  { normalStencil: 'centered', stagnationMotion: 'interpolated' },
  { normalStencil: 'body-stations', stagnationMotion: 'walls-only' },
]) test(`additive-density Euler Newton (${motion.stagnationMotion}) closes a two-element case in a WASM worker and publishes accepted grids`, async ({ page }) => {
  await page.goto('/');
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), ...motion };
  const { result: r, meshes } = await page.evaluate(input => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-body.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createStreamtubeBodySystem,solveStreamtubeBody}=await import(${JSON.stringify(module)});
      const s=createStreamtubeBodySystem(data),meshes=[];
      const result=solveStreamtubeBody(s,{stepMethod:'density-newton',maxIterations:16,tolerance:1e-11,
        onMesh:m=>meshes.push({iteration:m.iteration,nodes:m.nodes})});
      self.postMessage({result,meshes});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(r.converged).toBe(true); expect(r.surfaces).toHaveLength(4);
  expect(r.projectedSteps).toBe(false); expect(r.secondOrderSteps).toBe(false);
  expect(r.diagnostics.residual).toBeLessThan(1e-11);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(meshes.map(m => m.iteration.iteration)).toEqual(r.history.slice(1).map(h => h.iteration));
  expect(meshes.at(-1).nodes).toEqual(r.nodes);
  expect(r.history.slice(1).every(h => h.stepKind === 'density-newton' && h.undampedUpdate)).toBe(true);
  for (let g = 0; g < r.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.nodes[g], sections: r.sections.map(row => row[g]), cells: r.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
  console.log(JSON.stringify({ ...motion, unknowns: r.x.length, iterations: r.history.length - 1, residual: r.diagnostics.residual }));
});
