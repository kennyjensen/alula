import { test, expect } from '@playwright/test';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directChannelConservation } from '../oracles/streamtube.js';

test('two-element prescribed displacement and separated wake banks solve in a real browser WASM worker', async ({ page }) => {
  await page.goto('/');
  const input = intrinsicBodyFixture({ elements: 2 });
  input.displacement = { surfaces: input.bodies.map(b => ({
    upper: Array(b.trailingIndex - b.leadingIndex + 1).fill(.0003), lower: Array(b.trailingIndex - b.leadingIndex + 1).fill(.0003) })),
    wakes: input.bodies.map(b => Array(input.outerLower.length - 1 - b.trailingIndex).fill(.0006)) };
  const { result, mesh } = await page.evaluate(input => new Promise((resolve, reject) => {
    const solver = new URL('/src/euler/streamtube-body.js', location.href).href;
    const preview = new URL('/src/euler/streamtube-mesh-preview.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createStreamtubeBodySystem,solveStreamtubeBody}=await import(${JSON.stringify(solver)});
      const {streamtubeMeshSnapshot}=await import(${JSON.stringify(preview)});
      const system=createStreamtubeBodySystem(data),result=solveStreamtubeBody(system,{tolerance:1e-11,maxIterations:20});
      self.postMessage({result,mesh:streamtubeMeshSnapshot({system,nodes:result.nodes})});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(result.converged).toBe(true); expect(result.diagnostics.residual).toBeLessThan(1e-11);
  expect(result.linearBackend).toBe('klu'); expect(result.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(result.surfaces).toHaveLength(4); expect(result.displacement.wakes).toHaveLength(2);
  expect(result.formulation).toContain('no BL coupling'); expect(result.forceStatus).toContain('not solid-wall aerodynamic forces');
  expect(mesh.quality.valid).toBe(true); expect(mesh.cells.every(c => c.length === 4)).toBe(true);
  for (let b = 0; b < input.bodies.length; b++) for (let i = input.bodies[b].trailingIndex + 1; i < input.outerLower.length; i++) {
    const lower = result.nodes[b][i].at(-1), upper = result.nodes[b + 1][i][0];
    expect(Math.hypot(upper.x - lower.x, upper.y - lower.y)).toBeCloseTo(.0006, 12);
  }
  for (let g = 0; g < result.nodes.length; g++) {
    const c = directChannelConservation({ nodes: result.nodes[g], sections: result.sections.map(row => row[g]), cells: result.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'external', 'internalCancellation']) expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
});
