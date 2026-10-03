import { test, expect } from '@playwright/test';
test('grid smoothing WASM loads and converges in a browser module worker', async ({ page }) => {
  await page.route('**/smoothing-check', route => route.fulfill({ contentType: 'text/html', body: '<!doctype html><title>Smoothing</title>' }));
  await page.goto('/smoothing-check');
  const result = await page.evaluate(() => new Promise((resolve, reject) => {
    const base = location.origin;
    const source = `
      import {createEllipticStreamtubeGrid} from '${base}/src/geometry/elliptic-streamtube-grid.js';
      import {smoothPairedBoundaryGrid} from '${base}/src/geometry/paired-boundary-slor.js';
      try {
        const nx=16,nt=7,xi=Array.from({length:nx+1},(_,i)=>i/nx);
        const nodes=xi.map(u=>Array.from({length:nt+1},(_,j)=>({x:u,y:j/nt+.01*Math.sin(Math.PI*u)*Math.sin(Math.PI*j/nt)})));
        const system=createEllipticStreamtubeGrid({nodes,massFlows:Array(nt).fill(1),streamwiseCoordinates:xi,
          discretization:'giles-1985',lineLinearization:'full-metrics',lineGrouping:'boundary-pairs',lineSearch:'armijo',
          orthogonalBoundaryControl:{background:nodes.map(row=>row.map(()=>0))}});
        const result=smoothPairedBoundaryGrid(system,{maxSweeps:200});
        self.postMessage({backend:result.backend,replays:result.referenceReplays,converged:result.converged,valid:result.quality.valid,residual:result.residual});
      }catch(error){self.postMessage({error:error.message});}
    `;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); resolve(data); };
  }));
  expect(result.error).toBeUndefined(); expect(result.backend).toBe('wasm'); expect(result.replays).toBe(0);
  expect(result.converged).toBe(true); expect(result.valid).toBe(true); expect(result.residual).toBeLessThanOrEqual(1e-9);
});
