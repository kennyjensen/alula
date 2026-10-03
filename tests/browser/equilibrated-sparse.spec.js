import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('browser WASM column equilibration solves a manufactured full multielement coupled Jacobian', async ({ page }) => {
  const sourceHashes = numericalSourceHashes(['tests/browser/equilibrated-sparse.spec.js',
    'scripts/validation/coupled-downstream-control.js', 'tests/fixtures/intrinsic-body.js']);
  await page.goto('/');
  const data = await page.evaluate(() => new Promise((resolve, reject) => {
    const code = `self.onmessage=async()=>{try{
      const {coupledDownstreamControl}=await import('${location.origin}/scripts/validation/coupled-downstream-control.js');
      const {solveEquilibratedSparse}=await import('${location.origin}/src/numerics/tests/equilibrated-sparse.js');
      const {sparseProduct}=await import('${location.origin}/src/numerics/sparse.js');
      const {mapped,result}=coupledDownstreamControl(),s=mapped.system,a=s.jacobian(result.x);
      const expected=Float64Array.from({length:s.n},(_,i)=>Math.sin(.37*i)),b=sparseProduct(a,expected);
      const r=solveEquilibratedSparse(a,b),ax=sparseProduct(a,r.x);
      self.postMessage({rootConverged:result.converged,unknowns:s.n,surfaces:s.bl.surfaces.length,wakes:s.bl.wakes.length,
        linear:{...r,x:undefined},forwardError:Math.max(...r.x.map((v,i)=>Math.abs(v-expected[i]))),
        independentRelativeResidual:Math.hypot(...ax.map((v,i)=>v-b[i]))/Math.hypot(...b)});
    }catch(error){self.postMessage({error:error.stack});}};`;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = event => { close(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({});
  }));
  expect(data.rootConverged).toBe(true); expect(data.unknowns).toBe(339);
  expect(data.surfaces).toBe(4); expect(data.wakes).toBe(2);
  expect(data.forwardError).toBeLessThan(2e-7); expect(data.independentRelativeResidual).toBeLessThan(1e-10);
  expect(data.linear.relativeResidual).toBeLessThan(1e-10); expect(data.linear.scaledRelativeResidual).toBeLessThan(1e-10);
  expect(changedSources(sourceHashes)).toEqual([]);
  fs.writeFileSync('docs/current-coupled-spatial-linear-browser.json', JSON.stringify({ date: new Date().toISOString(),
    physicalAcceptance: false, sourceHashes,
    scope: 'Actual worker/WASM linear solve of a manufactured RHS on a complete coupled two-element Jacobian; not a new nonlinear refinement or physical benchmark.',
    ...data }) + '\n');
  console.log(JSON.stringify(data));
});
