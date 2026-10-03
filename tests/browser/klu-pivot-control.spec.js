import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';

test('browser WASM maximum pivoting resolves the retained coupled linear failure with the same ordering', async ({ page }) => {
  const buildPath = 'third_party/klu/BUILD.json', fixturePath = 'tests/fixtures/klu-refined-coupled.json';
  const sourceHashes = numericalSourceHashes(['tests/browser/klu-pivot-control.spec.js']);
  const buildHash = sha256(buildPath), fixtureHash = sha256(fixturePath);
  await page.goto('/');
  const data = await page.evaluate(() => new Promise((resolve, reject) => {
    const code = `self.onmessage=async()=>{try{
      const {solveSparseDirect}=await import('${location.origin}/src/numerics/klu.js');
      const {sparseProduct}=await import('${location.origin}/src/numerics/sparse.js');
      const f=await fetch('${location.origin}/tests/fixtures/klu-refined-coupled.json').then(r=>r.json());
      const a={...f.matrix,rowPtr:Int32Array.from(f.matrix.rowPtr),colIndex:Int32Array.from(f.matrix.colIndex),values:Float64Array.from(f.matrix.values)};
      const b=Float64Array.from(f.rhs);let originalFailure;
      try{solveSparseDirect(a,b,{ordering:'amd',pivotFallback:false});}catch(error){originalFailure=error.code;}
      const r=solveSparseDirect(a,b,{ordering:'amd'}),ax=sparseProduct(a,r.x);
      const expected=Float64Array.from({length:a.n},(_,i)=>Math.sin(.37*i));
      const known=solveSparseDirect(a,sparseProduct(a,expected),{ordering:'amd',pivotTolerance:1});
      self.postMessage({unknowns:a.n,originalFailure,linear:{...r,x:undefined},
        originalRelativeResidual:Math.hypot(...b.map((v,i)=>v-ax[i]))/Math.hypot(...b),
        manufacturedForwardError:Math.max(...known.x.map((v,i)=>Math.abs(v-expected[i])))});
    }catch(error){self.postMessage({error:error.stack});}};`;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = event => { close(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({});
  }));
  expect(data.originalFailure).toBe('KLU_RESIDUAL_LIMIT'); expect(data.linear.ordering).toBe('amd');
  expect(data.linear.pivotTolerance).toBe(1); expect(data.linear.attempts).toHaveLength(2);
  expect(data.linear.relativeResidual).toBeLessThan(1e-10); expect(data.originalRelativeResidual).toBeLessThan(1e-10);
  expect(data.manufacturedForwardError).toBeLessThan(2e-7);
  expect(changedSources(sourceHashes)).toEqual([]); expect(sha256(buildPath)).toBe(buildHash); expect(sha256(fixturePath)).toBe(fixtureHash);
  fs.writeFileSync('docs/current-coupled-pivot-browser.json', JSON.stringify({ date: new Date().toISOString(), physicalAcceptance: false,
    sourceHashes, build: { path: buildPath, sha256: buildHash }, fixture: { path: fixturePath, sha256: fixtureHash },
    scope: 'Actual browser worker, shipped sparse backend with automatic maximum-pivot fallback, fixed coupled Jacobian and manufactured solution; not a nonlinear flow acceptance result.', ...data }) + '\n');
  console.log(JSON.stringify(data));
});
