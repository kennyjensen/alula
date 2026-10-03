import { test, expect } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { sparseProduct } from '../../src/numerics/sparse.js';

const fixture = JSON.parse(readFileSync(new URL('../fixtures/klu-refined-coupled.json', import.meta.url)));

test('browser WASM retries an inaccurate factor ordering and certifies the original coupled matrix', async ({ page }) => {
  await page.goto('/');
  const r = await page.evaluate(fixture => new Promise((resolve, reject) => {
    const module = new URL('/src/numerics/klu.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {solveSparseDirect}=await import(${JSON.stringify(module)});
      self.postMessage(solveSparseDirect(data.matrix,Float64Array.from(data.rhs)));
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(fixture);
  }), fixture);
  expect(r.backend).toBe('klu-wasm'); expect(r.ordering).toBe('colamd');
  expect(r.attempts).toHaveLength(2); expect(r.attempts[0].relativeResidual).toBeGreaterThan(1e-4);
  expect(r.relativeResidual).toBeLessThan(1e-10);
  const ax = sparseProduct(fixture.matrix, r.x), b = fixture.rhs;
  expect(Math.hypot(...b.map((v, i) => v - ax[i])) / Math.hypot(...b)).toBeLessThan(1e-10);
  console.log(JSON.stringify({ unknowns: fixture.matrix.n, attempts: r.attempts }));
});
