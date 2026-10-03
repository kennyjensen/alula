import { test, expect } from '@playwright/test';
import fs from 'node:fs';

test('browser inverse BL initialization escapes a clipped closure while matching native residuals', async ({ page }) => {
  const fixture = JSON.parse(fs.readFileSync(new URL('../fixtures/streamtube-laminar-inverse.json', import.meta.url)));
  await page.goto('/');
  const result = await page.evaluate(data => new Promise((resolve, reject) => {
    const source = `self.onmessage=async({data})=>{try{
      const {createIntegralKernel}=await import('/src/viscous/integral.js');
      const {initializeStreamtubeBLStation}=await import('/src/euler/streamtube-boundary-layers.js');
      const kernel=createIntegralKernel({...data.parameters,exactJacobian:true});
      const before=kernel.station(data.downstream);
      const r=initializeStreamtubeBLStation({upstream:data.upstream,s:data.downstream.s,ue:data.downstream.ue,
        regime:'laminar',reynolds:data.parameters.reynolds,interval:kernel.interval,properties:kernel.station});
      const block=kernel.interval({upstream:data.upstream,downstream:r.state,regime:'laminar'});
      self.postMessage({before,result:r,block});
    }catch(error){self.postMessage({error:error.message});}};`;
    // Blob workers need absolute module URLs, since their own base is blob:.
    const code = source.replaceAll("import('/src/", `import('${location.origin}/src/`);
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(data);
  }), fixture);
  expect(result.before.hk).toBe(1.05); expect(result.before.rawHk).toBeGreaterThan(1);
  expect(result.before.rawHk).toBeLessThan(1.000001);
  expect(result.result.converged).toBe(true); expect(result.result.mode).toBe('inverse');
  expect(Math.abs(result.block.properties.rawHk - result.result.targetHK)).toBeLessThan(1e-10);
  const native = fixture.inverse.native.cases.find(c => c.name === 'existing inverse initializer');
  for (const key of ['aux', 'theta', 'deltaStar', 'ue']) expect(Math.abs(result.result.state[key] - native.downstream[key])).toBeLessThan(1e-10);
  result.block.residual.forEach((v, i) => expect(Math.abs(v - native.nativeResidual[i])).toBeLessThan(1e-10));
  console.log(JSON.stringify({ mode: result.result.mode, updates: result.result.history.length - 1,
    state: result.result.state, residual: result.block.residual }));
});
