import { test, expect } from '@playwright/test';
import fs from 'node:fs';

test('browser worker constructs an admissible nested coupled grid while conserving parent streamtube masses', async ({ page }) => {
  const fixture = JSON.parse(fs.readFileSync(new URL('../fixtures/streamtube-isentropic-coupled-root.json', import.meta.url)));
  await page.goto('/');
  const result = await page.evaluate(fixture => new Promise((resolve, reject) => {
    const source = `self.onmessage=async({data})=>{try{
      const {createCoupledStreamtubeBody}=await import('/src/euler/streamtube-coupled.js');
      const {refineCoupledStreamtubeBody}=await import('/src/euler/streamtube-coupled-refinement.js');
      const parent=createCoupledStreamtubeBody(data.input,{...data.options,initialEuler:data.initialEuler,initialBL:data.initialBL});
      const before=parent.evaluate(parent.initial), r=refineCoupledStreamtubeBody(data.input,parent), after=r.system.evaluate(r.system.initial);
      let massError=0,nodeError=0;
      before.outer.allocation.groups.forEach((group,g)=>group.forEach((cell,j)=>{
        const a=after.outer.allocation.groups[g][2*j].massFlow,b=after.outer.allocation.groups[g][2*j+1].massFlow;
        massError=Math.max(massError,Math.abs((a+b)/cell.massFlow-1));
      }));
      before.outer.nodes.forEach((grid,g)=>grid.forEach((row,i)=>row.forEach((p,j)=>{
        const q=after.outer.nodes[g][2*i][2*j];nodeError=Math.max(nodeError,Math.hypot(q.x-p.x,q.y-p.y));
      })));
      self.postMessage({diagnostics:r.diagnostics,massError,nodeError,admissible:r.system.admissible(r.system.initial),families:after.families});
    }catch(error){self.postMessage({error:error.message});}};`;
    const code = source.replaceAll("import('/src/", `import('${location.origin}/src/`);
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(fixture);
  }), fixture);
  expect(result.admissible).toBe(true); expect(result.diagnostics.unknowns).toBe(11929);
  expect(result.diagnostics.quality.valid).toBe(true); expect(result.massError).toBeLessThan(1e-14);
  expect(result.nodeError).toBeLessThan(2e-12); expect(result.families.euler).toBeGreaterThan(1e-4);
  console.log(JSON.stringify(result));
});
