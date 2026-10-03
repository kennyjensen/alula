import { test, expect } from '@playwright/test';

test('research variational grid and physical boundary refinement run in a module worker and approach an independent harmonic map', async ({ page }) => {
  await page.goto('/');
  const n = 8, coordinates = Array.from({ length: n + 1 }, (_, i) => Array.from({ length: n + 1 }, (_, j) => ({ x: -1 + 2 * i / n, y: j / n })));
  const nodes = coordinates.map(row => row.map(p => ({ x: p.x,
    y: 2 * (p.y - .12 * p.x * p.x) / (1 + Math.sqrt(1 - .48 * (p.y - .12 * p.x * p.x))) })));
  const result = await page.evaluate(input => new Promise((resolve, reject) => {
    const module = new URL('/src/geometry/tests/variational-streamtube-grid.js', location.href).href;
    const refiner = new URL('/src/geometry/tests/refine-potential-grid.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createVariationalStreamtubeGrid,smoothVariationalStreamtubeGrid}=await import(${JSON.stringify(module)});
      const {refinePotentialGrid}=await import(${JSON.stringify(refiner)});
      const coarse=smoothVariationalStreamtubeGrid(createVariationalStreamtubeGrid(data));
      if(!coarse.converged)throw new Error(coarse.reason);
      const refined=refinePotentialGrid({nodes:coarse.nodes,coordinates:data.coordinates},{boundaryAt:(i,j)=>{
        const u=-1+2*i/8,v=j/8;
        return {point:{x:u,y:2*(v-.12*u*u)/(1+Math.sqrt(1-.48*(v-.12*u*u)))},coordinate:{x:u,y:v}};
      }});
      const result=smoothVariationalStreamtubeGrid(createVariationalStreamtubeGrid(refined),{
        onIteration:history=>self.postMessage({kind:'progress',history})
      });self.postMessage({kind:'result',result:{...result,coarse:coarse.nodes,refinedSeed:refined.nodes,coordinates:refined.coordinates}});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const progress = [], close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => {
      if (data.error) { close(); reject(new Error(data.error)); }
      else if (data.kind === 'progress') progress.push(data.history);
      else { close(); resolve({ ...data.result, progress }); }
    };
    worker.postMessage(input);
  }), { nodes, coordinates });
  expect(result.converged).toBe(true); expect(result.quality.valid).toBe(true);
  expect(result.physicsValidated).toBe(false); expect(result.progress.length).toBeGreaterThan(1);
  let error = 0;
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    const w = result.coordinates[i][j];
    error = Math.max(error, Math.hypot(p.x - w.x, p.y + .12 * (p.x * p.x - p.y * p.y) - w.y));
    if (!i || i === 2 * n || !j || j === 2 * n) expect(p).toEqual(result.refinedSeed[i][j]);
    if (i % 2 === 0 && j % 2 === 0) expect(result.refinedSeed[i][j]).toEqual(result.coarse[i / 2][j / 2]);
  }));
  expect(error).toBeLessThan(3.4e-6); expect(error).toBeGreaterThan(3.2e-6);
});
