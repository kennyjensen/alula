import { test, expect } from '@playwright/test';

test('Drela elliptic grid SLOR executes and reports intermediate geometry in a browser module worker', async ({ page }) => {
  await page.goto('/');
  const nx = 16, nt = 8, eta = Array.from({ length: nt + 1 }, (_, j) => Math.expm1(.8 * j / nt) / Math.expm1(.8));
  const exact = Array.from({ length: nx + 1 }, (_, i) => eta.map(e => {
    const r = 2 * Math.exp(-.5 * e), theta = -.4 + .8 * i / nx;
    return { x: r * Math.cos(theta), y: r * Math.sin(theta) };
  }));
  const nodes = exact.map((row, i) => row.map((p, j) => !i || i === nx || !j || j === nt ? { ...p }
    : { x: (1 - eta[j]) * row[0].x + eta[j] * row[nt].x, y: (1 - eta[j]) * row[0].y + eta[j] * row[nt].y }));
  const result = await page.evaluate(input => new Promise((resolve, reject) => {
    const module = new URL('/src/geometry/elliptic-streamtube-grid.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createEllipticStreamtubeGrid,smoothEllipticStreamtubeGrid}=await import(${JSON.stringify(module)});
      const r=smoothEllipticStreamtubeGrid(createEllipticStreamtubeGrid(data),{
        tolerance:1e-10,onSweep:(h,nodes)=>{if(h.iteration%10===0)self.postMessage({kind:'grid',history:h,nodes});}
      });self.postMessage({kind:'result',result:r});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const previews = [], close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => {
      if (data.error) { close(); reject(new Error(data.error)); }
      else if (data.kind === 'grid') previews.push(data);
      else { close(); resolve({ ...data.result, previews }); }
    };
    worker.postMessage(input);
  }), { nodes, massFlows: eta.slice(1).map((v, j) => v - eta[j]) });
  expect(result.converged).toBe(true); expect(result.quality.valid).toBe(true);
  expect(result.history.at(-1).residual).toBeLessThan(1e-10);
  expect(result.previews.length).toBeGreaterThan(1); expect(result.previews[0].history.iteration).toBe(0);
  expect(result.previews[0].nodes).toEqual(nodes);
  let error = 0;
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    error = Math.max(error, Math.abs((Math.atan2(p.y, p.x) + .4) / .8 - i / nx),
      Math.abs(-2 * Math.log(Math.hypot(p.x, p.y) / 2) - eta[j]));
    if (!i || i === nx || !j || j === nt) expect(p).toEqual(nodes[i][j]);
  }));
  expect(error).toBeLessThan(.00021);
});
