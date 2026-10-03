import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';

const cases = [
  ...[1, 2].map(streamwiseFactor => ({ streamwiseFactor, normalInterpolation: 'linear' })),
  { streamwiseFactor: 2, normalInterpolation: 'streamfunction-quadratic' }
];
for (const controls of cases) test(`default independent-bank refinement retains the complete fluid boundary in a worker, streamwise factor ${controls.streamwiseFactor}, ${controls.normalInterpolation}`, async ({ page }) => {
  await page.goto('/');
  const saved = JSON.parse(fs.readFileSync(new URL('../../docs/default-wake-initialization.json', import.meta.url)));
  const fixture = saved.cases.find(c => c.wakeGeometry === 'independent-banks').restart;
  const r = await page.evaluate(({ fixture, controls }) => new Promise((resolve, reject) => {
    const coupled = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const refine = new URL('/src/euler/streamtube-coupled-refinement.js', location.href).href;
    const code = `self.onmessage=async({data:{fixture:f,controls:{streamwiseFactor,normalInterpolation}}})=>{try{
      const {createCoupledStreamtubeBody}=await import(${JSON.stringify(coupled)});
      const {refineCoupledStreamtubeBody}=await import(${JSON.stringify(refine)});
      const source=createCoupledStreamtubeBody(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:Float64Array.from(f.initialBL)});
      const before=source.evaluate(source.initial),counts=source.euler.layout.tubes.map(n=>Array(n).fill(1));
      counts[0][counts[0].length-1]=normalInterpolation==='linear'?2:4;
      if(normalInterpolation!=='linear')counts[1][0]=4;
      const child=refineCoupledStreamtubeBody(f.input,source,{streamwiseFactor,normalInterpolation,normalSubdivisions:counts});
      const v=child.system.evaluate(child.system.initial);
      let nodeError=0,massError=0;
      before.outer.nodes.forEach((grid,g)=>{const positions=[0];for(const n of counts[g])positions.push(positions.at(-1)+n);
        grid.forEach((row,i)=>row.forEach((p,j)=>{const q=v.outer.nodes[g][streamwiseFactor*i][positions[j]];
          nodeError=Math.max(nodeError,Math.hypot(p.x-q.x,p.y-q.y));}));
        before.outer.allocation.groups[g].forEach((p,j)=>{const total=v.outer.allocation.groups[g].slice(positions[j],positions[j+1]).reduce((s,q)=>s+q.massFlow,0);
          massError=Math.max(massError,Math.abs(total/p.massFlow-1));});
      });
      self.postMessage({nodeError,massError,unknowns:child.system.n,nodes:v.outer.nodes,diagnostics:child.diagnostics,
        families:v.families,surfaces:child.system.bl.surfaces.length,wakes:child.system.bl.wakes.length});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({ fixture, controls });
  }), { fixture, controls });
  expect(r.nodeError).toBeLessThan(2e-12); expect(r.massError).toBeLessThan(1e-14);
  expect(r.surfaces).toBe(4); expect(r.wakes).toBe(2); expect(r.diagnostics.quality.valid).toBe(true);
  expect(directStreamtubeVolumeGeometry(r.nodes).valid).toBe(true);
  expect(Math.max(...Object.values(r.families))).toBeGreaterThan(1e-4);
  expect(r.diagnostics.normalInterpolation).toBe(controls.normalInterpolation);
  if (controls.normalInterpolation === 'streamfunction-quadratic') {
    expect(r.unknowns).toBe(19859);
    expect(r.diagnostics.massInterpolation.some(d => d.maximumFractionChange > .1)).toBe(true);
  }
  console.log(JSON.stringify({ ...controls, unknowns: r.unknowns, nodeError: r.nodeError, massError: r.massError, families: r.families }));
});
