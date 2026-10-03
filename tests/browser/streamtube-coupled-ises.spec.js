import { test, expect } from '@playwright/test';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directChannelConservation } from '../oracles/streamtube.js';
import fs from 'node:fs';

test('coupled ISES runs all BLs and wakes in a browser worker and retains strict final grid acceptance', async ({ page }) => {
  await page.goto('/');
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const fixture = JSON.parse(fs.readFileSync(new URL('../fixtures/default-coupled-ises-wake-fold.json', import.meta.url)));
  const r = await page.evaluate(({ input, fixture }) => new Promise((resolve, reject) => {
    const coupled = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const ises = new URL('/src/euler/streamtube-coupled-ises.js', location.href).href;
    const source = `self.onmessage=async({data:{input,fixture:f}})=>{try{
      const {solveCoupledStreamtubeIses}=await import(${JSON.stringify(ises)});
      const {createCoupledStreamtubeBody,coupledStreamtubeResult}=await import(${JSON.stringify(coupled)});
      const meshes=[];
      const result=solveCoupledStreamtubeIses(input,{edgeMatching:'section-velocity',maxIterations:12,tolerance:1e-10,
        stepAcceptance:'admissible',onMesh:m=>meshes.push({nodes:m.nodes,iteration:m.iteration.iteration})});
      const system=createCoupledStreamtubeBody(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:Float64Array.from(f.initialBL)});
      const rejected=coupledStreamtubeResult(system,system.initial,{tolerance:1});
      const chart=system.euler.geometryChart(),nodes=rejected.flow.nodes;
      const wakeNormalError=Math.max(...system.euler.layout.bodies.map((b,k)=>{
        const i=b.trailingIndex+1,col=system.euler.layout.nodes[k][i].at(-1).column;
        const n=chart.find(p=>p.column===col).normal;
        const tx=nodes[k][i+1].at(-1).x-nodes[k][i-1].at(-1).x+nodes[k+1][i+1][0].x-nodes[k+1][i-1][0].x;
        const ty=nodes[k][i+1].at(-1).y-nodes[k][i-1].at(-1).y+nodes[k+1][i+1][0].y-nodes[k+1][i-1][0].y;
        return Math.abs((tx*n.x+ty*n.y)/Math.hypot(tx,ty));
      }));
      self.postMessage({result,meshes,endpoint:{wakeNormalError,unknowns:rejected.x.length,converged:rejected.converged,
        quality:rejected.mesh.quality,flowSolved:rejected.mesh.initialization.flowSolved,
        strict:system.admissible(system.initial),sampled:system.admissible(system.initial,{requireConvex:false})}});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({ input, fixture });
  }), { input, fixture });
  expect(r.result.converged).toBe(true); expect(r.result.x.length).toBe(349);
  expect(r.result.boundaryLayer.surfaces).toHaveLength(4); expect(r.result.boundaryLayer.wakes).toHaveLength(2);
  expect(r.result.initialRedistribution.accepted).toBe(true);
  expect(r.result.linearDiagnostics.solves).toBe(r.result.history.length - 1);
  expect(r.meshes.map(m => m.iteration)).toEqual(r.result.history.map(h => h.iteration));
  expect(r.meshes.at(-1).nodes).toEqual(r.result.flow.nodes);
  const flow = r.result.flow;
  for (let g = 0; g < flow.nodes.length; g++) {
    const balance = directChannelConservation({ nodes: flow.nodes[g], sections: flow.sections.map(row => row[g]), cells: flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) expect(Math.max(...balance[key].map(Math.abs))).toBeLessThan(2e-9);
  }
  expect(r.endpoint.wakeNormalError).toBeLessThan(1e-12); expect(r.endpoint.unknowns).toBe(8195); expect(r.endpoint.sampled).toBe(true); expect(r.endpoint.strict).toBe(false);
  expect(r.endpoint.quality.valid).toBe(false); expect(r.endpoint.converged).toBe(false); expect(r.endpoint.flowSolved).toBe(false);
  console.log(JSON.stringify({ families: r.result.families, updates: r.result.history.length - 1, endpoint: r.endpoint }));
});
