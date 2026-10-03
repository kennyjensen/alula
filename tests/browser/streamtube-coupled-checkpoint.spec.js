import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('a serialized coupled checkpoint survives worker replacement and resumes every BL and wake', async ({ page }) => {
  const sourceHashes = numericalSourceHashes(['tests/browser/streamtube-coupled-checkpoint.spec.js', 'tests/fixtures/intrinsic-body.js']);
  await page.goto('/');
  const data = await page.evaluate(async () => {
    const code = `self.onmessage=async({data})=>{try{
      const {solveCoupledStreamtubeIses}=await import('${location.origin}/src/euler/streamtube-coupled-ises.js');
      const {intrinsicBodyFixture}=await import('${location.origin}/tests/fixtures/intrinsic-body.js');
      const controls={edgeMatching:'section-velocity',maxIterations:12,tolerance:1e-10,stepAcceptance:'admissible'};
      if(data.mode==='stop'){
        let checkpoint;
        try{solveCoupledStreamtubeIses(intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2}),{...controls,
          onCheckpoint:(c,d)=>{checkpoint=c;if(d.history.at(-1).iteration===3)throw new Error('stop after checkpoint');}});}
        catch(error){if(error.message!=='stop after checkpoint')throw error;}
        self.postMessage({checkpoint});
      }else{
        const r=data.mode==='resume'?solveCoupledStreamtubeIses(undefined,{...controls,resume:data.checkpoint})
          :solveCoupledStreamtubeIses(intrinsicBodyFixture({elements:2,bodySegments:4,tubes:2}),controls);
        self.postMessage({converged:r.converged,x:r.x,nodes:r.flow.nodes,families:r.families,
          surfaces:r.boundaryLayer.surfaces.length,wakes:r.boundaryLayer.wakes.length,
          iterations:r.history.length-1,redistribution:r.initialRedistribution});
      }
    }catch(error){self.postMessage({error:error.stack});}};`;
    const run = data => new Promise((resolve, reject) => {
      const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), w = new Worker(url, { type: 'module' });
      const close = () => { w.terminate(); URL.revokeObjectURL(url); };
      w.onerror = e => { close(); reject(new Error(e.message)); };
      w.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
      w.postMessage(data);
    });
    const first = await run({ mode: 'stop' });
    const checkpoint = JSON.parse(JSON.stringify(first.checkpoint, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
    const resumed = await run({ mode: 'resume', checkpoint }), full = await run({ mode: 'full' });
    let maximumStateDifference = 0, maximumNodeDifference = 0;
    resumed.x.forEach((v, i) => maximumStateDifference = Math.max(maximumStateDifference, Math.abs(v - full.x[i])));
    resumed.nodes.forEach((g, gi) => g.forEach((row, i) => row.forEach((p, j) => {
      const q = full.nodes[gi][i][j]; maximumNodeDifference = Math.max(maximumNodeDifference, Math.hypot(p.x - q.x, p.y - q.y));
    })));
    return { checkpoint, resumed: { ...resumed, x: undefined, nodes: undefined }, full: { ...full, x: undefined, nodes: undefined },
      maximumStateDifference, maximumNodeDifference };
  });
  expect(data.resumed.converged).toBe(true); expect(data.full.converged).toBe(true);
  expect(data.resumed.surfaces).toBe(4); expect(data.resumed.wakes).toBe(2);
  expect(data.resumed.redistribution.resumed).toBe(true); expect(data.resumed.redistribution.passages).toEqual([]);
  expect(data.full.iterations - data.resumed.iterations).toBe(3);
  expect(data.maximumStateDifference).toBeLessThan(1e-12); expect(data.maximumNodeDifference).toBeLessThan(1e-12);
  expect(changedSources(sourceHashes)).toEqual([]);
  fs.writeFileSync('docs/current-coupled-checkpoint-browser.json', JSON.stringify({ date: new Date().toISOString(), physicalAcceptance: false,
    sourceHashes, scope: 'Complete small two-element worker interruption/serialization/replacement/recovery, compared with an uninterrupted solve.', ...data }) + '\n');
  console.log(JSON.stringify({ maximumStateDifference: data.maximumStateDifference, maximumNodeDifference: data.maximumNodeDifference,
    resumed: data.resumed, fullIterations: data.full.iterations }));
});
