import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';
import { directChannelConservation } from '../oracles/streamtube.js';

test('ISES maintenance recovers the actual refined two-element grid with every BL and wake in a browser worker', async ({ page }) => {
  test.setTimeout(120000);
  const source = 'docs/current-coupled-checkpoint-default-failed-import.json';
  const referencePath = 'docs/current-coupled-checkpoint-default-root-import.json';
  const saved = JSON.parse(fs.readFileSync(source)), reference = JSON.parse(fs.readFileSync(referencePath));
  expect(changedSources(saved.sourceHashes)).toEqual([]); expect(changedSources(reference.sourceHashes)).toEqual([]);
  expect(saved.result.converged).toBe(false); expect(reference.result.converged).toBe(true);
  const sourceHashes = numericalSourceHashes(['tests/browser/streamtube-streamwise-recovery.spec.js',
    'tests/oracles/streamtube-control-volume-geometry.js', 'tests/oracles/streamtube.js']);
  const sourceHash = sha256(source), referenceHash = sha256(referencePath);
  await page.goto('/');
  const data = await page.evaluate(f => new Promise((resolve, reject) => {
    const code = `self.onmessage=async({data:f})=>{try{
      const {solveCoupledStreamtubeIses}=await import('${location.origin}/src/euler/streamtube-coupled-ises.js');
      const updates=[],r=solveCoupledStreamtubeIses(f.input,{...f.options,initialEuler:f.initialEuler,initialBL:f.initialBL,
        maxIterations:8,tolerance:1e-10,iterationGeometry:'ises-sampled',stepAcceptance:'admissible',
        onMesh:m=>updates.push({iteration:m.iteration.iteration,families:m.coupledFamilies})});
      self.postMessage({result:r,updates});
    }catch(error){self.postMessage({error:error.stack});}};`;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = event => { close(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(f);
  }), saved.restart);
  const r = data.result;
  expect(r.converged).toBe(true); expect(r.x.length).toBe(19433);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(r.history).toHaveLength(7); expect(r.history.slice(1).every(h => h.step === 1 && h.backtracks === 0)).toBe(true);
  expect(data.updates.map(m => m.iteration)).toEqual(r.history.map(h => h.iteration));
  expect(r.initialRedistribution.accepted).toBe(true);
  expect(r.initialRedistribution.passages.every(p => p.pairs === 5)).toBe(true);
  const geometry = directStreamtubeVolumeGeometry(r.flow.nodes);
  expect(geometry.valid).toBe(true); expect(geometry.concavePrimal).toHaveLength(0);
  const conservation = r.flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) }));
  for (const c of conservation) for (const key of ['maxLocal', 'total', 'internalCancellation'])
    for (const k of key === 'internalCancellation' ? [0, 1, 2, 3] : [0, 3]) expect(Math.abs(c[key][k])).toBeLessThan(2e-9);
  let maximumNodeDifference = 0, maximumStateDifference = 0;
  r.flow.nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    const q = reference.restart.initialEuler.nodes[g][i][j];
    maximumNodeDifference = Math.max(maximumNodeDifference, Math.hypot(p.x - q.x, p.y - q.y));
  })));
  const expected = [...reference.restart.initialEuler.x, ...reference.restart.initialBL];
  r.x.forEach((v, i) => { maximumStateDifference = Math.max(maximumStateDifference, Math.abs(v - expected[i])); });
  expect(maximumNodeDifference).toBeLessThan(1e-10); expect(maximumStateDifference).toBeLessThan(1e-9);
  expect(changedSources(sourceHashes)).toEqual([]); expect(sha256(source)).toBe(sourceHash); expect(sha256(referencePath)).toBe(referenceHash);
  const ne = r.x.length - 4 * r.boundaryLayer.stations.length;
  fs.writeFileSync('docs/current-coupled-checkpoint-default-browser.json', JSON.stringify({ date: new Date().toISOString(), physicalAcceptance: false,
    source: { path: source, sha256: sourceHash }, reference: { path: referencePath, sha256: referenceHash }, sourceHashes,
    scope: 'Actual browser worker recovery of the complete default streamwise-refined coupled state, compared with the audited Node root; not a new cold GUI startup or independent physical reference.',
    result: { converged: r.converged, families: r.families, unknowns: r.x.length, history: r.history, quality: r.mesh.quality,
      surfaces: r.boundaryLayer.surfaces, wakes: r.boundaryLayer.wakes, initialRedistribution: r.initialRedistribution },
    restart: { input: r.solverInput, options: r.coupledOptions, initialEuler: { x: r.x.slice(0, ne), nodes: r.flow.nodes,
      undisplacedNodes: r.flow.undisplacedNodes }, initialBL: r.x.slice(ne) },
    geometry, conservation, updates: data.updates, maximumNodeDifference, maximumStateDifference },
  (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v) + '\n');
  console.log(JSON.stringify({ unknowns: r.x.length, families: r.families, updates: r.history.length - 1,
    maximumNodeDifference, maximumStateDifference }));
});
