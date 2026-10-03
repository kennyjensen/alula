import { test, expect } from '@playwright/test';
import fs from 'node:fs';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';
import { directChannelConservation } from '../oracles/streamtube.js';
import { numericalSourceHashes, changedSources } from '../../scripts/validation/provenance.js';

test('a browser worker redistributes surface stations and reconverges all four BLs and both wakes', async ({ page }) => {
  const sourceHashes = numericalSourceHashes(['tests/browser/streamtube-coupled-redistribution.spec.js',
    'scripts/validation/coupled-station-control.js', 'tests/fixtures/intrinsic-body.js']);
  await page.goto('/');
  const data = await page.evaluate(() => new Promise((resolve, reject) => {
    const code = `self.onmessage=async()=>{try{
        const {coupledStationControl}=await import('${location.origin}/scripts/validation/coupled-station-control.js');
        const {source,mapped,seed,result}=coupledStationControl(),before=source.evaluate(source.initial);
        self.postMessage({result,unknowns:source.n,diagnostics:mapped.diagnostics,
          originalAllocation:before.outer.allocation,transferredAllocation:seed.outer.allocation,
          originalTrips:source.bl.surfaces.map(b=>b.tripParameter),transferredTrips:mapped.system.bl.surfaces.map(b=>b.tripParameter),
          originalWakeStates:source.bl.wakes.map(w=>w.ids.map(id=>before.layers.states[id])),
          transferredWakeStates:mapped.system.bl.wakes.map(w=>w.ids.map(id=>seed.layers.states[id]))});
      }catch(error){self.postMessage({error:error.stack});}};`;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = event => { close(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({});
  }));
  const r = data.result;
  expect(r.converged).toBe(true); expect(r.x.length).toBe(data.unknowns);
  expect(r.boundaryLayer.surfaces).toHaveLength(4); expect(r.boundaryLayer.wakes).toHaveLength(2);
  expect(Math.max(...Object.values(r.families))).toBeLessThan(1e-10);
  expect(data.transferredAllocation).toEqual(data.originalAllocation);
  expect(data.transferredTrips).toEqual(data.originalTrips);
  expect(data.transferredWakeStates).toEqual(data.originalWakeStates);
  expect(data.diagnostics.maxNodeMovement).toBeGreaterThan(0);
  const geometry = directStreamtubeVolumeGeometry(r.flow.nodes);
  expect(geometry.valid).toBe(true); expect(geometry.concavePrimal).toHaveLength(0);
  const conservation = r.flow.nodes.map((nodes, g) => directChannelConservation({ nodes,
    sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) }));
  for (const c of conservation) for (const key of ['maxLocal', 'total', 'internalCancellation'])
    for (const k of key === 'internalCancellation' ? [0, 1, 2, 3] : [0, 3]) expect(Math.abs(c[key][k])).toBeLessThan(2e-9);
  expect(changedSources(sourceHashes)).toEqual([]);
  fs.writeFileSync('docs/current-coupled-surface-spacing-browser.json', JSON.stringify({ date: new Date().toISOString(),
    physicalAcceptance: false, sourceHashes, scope: 'Actual browser worker, complete controlled two-element solve; no airfoil-accuracy claim.',
    geometry, conservation, ...data }, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v) + '\n');
  console.log(JSON.stringify({ unknowns: data.unknowns, families: r.families, iterations: r.history.length - 1 }));
});
