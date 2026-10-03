import { test, expect } from '@playwright/test';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';
import { directChannelConservation } from '../oracles/streamtube.js';
import { directStreamtubeVolumeGeometry } from '../oracles/streamtube-control-volume-geometry.js';
import { directBodyConservation } from '../oracles/streamtube-body.js';
import { readFileSync } from 'node:fs';

test('ISES research Newton and grid-maintenance sequence runs in a browser WASM worker with independent conservation', async ({ page }) => {
  await page.goto('/');
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const { result: r, sampled, meshes } = await page.evaluate(input => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-ises-update.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {solveStreamtubeIses}=await import(${JSON.stringify(module)}),meshes=[];
      const result=solveStreamtubeIses(data,{tolerance:1e-11,onMesh:m=>meshes.push({iteration:m.iteration,nodes:m.nodes})});
      const sampled=solveStreamtubeIses(data,{tolerance:1e-11,iterationGeometry:'ises-sampled'});
      self.postMessage({result,sampled,meshes});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(r.converged).toBe(true); expect(r.surfaces).toHaveLength(4); expect(r.lastRejectedStep).toBeNull();
  expect(r.finalQuality.valid).toBe(true); expect(r.residualConverged).toBe(true);
  expect(sampled.converged).toBe(true); expect(sampled.finalQuality.valid).toBe(true);
  expect(sampled.nodes).toEqual(r.nodes); expect(sampled.x).toEqual(r.x);
  expect(r.diagnostics.residual).toBeLessThan(1e-11); expect(r.initialRedistribution.accepted).toBe(true);
  expect(r.initialRedistribution.passages.map(p => p.fixedBanks)).toEqual([[false, true], [true, true], [true, false]]);
  expect(r.initialRedistribution.passages.every(p => p.pairs === 5)).toBe(true);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  expect(meshes.map(m => m.iteration.iteration)).toEqual(r.history.map(h => h.iteration));
  expect(meshes.at(-1).nodes).toEqual(r.nodes);
  for (let g = 0; g < r.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.nodes[g], sections: r.sections.map(row => row[g]), cells: r.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
  console.log(JSON.stringify({ unknowns: r.x.length, iterations: r.history.length - 1, residual: r.diagnostics.residual, initialRedistribution: r.initialRedistribution }));
});

test('the default wall-refined isentropic Euler case reaches a positive root in a browser worker', async ({ page }) => {
  test.setTimeout(75000);
  const cold = JSON.parse(readFileSync(new URL('../../docs/default-euler-density-newton.json', import.meta.url)));
  const seed = cold.cases.find(c => c.streamwiseMode === 'momentum' && c.stepMethod === 'density-newton');
  const input = { ...seed.restart.input, streamwiseMode: 'isentropic', normalStencil: 'body-stations', stagnationMotion: 'walls-only', geometryDomain: 'positive-simple' };
  await page.goto('/'); await expect(page.locator('#status')).toContainText('Solved');
  const data = await page.evaluate(({ input, initialEuler }) => new Promise((resolve, reject) => {
    const base = new URL('/src/euler/', location.href).href;
    const code = `self.onmessage=async({data})=>{try{
      const {createStreamtubeBodySystem}=await import(${JSON.stringify(base + 'streamtube-body.js')});
      const {refineStreamtubeBody}=await import(${JSON.stringify(base + 'streamtube-refinement.js')});
      const {solveStreamtubeIses}=await import(${JSON.stringify(base + 'streamtube-ises-update.js')});
      const system=createStreamtubeBodySystem(data.input);
      const initial=system.adoptGeometry(Float64Array.from(data.initialEuler.x),data.initialEuler.nodes);
      const counts=system.layout.tubes.map((n,g)=>Array.from({length:n},(_,j)=>g>0&&j===0||g<system.layout.elements&&j===n-1?2:1));
      const refined=refineStreamtubeBody(data.input,system,{initial,streamwiseFactor:1,normalSubdivisions:counts});
      const r=solveStreamtubeIses(refined.input,{initialEuler:refined.initialEuler,maxIterations:16,
        iterationGeometry:'ises-sampled',stepAcceptance:'admissible'});
      self.postMessage({converged:r.converged,reason:r.reason,diagnostics:r.diagnostics,quality:r.finalQuality,
        history:r.history,conditions:refined.system.conditions,bodies:refined.input.bodies,
        flow:{nodes:r.nodes,sections:r.sections,cells:r.cells,diagnosticForces:r.diagnosticForces},unknowns:r.x.length});
    }catch(e){self.postMessage({error:e.message});}};`;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage({ input, initialEuler });
  }), { input, initialEuler: seed.restart.initialEuler });
  expect(data.converged, data.reason).toBe(true); expect(data.unknowns).toBe(7015);
  expect(data.diagnostics.residual).toBeLessThanOrEqual(1e-10);
  expect(data.quality.valid).toBe(true); expect(data.quality.minCornerSine).toBeGreaterThan(.04);
  expect(directStreamtubeVolumeGeometry(data.flow.nodes).valid).toBe(true);
  const conservation = directBodyConservation(data.flow, data.bodies, data.conditions);
  expect(Math.abs(conservation.balance[0])).toBeLessThan(2e-12);
  expect(Math.abs(conservation.balance[3])).toBeLessThan(2e-11);
  expect(Math.max(...conservation.cutTraction.map(Math.abs))).toBeLessThan(2e-9);
  expect(data.diagnostics.maxStagnationPressureError).toBeLessThan(1e-9);
  // Isentropy replaces streamwise momentum. Report its finite defect; this
  // browser equation/geometry check does not certify physical force accuracy.
  console.log(JSON.stringify({ unknowns: data.unknowns, iterations: data.history.length - 1,
    residual: data.diagnostics.residual, quality: data.quality, conservation: conservation.balance }));
});

test('refined Euler common-step retry runs in a browser worker and publishes only valid volumes', async ({ page }) => {
  await page.goto('/');
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'walls-only', normalStencil: 'body-stations' };
  const { result, meshes, parentMasses } = await page.evaluate(input => new Promise((resolve, reject) => {
    const base = new URL('/src/euler/', location.href).href;
    const code = `self.onmessage=async({data})=>{try{
      const {createStreamtubeBodySystem,solveStreamtubeBody}=await import(${JSON.stringify(base + 'streamtube-body.js')});
      const {refineStreamtubeBody}=await import(${JSON.stringify(base + 'streamtube-refinement.js')});
      const {solveStreamtubeIses}=await import(${JSON.stringify(base + 'streamtube-ises-update.js')});
      const system=createStreamtubeBodySystem(data),parent=solveStreamtubeBody(system,{maxIterations:8});
      if(!parent.converged)throw new Error('Parent equations did not converge');
      const refined=refineStreamtubeBody(data,system,{initial:parent.x}),meshes=[];
      const result=solveStreamtubeIses(refined.input,{initialEuler:refined.initialEuler,maxIterations:1,
        stepAcceptance:'admissible',onMesh:m=>meshes.push(m.nodes)});
      self.postMessage({result,meshes,parentMasses:parent.allocation.groups.map(g=>g.map(t=>t.massFlow))});
    }catch(e){self.postMessage({error:e.message});}};`;
    const url = URL.createObjectURL(new Blob([code], { type: 'text/javascript' })), worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(result.converged).toBe(false); expect(result.history).toHaveLength(2);
  expect(result.history[1].backtracks).toBeGreaterThan(0);
  expect(result.linearDiagnostics.solves).toBe(1); expect(meshes).toHaveLength(2);
  for (const nodes of meshes) expect(directStreamtubeVolumeGeometry(nodes).valid).toBe(true);
  // Capture can move in the Newton update. Child fractions within each
  // passage must still follow the parent subdivision proportions.
  result.allocation.groups.forEach((group, g) => {
    const sum = group.reduce((s, t) => s + t.massFlow, 0), parentSum = parentMasses[g].reduce((a, b) => a + b, 0);
    for (let j = 0; j < parentMasses[g].length; j++)
      expect(Math.abs((group[2 * j].massFlow + group[2 * j + 1].massFlow) / sum - parentMasses[g][j] / parentSum)).toBeLessThan(1e-13);
  });
});
