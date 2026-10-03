import { test, expect } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { directChannelConservation } from '../oracles/streamtube.js';
import { intrinsicBodyFixture } from '../fixtures/intrinsic-body.js';

for (const [name, file, unknowns, stations] of [
  ['coarse', 'streamtube-matched-coupled-root.json', 3185, 79],
  ['refined', 'streamtube-matched-refined-root.json', 7919, 137],
  ['isentropic coarse', 'streamtube-isentropic-coupled-root.json', 3185, 79],
  ['isentropic refined', 'streamtube-isentropic-refined-root.json', 7919, 137],
  ['isentropic nested', 'streamtube-isentropic-nested-root.json', 11929, 157],
  ['isentropic stagnation-scaled', 'streamtube-isentropic-stagnation-root.json', 23097, 157],
  ['isentropic section-speed', 'streamtube-section-velocity-root.json', 3185, 79],
]) test(`browser WASM worker restores the ${name} attached coupled root after a BL perturbation`, async ({ page }) => {
  const fixture = JSON.parse(readFileSync(new URL(`../fixtures/${file}`, import.meta.url)));
  await page.goto('/');
  const r = await page.evaluate(fixture => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createCoupledStreamtubeBody,solveCoupledStreamtubeBody}=await import(${JSON.stringify(module)});
      const system=createCoupledStreamtubeBody(data.input,{...data.options,
        initialEuler:data.initialEuler,initialBL:data.initialBL});
      const initial=system.initial.map((v,i)=>i<system.ne?v:v*(1+1e-4*Math.sin(i)));
      self.postMessage(solveCoupledStreamtubeBody(system,{initial,maxIterations:4,tolerance:1e-10,
        stepMethod:data.input.streamwiseMode==='isentropic'?'dogleg':'newton'}));
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(fixture);
  }), fixture);
  expect(r.converged).toBe(true); expect(r.x).toHaveLength(unknowns);
  expect(r.history[0].residual).toBeGreaterThan(1e-7);
  expect(r.linearDiagnostics.solves).toBeGreaterThan(0); expect(r.linearDiagnostics.solves).toBeLessThanOrEqual(4);
  expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  for (const error of Object.values(r.families)) expect(error).toBeLessThan(1e-10);
  expect(r.boundaryLayer.stations).toHaveLength(stations);
  if (name === 'refined') expect(r.linearDiagnostics.orderingFallbacks).toBe(1);
  expect(r.boundaryLayer.surfaces).toHaveLength(2); expect(r.boundaryLayer.wakes).toHaveLength(1);
  expect(r.mesh.quality.valid).toBe(true); expect(r.mesh.cells.every(c => c.length === 4)).toBe(true);
  expect(r.cl).toBeNull(); expect(r.cd).toBeNull(); expect(r.cm).toBeNull();
  for (const [g, nodes] of r.flow.nodes.entries()) {
    const c = directChannelConservation({ nodes, sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) });
    if (fixture.input.streamwiseMode === 'isentropic') {
      for (const key of ['maxLocal', 'total']) for (const i of [0, 3]) expect(Math.abs(c[key][i])).toBeLessThan(2e-9);
      expect(Math.max(...c.internalCancellation.map(Math.abs))).toBeLessThan(2e-9);
    } else for (const key of ['maxLocal', 'total', 'internalCancellation']) expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
  if (fixture.input.streamwiseMode === 'isentropic') {
    const { gamma = 1.4, mach } = fixture.input, pInf = 1 / (gamma * mach ** 2);
    const p0Inf = pInf * (1 + .5 * (gamma - 1) * mach ** 2) ** (gamma / (gamma - 1));
    for (const s of r.flow.sections.flat(2)) {
      const m2 = s.q ** 2 / (gamma * s.p / s.rho);
      expect(Math.abs(s.p * (1 + .5 * (gamma - 1) * m2) ** (gamma / (gamma - 1)) / p0Inf - 1)).toBeLessThan(1e-10);
    }
  }
  console.log(JSON.stringify({ unknowns: r.x.length, iterations: r.history.length - 1, families: r.families, linear: r.linearDiagnostics }));
});

test('browser WASM worker selectively refines and solves all four surfaces and both wakes', async ({ page }) => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  await page.goto('/');
  const r = await page.evaluate(input => new Promise((resolve, reject) => {
    const coupled = new URL('/src/euler/streamtube-coupled.js', location.href).href;
    const refinement = new URL('/src/euler/streamtube-coupled-refinement.js', location.href).href;
    const source = `self.onmessage=async({data})=>{try{
      const {createCoupledStreamtubeBody,solveCoupledStreamtubeBody}=await import(${JSON.stringify(coupled)});
      const {refineCoupledStreamtubeBody}=await import(${JSON.stringify(refinement)});
      const parent=createCoupledStreamtubeBody(data,{edgeMatching:'section-velocity'});
      const first=solveCoupledStreamtubeBody(parent,{tolerance:1e-10});
      if(!first.converged)throw new Error('Parent did not converge');
      const refined=refineCoupledStreamtubeBody(data,parent,{initial:first.x,streamwiseFactor:1,
        normalSubdivisions:[[1,4],[4,4],[4,1]]});
      self.postMessage({parentUnknowns:parent.n,diagnostics:refined.diagnostics,
        result:solveCoupledStreamtubeBody(refined.system,{maxIterations:12,tolerance:1e-10})});
    }catch(error){self.postMessage({error:error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = e => { close(); reject(new Error(e.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(input);
  }), input);
  expect(r.parentUnknowns).toBe(349); expect(r.result.x).toHaveLength(769);
  expect(r.result.converged).toBe(true); expect(r.result.mesh.quality.valid).toBe(true);
  expect(r.result.boundaryLayer.surfaces).toHaveLength(4); expect(r.result.boundaryLayer.wakes).toHaveLength(2);
  expect(r.result.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
  for (const error of Object.values(r.result.families)) expect(error).toBeLessThan(1e-10);
  for (const [g, nodes] of r.result.flow.nodes.entries()) {
    const c = directChannelConservation({ nodes, sections: r.result.flow.sections.map(row => row[g]),
      cells: r.result.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation'])
      expect(Math.max(...c[key].map(Math.abs))).toBeLessThan(2e-9);
  }
  console.log(JSON.stringify({ unknowns: r.result.x.length, iterations: r.result.history.length - 1,
    families: r.result.families, linear: r.result.linearDiagnostics }));
});
