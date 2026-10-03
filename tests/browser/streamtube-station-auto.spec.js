// SPDX-License-Identifier: GPL-2.0-or-later
// One actual supplied-state update per policy, in a real browser worker.
// The blank page prevents the application's automatic GUI solve from running.
import fs from 'node:fs';
import { test, expect } from '@playwright/test';
import { numericalSourceHashes, changedSources, sha256 } from '../../scripts/validation/provenance.js';

const fixturePath = 'docs/coupled-current-profile-preparation/two-element-six-update/initial.json';
const fixture = JSON.parse(fs.readFileSync(fixturePath)).checkpoint;

test('fresh coupled station-auto preserves a two-body Newton update in a browser worker', async ({ page }) => {
  const sources = numericalSourceHashes(['tests/browser/streamtube-station-auto.spec.js']);
  const fixtureHash = sha256(fixturePath), buildHash = sha256('third_party/klu/BUILD.json');
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.route('**/station-ordering-check', route => route.fulfill({
    contentType: 'text/html', body: '<!doctype html><title>Coupled ordering check</title>',
  }));
  await page.goto('/station-ordering-check');
  const result = await page.evaluate(checkpoint => new Promise((resolve, reject) => {
    const module = new URL('/src/euler/streamtube-coupled-ises.js', location.href).href;
    const source = `self.onmessage=async({data:checkpoint})=>{try{
      const {solveCoupledStreamtubeIses}=await import(${JSON.stringify(module)});
      const before=JSON.stringify(checkpoint),f=checkpoint.restart,c=checkpoint.continuation;
      const options={...f.options,initialEuler:f.initialEuler,initialBL:f.initialBL,
        iterationGeometry:c.iterationGeometry,stepAcceptance:c.stepAcceptance,stagnationLimiter:c.stagnationLimiter,
        ...(c.blUpdate===undefined?{}:{blUpdate:c.blUpdate}),
        ...(c.projectionGeometry===undefined?{}:{projectionGeometry:c.projectionGeometry}),
        maxIterations:1,tolerance:1e-10};
      const published=[],started=performance.now();
      const fresh=solveCoupledStreamtubeIses(f.input,{...options,onCheckpoint:cp=>published.push(cp.continuation.linearOrdering)});
      const automatic=solveCoupledStreamtubeIses(f.input,{...options,linearOrdering:'auto'});
      const maxDifference=(a,b)=>Math.max(0,...a.map((v,i)=>Math.abs(v-b[i])));
      const nodes=r=>r.flow.nodes.flat(2).flatMap(p=>[p.x,p.y]);
      const pressure=r=>r.flow.surfaces.flatMap(s=>s.points.map(p=>p.cp));
      const forces=r=>r.flow.diagnosticForces.flatMap(f=>[f.cl,f.cd,f.cm]);
      const summarize=r=>({unknowns:r.x.length,finite:r.x.every(Number.isFinite),history:r.history,
        families:r.families,diagnosticForces:r.flow.diagnosticForces,pressure:pressure(r),
        surfaceCount:r.boundaryLayer.surfaces.length,wakeCount:r.boundaryLayer.wakes.length,
        transitionState:r.boundaryLayer.transitionState,validGrid:r.mesh.quality.valid,
        gridAcceptance:r.gridAcceptance,initialRedistribution:r.initialRedistribution,
        linearOrdering:r.linearOrdering,linearDiagnostics:r.linearDiagnostics,
        checkpointOrdering:r.checkpoint.continuation.linearOrdering});
      const elapsedMilliseconds=performance.now()-started;
      self.postMessage({fresh:summarize(fresh),automatic:summarize(automatic),published,elapsedMilliseconds,
        maximumStateDifference:maxDifference(fresh.x,automatic.x),
        maximumNodeDifference:maxDifference(nodes(fresh),nodes(automatic)),
        maximumPressureDifference:maxDifference(pressure(fresh),pressure(automatic)),
        maximumDiagnosticForceDifference:maxDifference(forces(fresh),forces(automatic)),
        unchangedInput:before===JSON.stringify(checkpoint),
        sameSolverInput:JSON.stringify(fresh.solverInput)===JSON.stringify(automatic.solverInput),
        sameCoupledOptions:JSON.stringify(fresh.coupledOptions)===JSON.stringify(automatic.coupledOptions),
        wasmResources:performance.getEntriesByType('resource').filter(r=>r.name.endsWith('/klu.wasm')).map(r=>r.name)});
    }catch(error){self.postMessage({error:error.stack??error.message});}};`;
    const url = URL.createObjectURL(new Blob([source], { type: 'text/javascript' }));
    const worker = new Worker(url, { type: 'module' });
    const close = () => { worker.terminate(); URL.revokeObjectURL(url); };
    worker.onerror = event => { close(); reject(new Error(event.message)); };
    worker.onmessage = ({ data }) => { close(); data.error ? reject(new Error(data.error)) : resolve(data); };
    worker.postMessage(checkpoint);
  }), fixture);
  expect(errors).toEqual([]);
  expect(result.wasmResources).toHaveLength(1);
  for (const r of [result.fresh, result.automatic]) {
    expect(r.unknowns).toBe(349); expect(r.finite).toBe(true);
    expect(r.history).toHaveLength(2); expect(r.surfaceCount).toBe(4); expect(r.wakeCount).toBe(2);
    expect(r.initialRedistribution.accepted).toBe(true); expect(r.initialRedistribution.resumed).not.toBe(true);
    expect(r.validGrid).toBe(true); expect(r.gridAcceptance).toBe('convex');
    expect(r.linearDiagnostics.solves).toBe(1);
    expect(r.linearDiagnostics.maxRelativeResidual).toBeLessThanOrEqual(1e-10);
    expect(Object.values(r.families).every(Number.isFinite)).toBe(true);
    expect(r.diagnosticForces).toHaveLength(2);
    expect(r.diagnosticForces.every(f => Object.values(f).every(Number.isFinite))).toBe(true);
    expect(r.pressure).toHaveLength(20); expect(r.pressure.every(Number.isFinite)).toBe(true);
  }
  expect(result.fresh.linearOrdering).toBe('station-auto');
  expect(result.fresh.checkpointOrdering).toBe('station-auto');
  expect(result.published).toEqual(['station-auto', 'station-auto']);
  expect(result.automatic.linearOrdering).toBeUndefined(); expect(result.automatic.checkpointOrdering).toBeUndefined();
  const d = result.fresh.linearDiagnostics.iterations[0];
  expect(d.ordering).toBe('given'); expect(d.pivotTolerance).toBe(.001);
  expect(d.attempts).toHaveLength(1); expect(d.attempts[0].btf).toBe(false);
  expect(d.stationOrdering.matchedNonzeroDiagonals).toBe(349);
  expect(d.stationPolicy).toMatchObject({ attempted: true, accepted: true, fallback: false,
    recommendation: 'station-auto', selected: { ordering: 'given', pivotTolerance: .001 } });
  expect(Object.values(d.stationPolicy.timings).every(v => Number.isFinite(v) && v >= 0)).toBe(true);
  expect(result.fresh.history[1].backtracks).toBe(result.automatic.history[1].backtracks);
  expect(Math.abs(result.fresh.history[1].step - result.automatic.history[1].step)).toBeLessThan(1e-12);
  expect(result.fresh.transitionState).toEqual(result.automatic.transitionState);
  expect(result.maximumStateDifference).toBeLessThan(1e-10);
  expect(result.maximumNodeDifference).toBeLessThan(1e-10);
  expect(result.maximumPressureDifference).toBeLessThan(1e-10);
  expect(result.maximumDiagnosticForceDifference).toBeLessThan(1e-10);
  expect(result.unchangedInput).toBe(true); expect(result.sameSolverInput).toBe(true); expect(result.sameCoupledOptions).toBe(true);
  expect(changedSources(sources)).toEqual([]);
  expect(sha256(fixturePath)).toBe(fixtureHash); expect(sha256('third_party/klu/BUILD.json')).toBe(buildHash);
  fs.writeFileSync('docs/linear-ordering-integration/station-auto-browser.json', JSON.stringify({
    date: new Date().toISOString(), sourceHashes: sources, fixture: { path: fixturePath, sha256: fixtureHash }, buildHash,
    scope: 'Real Chromium worker, shipped WASM, fresh supplied two-body state: one coupled Newton update for each ordering; four surface boundary layers and two wakes. No panel/mesh initialization, full convergence claim, or large trajectory.',
    ...result,
  }, null, 2) + '\n');
  console.log(JSON.stringify({ unknowns: result.fresh.unknowns, elapsedMilliseconds: result.elapsedMilliseconds,
    maximumStateDifference: result.maximumStateDifference, maximumNodeDifference: result.maximumNodeDifference,
    maximumPressureDifference: result.maximumPressureDifference,
    linearRelativeResidual: result.fresh.linearDiagnostics.maxRelativeResidual, stationPolicy: d.stationPolicy }));
});
