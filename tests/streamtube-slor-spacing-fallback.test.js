import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { createEllipticStreamtubeGrid, smoothEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { createBoundaryStretchControl } from '../src/geometry/boundary-stretch-control.js';

// Control just the paired solver's terminal status. The fallback, coordinates,
// mass fractions and scalar tridiagonal SLOR below are the actual implementation.
async function adapterWithTerminal(reason, termination) {
  const url = new URL('../src/euler/streamtube-elliptic-initializer.js', import.meta.url);
  const encode = text => 'data:text/javascript;base64,' + Buffer.from(text).toString('base64');
  const bridge = encode(`export function smoothPairedBoundaryGrid(system) {
    return {nodes:structuredClone(system.initial),converged:false,reason:${JSON.stringify(reason)},
      history:[{iteration:0}],coordinateEquations:system.coordinateEquations,
      termination:${JSON.stringify(termination)}};
  }`);
  const text = fs.readFileSync(url,'utf8').replace(/from '(\.\.?\/[^']+)'/g, (_, spec) =>
    `from '${spec.endsWith('/paired-boundary-slor.js') ? bridge : new URL(spec,url).href}'`);
  return (await import(encode(text))).relaxStreamtubeInitialGrid;
}

test('a rejected wall-angle line retries spacing and retains the harmonic fallback when needed', async () => {
  const relax = await adapterWithTerminal('No admissible decreasing Newton step for eta rows 1, 2.',
    { origin: 'solver', termination: 'no-admissible-decreasing-line-step' });
  const system=createStreamtubeBodySystem(intrinsicBodyFixture({elements:2,bodySegments:8,tubes:3}));
  const before=system.decode(system.initial), chart=system.geometryChart(), snapshots=[];
  const result=relax({system},{seed:'supplied',boundaryControl:'wall-angle',maxSweeps:2,
    onSweep:(h,nodes)=>snapshots.push({h,nodes})});
  assert.deepEqual(system.decode(system.initial),before);assert.deepEqual(system.geometryChart(),chart);
  for(let g=0;g<before.nodes.length;g++) {
    const report=result.regions[g], seed=before.nodes[g], massFlows=before.allocation.groups[g].map(t=>t.massFlow);
    assert.equal(report.spacingFallback.attempted,true);
    assert.equal(report.spacingFallback.pointwiseWallAngleEnforced,false);
    assert.equal(report.spacingFallback.originalAttempt.converged,false);
    if (report.harmonicFallback) {
      assert.equal(report.harmonicFallback.attempted,true);
      assert.equal(report.coordinateEquations.xi,'Laplace');
    } else assert.match(report.coordinateEquations.xi,/prescribed metric-scaled streamwise stretch/);
    assert.equal(report.coordinateEquations.eta,'Laplace');
    const xi=result.stationCoordinate.xi;
    const base=createEllipticStreamtubeGrid({nodes:seed,massFlows,streamwiseCoordinates:xi,discretization:'giles-1985'});
    const fixed=createEllipticStreamtubeGrid({nodes:seed,massFlows,streamwiseCoordinates:xi,discretization:'giles-1985',
      streamwiseStretch:createBoundaryStretchControl({nodes:seed,xi,eta:base.eta,metric:'polygon-arc'}).values});
    const direct=smoothEllipticStreamtubeGrid(report.harmonicFallback ? base : fixed,
      {maxSweeps:2,omega:1,...(report.harmonicFallback ? {requireConvex:true} : {})});
    assert.deepEqual(result.nodes[g],direct.nodes);
    for(let i=0;i<seed.length;i++)for(let j=0;j<seed[i].length;j++)
      if(!i||i===seed.length-1||!j||j===seed[i].length-1)assert.deepEqual(result.nodes[g][i][j],seed[i][j]);
    const first=snapshots.find(s=>s.h.region===g);
    assert(first.h.spacingFallback);assert.deepEqual(first.nodes,seed);
    first.nodes[1][1].x+=100;assert.notEqual(first.nodes[1][1].x,result.nodes[g][1][1].x);
  }
});

test('an iteration limit or unrelated error does not start a different smoothing solve', async () => {
  const system=createStreamtubeBodySystem(intrinsicBodyFixture({elements:2,bodySegments:8,tubes:3}));
  for(const reason of ['sweep limit','Invalid starting state: folded grid','observer failed']) {
    const relax=await adapterWithTerminal(reason);
    const r=relax({system},{seed:'supplied',boundaryControl:'wall-angle',maxSweeps:0});
    assert.equal(r.converged,false);
    assert(r.regions.every(region=>region.reason===reason&&region.spacingFallback===undefined));
    assert.deepEqual(r.nodes,system.decode(system.initial).nodes);
  }
});
