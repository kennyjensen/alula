import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
const maximum = a => Math.max(0, ...Array.from(a, Math.abs));

for (const streamwiseMode of ['momentum','isentropic']) for (const wakeGeometry of ['centerline', 'independent-banks'])
for (const edgeMatching of streamwiseMode==='isentropic'?['section-velocity','section-velocity-distance']:['section-velocity'])
test(`${streamwiseMode}, ${wakeGeometry}, ${edgeMatching}: local coupled refinement verifies station/BL/wake transfer and the full Jacobian`, async t => {
  // A controlled attached, fixed-trip fixture. Its lower Reynolds number
  // isolates station transfer from unsupported free-transition migration.
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), wakeGeometry, streamwiseMode };
  const r = solveCoupledStreamtubeIses(input, { reynolds: 1e5, edgeMatching, tolerance: 1e-10, stepAcceptance: 'admissible' });
  assert.equal(r.converged, true, r.reason);
  const ne = r.x.length - 4 * r.boundaryLayer.stations.length;
  const source = createCoupledStreamtubeBody(r.solverInput, { ...r.coupledOptions, initialEuler: r.flow, initialBL: r.x.slice(ne) });
  const x = source.initial.slice(), inputBefore = structuredClone(r.solverInput), old = source.evaluate(x);
  const counts = Array(source.euler.layout.nx).fill(1);
  for (const body of source.euler.layout.bodies) {
    counts[body.leadingIndex - 1] = 2;
    counts[body.leadingIndex] = 2;
    counts[body.leadingIndex + 1] = 3;
    counts[body.trailingIndex] = 2;
  }
  const retained = [0]; for (const n of counts) retained.push(retained.at(-1) + n);
  const normal = [[1, 4], [4, 4], [4, 1]], positions = [[0, 1, 5], [0, 4, 8], [0, 4, 5]];
  const refined = refineCoupledStreamtubeBody(r.solverInput, source, {
    streamwiseSubdivisions: counts, normalSubdivisions: normal, normalInterpolation: 'streamfunction-quadratic' });
  const child = refined.system, value = child.evaluate(child.initial);
  assert.deepEqual(r.solverInput, inputBefore); assert.deepEqual(source.initial, x);
  assert.equal(child.euler.layout.nx, retained.at(-1));
  assert.deepEqual(refined.diagnostics.retainedStreamwiseStations, retained);
  assert.equal(refined.diagnostics.streamwiseFactor, null);
  assert.equal(child.bl.surfaces.length, 4); assert.equal(child.bl.wakes.length, 2);
  assert.equal(child.conditions.edgeMatching,edgeMatching);
  assert.equal(child.admissible(child.initial), true);
  assert.equal(directStreamtubeVolumeGeometry(value.outer.nodes).valid, true);
  assert.deepEqual(old.outer.captured, value.outer.captured);
  for (const [name, col] of Object.entries(source.euler.layout.globals)) {
    const after = child.euler.layout.globals[name];
    if (Array.isArray(col)) col.forEach((c,b) => { if (c !== null) assert.equal(child.initial[after[b]], x[c]); });
    else assert.equal(child.initial[after], x[col]);
  }
  for (let g = 0; g < old.outer.nodes.length; g++) {
    old.outer.nodes[g].forEach((row,i) => row.forEach((p,j) => {
      const q = value.outer.nodes[g][retained[i]][positions[g][j]];
      // Unequal wake subdivisions change the centerline mode's secant
      // normal. Its prescribed-offset banks must be reconstructed; the
      // independent-bank mode retains the physical nodes themselves.
      const wakeAdjacent = g > 0 && i > source.euler.layout.bodies[g-1].trailingIndex
        || g < source.euler.layout.elements && i > source.euler.layout.bodies[g].trailingIndex;
      const allowed = wakeGeometry === 'centerline' && wakeAdjacent ? Math.max(...[0, row.length-1].map(k => {
        const a=row[k], b=value.outer.nodes[g][retained[i]][positions[g][k]]; return Math.hypot(a.x-b.x,a.y-b.y);
      })) : 0;
      assert.ok(Math.hypot(q.x-p.x,q.y-p.y) < allowed+2e-12, `retained node ${g}/${i}/${j}`);
    }));
    old.outer.allocation.groups[g].forEach((tube,j) => {
      const sum = value.outer.allocation.groups[g].slice(positions[g][j],positions[g][j+1]).reduce((a,t) => a+t.massFlow,0);
      assert.ok(Math.abs(sum/tube.massFlow-1) < 1e-14);
    });
  }
  if (wakeGeometry === 'centerline') for (const [b,body] of source.euler.layout.bodies.entries()) {
    const center = i => {
      const a=value.outer.nodes[b][i].at(-1),q=value.outer.nodes[b+1][i][0];return {x:.5*(a.x+q.x),y:.5*(a.y+q.y)};
    };
    for (let i=body.trailingIndex+1;i<=source.euler.layout.nx;i++) {
      const next=retained[i],p=center(next),a=old.outer.nodes[b][i].at(-1),q=old.outer.nodes[b+1][i][0];
      assert.ok(Math.hypot(p.x-.5*(a.x+q.x),p.y-.5*(a.y+q.y))<2e-12);
      const left=center(next-1),right=center(Math.min(child.euler.layout.nx,next+1));
      const lower=value.outer.nodes[b][next].at(-1),upper=value.outer.nodes[b+1][next][0];
      const dx=upper.x-lower.x,dy=upper.y-lower.y,tx=right.x-left.x,ty=right.y-left.y;
      assert.ok(Math.abs(dx*tx+dy*ty)<1e-12);
      assert.ok(Math.abs(Math.hypot(dx,dy)-Math.hypot(a.x-q.x,a.y-q.y))<2e-12);
    }
  }
  for (const station of source.bl.stations) {
    const next = child.bl.stations.find(s => s.kind === station.kind && s.body === station.body && s.side === station.side && s.i === retained[station.i]);
    assert.ok(next);
    for (const key of ['aux','theta','deltaStar','ue']) assert.equal(value.layers.states[next.id][key],old.layers.states[station.id][key]);
  }
  assert.deepEqual(child.bl.surfaces.map(s=>s.tripParameter),source.bl.surfaces.map(s=>s.tripParameter));
  for (const [b,body] of source.euler.layout.bodies.entries()) {
    assert.equal(child.euler.layout.bodies[b].leadingIndex,retained[body.leadingIndex]);
    assert.equal(child.euler.layout.bodies[b].trailingIndex,retained[body.trailingIndex]);
  }
  const jac = child.jacobian(child.initial), h = 2e-6;
  const d = child.initial.map((v,i)=>(i < child.ne ? .001 : Math.max(.01,Math.abs(v)))*Math.sin(.43*i+.2));
  const exact = sparseProduct(jac,d), plus = child.residual(child.initial.map((v,i)=>v+h*d[i])), minus = child.residual(child.initial.map((v,i)=>v-h*d[i]));
  const error = maximum(exact.map((v,i)=> { const fd=(plus[i]-minus[i])/(2*h); return (fd-v)/Math.max(1,Math.abs(fd),Math.abs(v)); }));
  assert.ok(error < 5e-6, `mixed Jacobian error ${error}`);
  await t.test('complete coupled recovery on this deliberately coarse fixture',
    streamwiseMode==='momentum' ? {todo:'The retained local-refinement control stalls at crossed half-volumes; full recovery remains an open robustness gate.'} : {}, () => {
  const solved = solveCoupledStreamtubeIses({...refined.input,geometryDomain:'positive-simple'},{...refined.options,
    initialEuler:refined.initialEuler,initialBL:refined.initialBL,maxIterations:16,tolerance:1e-10,
    stepAcceptance:'admissible',iterationGeometry:'ises-sampled'});
  t.diagnostic(JSON.stringify({parentUnknowns:source.n,unknowns:child.n,jacobianError:error,converged:solved.converged,reason:solved.reason,families:solved.families,lastRejected:solved.lastRejectedStep}));
  if (!solved.converged) {
    const end = createCoupledStreamtubeBody(solved.solverInput,{...solved.coupledOptions,
      initialEuler:solved.flow,initialBL:solved.x.slice(child.ne)}), z=end.initial, base=end.evaluate(z);
    const direction=solveSparseDirect(end.jacobian(z),base.residual.map(v=>-v)).x;
    const trials=[1,.01,1e-4,1e-6].map(h=> {
      const candidate=z.map((v,i)=>v+h*direction[i]);
      try { const v=end.evaluate(candidate); return {h,admissible:end.admissible(candidate),families:v.families,
        minimumConstraint:Math.min(...end.constraintValues(candidate))}; }
      catch(error) {return {h,error:error.message};}
    });
    t.diagnostic(JSON.stringify({frozenOrdinaryNewton:trials,history:solved.history.map(h=>({iteration:h.iteration,residual:h.residual,step:h.step}))}));
  }
  assert.equal(solved.converged,true,solved.reason);
  assert.equal(directStreamtubeVolumeGeometry(solved.flow.nodes).valid,true);
  solved.flow.nodes.forEach((nodes,g)=> {
    const c=directChannelConservation({nodes,sections:solved.flow.sections.map(row=>row[g]),cells:solved.flow.cells.map(row=>row[g])});
    for(const key of ['maxLocal','total','internalCancellation'])
      for(const k of streamwiseMode==='momentum' || key==='internalCancellation' ? [0,1,2,3] : [0,3]) assert.ok(Math.abs(c[key][k]) < 2e-9);
  });
  });
});
