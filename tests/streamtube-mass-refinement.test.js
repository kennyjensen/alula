import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { sparseProduct } from '../src/numerics/sparse.js';
const maximum = a => Math.max(0, ...Array.from(a, Math.abs));

for (const wakeGeometry of ['centerline', 'independent-banks']) test(`${wakeGeometry}: physical mass interpolation retains all coupled parent data and passes full Jacobian/conservation checks`, t => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), wakeGeometry };
  const r = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', tolerance: 1e-10, stepAcceptance: 'admissible' });
  assert.equal(r.converged, true);
  const ne = r.x.length - 4 * r.boundaryLayer.stations.length;
  const source = createCoupledStreamtubeBody(r.solverInput, { ...r.coupledOptions, initialEuler: r.flow, initialBL: r.x.slice(ne) });
  const x = source.initial.slice(), before = source.evaluate(x), counts = [[1, 4], [4, 1], [1, 1]], positions = [[0, 1, 5], [0, 4, 5], [0, 1, 2]];
  const refined = refineCoupledStreamtubeBody(r.solverInput, source, { streamwiseFactor: 1, normalSubdivisions: counts, normalInterpolation: 'streamfunction-quadratic' });
  const child = refined.system, after = child.evaluate(child.initial);
  assert.deepEqual(source.initial, x); assert.equal(child.admissible(child.initial), true);
  assert.equal(refined.diagnostics.normalInterpolation, 'streamfunction-quadratic');
  assert.ok(refined.diagnostics.massInterpolation.some(d => d.maximumFractionChange > .01));
  assert.equal(directStreamtubeVolumeGeometry(after.outer.nodes).valid, true);
  for (let g = 0; g < before.outer.nodes.length; g++) {
    before.outer.nodes[g].forEach((row, i) => row.forEach((p, j) => {
      const q = after.outer.nodes[g][i][positions[g][j]]; assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 2e-12);
    }));
    before.outer.allocation.groups[g].forEach((tube, j) => {
      const sum = after.outer.allocation.groups[g].slice(positions[g][j], positions[g][j + 1]).reduce((s, t) => s + t.massFlow, 0);
      assert.ok(Math.abs(sum / tube.massFlow - 1) < 1e-14);
    });
  }
  assert.deepEqual(before.layers.states, after.layers.states);
  assert.deepEqual(source.bl.surfaces.map(s => s.tripParameter), child.bl.surfaces.map(s => s.tripParameter));
  const jac = child.jacobian(child.initial), h = 2e-6, d = child.initial.map((v, i) => (i < child.ne ? .001 : Math.max(.01, Math.abs(v))) * Math.sin(.43 * i + .2));
  const exact = sparseProduct(jac, d), plus = child.residual(child.initial.map((v, i) => v + h * d[i])), minus = child.residual(child.initial.map((v, i) => v - h * d[i]));
  const error = maximum(exact.map((v, i) => { const fd = (plus[i] - minus[i]) / (2 * h); return (fd - v) / Math.max(1, Math.abs(fd), Math.abs(v)); }));
  assert.ok(error < 5e-6);
  const solved = solveCoupledStreamtubeBody(child, { maxIterations: 12, tolerance: 1e-10 });
  assert.equal(solved.converged, true, solved.reason);
  solved.flow.nodes.forEach((nodes, g) => {
    const c = directChannelConservation({ nodes, sections: solved.flow.sections.map(row => row[g]), cells: solved.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(maximum(c[key]) < 2e-9);
  });
  t.diagnostic(JSON.stringify({ unknowns: child.n, interpolation: refined.diagnostics.massInterpolation, jacobianError: error, families: solved.families }));
});
