import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const maximum = a => Math.max(0, ...Array.from(a, Math.abs));
const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), wakeGeometry: 'independent-banks' };
const r = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', tolerance: 1e-10, stepAcceptance: 'admissible' });
assert.equal(r.converged, true);
const ne = r.x.length - 4 * r.boundaryLayer.stations.length;
const parent = () => createCoupledStreamtubeBody(r.solverInput, { ...r.coupledOptions, initialEuler: r.flow, initialBL: r.x.slice(ne) });

for (const streamwiseFactor of [1, 2]) test(`independent-bank refinement preserves both physical banks with streamwise factor ${streamwiseFactor}`, t => {
  const source = parent(), x = source.initial.slice(), before = source.evaluate(x);
  const counts = [[1, 2], [2, 1], [1, 1]], positions = counts.map(row => [0, ...row.map((_, j) => row.slice(0, j + 1).reduce((a, b) => a + b))]);
  const refined = refineCoupledStreamtubeBody(r.solverInput, source, { streamwiseFactor, normalSubdivisions: counts });
  const child = refined.system, after = child.evaluate(child.initial);
  assert.deepEqual(source.initial, x); assert.equal(child.admissible(child.initial), true);
  assert.equal(child.euler.conditions.wakeGeometry, 'independent-banks');
  assert.deepEqual(after.outer.captured, before.outer.captured);
  for (let g = 0; g < before.outer.nodes.length; g++) {
    before.outer.nodes[g].forEach((row, i) => row.forEach((p, j) => {
      const q = after.outer.nodes[g][streamwiseFactor * i][positions[g][j]];
      assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 2e-12, `parent node ${g}/${i}/${j}`);
    }));
    before.outer.allocation.groups[g].forEach((tube, j) => {
      const sum = after.outer.allocation.groups[g].slice(positions[g][j], positions[g][j + 1]).reduce((s, q) => s + q.massFlow, 0);
      assert.ok(Math.abs(sum / tube.massFlow - 1) < 1e-14);
    });
  }
  for (const station of source.bl.stations) {
    const next = child.bl.stations.find(s => s.kind === station.kind && s.side === station.side && s.body === station.body && s.i === streamwiseFactor * station.i);
    for (const key of ['aux', 'theta', 'deltaStar', 'ue']) assert.equal(after.layers.states[next.id][key], before.layers.states[station.id][key]);
  }
  assert.deepEqual(child.bl.surfaces.map(b => b.tripParameter), source.bl.surfaces.map(b => b.tripParameter));
  if (streamwiseFactor === 1) assert.ok(after.outer.diagnostics.residualByFamily.wakeGap < 1e-12);
  assert.equal(directStreamtubeVolumeGeometry(after.outer.nodes).valid, true);
  // Streamwise refinement introduces new BL stations; a prolonged state
  // must not inherit convergence. At these material trips the attempted
  // full recovery encounters the unsupported natural-transition branch.
  // Keep transfer/Jacobian coverage for it, and separately verify the
  // supported normal-refinement root's complete regional conservation.
  const solved = solveCoupledStreamtubeBody(child, { maxIterations: streamwiseFactor === 1 ? 12 : 0, tolerance: 1e-10 });
  if (streamwiseFactor === 1) {
    assert.equal(solved.converged, true, solved.reason); assert.ok(maximum(solved.residual) < 1e-10);
    for (let g = 0; g < solved.flow.nodes.length; g++) {
      const c = directChannelConservation({ nodes: solved.flow.nodes[g], sections: solved.flow.sections.map(row => row[g]), cells: solved.flow.cells.map(row => row[g]) });
      for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(maximum(c[key]) < 2e-9);
    }
  } else {
    assert.equal(solved.converged, false); assert.equal(solved.mesh.initialization.flowSolved, false);
    assert.ok(solved.families.euler > 1e-4 && solved.families.boundaryLayer > 1e-4);
  }
  const state = child.rebase(solved.x), jac = child.jacobian(state), h = 2e-6;
  const d = state.map((v, i) => (i < child.ne ? .001 : Math.max(.01, Math.abs(v))) * Math.sin(.43 * i + .2));
  const exact = sparseProduct(jac, d), a = child.residual(state.map((v, i) => v + h * d[i])), b = child.residual(state.map((v, i) => v - h * d[i]));
  const error = maximum(exact.map((v, i) => { const fd = (a[i] - b[i]) / (2 * h); return (fd - v) / Math.max(1, Math.abs(fd), Math.abs(v)); }));
  assert.ok(error < 5e-6);
  t.diagnostic(JSON.stringify({ unknowns: child.n, updates: solved.history.length - 1, families: solved.families, jacobianError: error }));
});
