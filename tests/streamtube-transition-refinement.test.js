import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { redistributeCoupledSurfaceStations } from '../src/euler/tests/streamtube-coupled-redistribution.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';

test('subdividing natural-transition intervals preserves physical fields and reconverges four BLs and two wakes', t => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const parent = createCoupledStreamtubeBody(input, { transitionMode: 'automatic', edgeMatching: 'section-velocity' });
  const root = solveCoupledStreamtubeBody(parent, { tolerance: 1e-10, maxIterations: 12 });
  assert.equal(root.converged, true, root.reason);
  const old = parent.evaluate(root.x), x = root.x.slice(), phases = parent.bl.snapshotActive(), beforeInput = structuredClone(input);
  const counts = Array(parent.euler.layout.nx).fill(1);
  for (const q of old.layers.transitions) {
    assert.equal(q.kind, 'natural'); counts[parent.bl.stations[q.id].i - 1] = 4;
  }
  const r = refineCoupledStreamtubeBody(input, parent, { initial: root.x, streamwiseSubdivisions: counts, normalFactor: 1 });
  const child = r.system, v = child.evaluate(child.initial), retained = r.diagnostics.retainedStreamwiseStations;
  assert.deepEqual(input, beforeInput); assert.deepEqual(root.x, x); assert.deepEqual(parent.bl.snapshotActive(), phases);
  assert.deepEqual(parent.evaluate(root.x).families, old.families);
  assert.deepEqual(v.outer.allocation, old.outer.allocation);
  assert.deepEqual(child.bl.trips, parent.bl.trips);
  for (const key of ['reynolds', 'ncrit', 'mach', 'edgeMatching']) assert.equal(child.conditions[key], parent.conditions[key]);
  for (const p of parent.bl.stations) {
    const q = child.bl.stations.find(q => q.body === p.body && q.side === p.side && q.kind === p.kind && q.i === retained[p.i]);
    assert.ok(q);
    for (const key of ['theta', 'deltaStar', 'ue']) assert.ok(Math.abs(old.layers.states[p.id][key] - v.layers.states[q.id][key]) < 1e-14);
  }
  old.outer.nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    // Solid offsets and the outer boundaries are retained exactly. In this
    // centerline-wake fixture, inserting unequal neighboring intervals changes
    // the wake normal and therefore its banks and nearby interior nodes.
    const solid = j === 0 && g > 0 && parent.euler.layout.active(g - 1, i)
      || j === row.length - 1 && g < parent.euler.layout.elements && parent.euler.layout.active(g, i);
    const outer = g === 0 && j === 0 || g === old.outer.nodes.length - 1 && j === row.length - 1;
    if (solid || outer) {
      const q = v.outer.nodes[g][retained[i]][j]; assert.ok(Math.hypot(q.x - p.x, q.y - p.y) < 2e-12);
    }
  })));
  for (const wake of parent.bl.wakes) for (const id of wake.ids.slice(1)) {
    const p = parent.bl.stations[id], i = retained[p.i], b = wake.body;
    const a = old.outer.nodes[b][p.i].at(-1), z = old.outer.nodes[b + 1][p.i][0];
    const l = v.outer.nodes[b][i].at(-1), u = v.outer.nodes[b + 1][i][0];
    assert.ok(Math.hypot(l.x + u.x - a.x - z.x, l.y + u.y - a.y - z.y) < 2e-12);
    assert.ok(Math.abs(Math.hypot(u.x - l.x, u.y - l.y) - old.layers.states[id].deltaStar * child.euler.conditions.lengthScale) < 2e-12);
  }
  assert.equal(child.bl.surfaces.length, 4); assert.equal(child.bl.wakes.length, 2);
  const matrix = child.jacobian(child.initial), errors = [];
  for (const direction of [0, 1]) {
    const d = child.initial.map((v, i) => (i < child.ne ? .001 : Math.max(.01, Math.abs(v))) * Math.sin((.43 + .17 * direction) * i + .2));
    const exact = sparseProduct(matrix, d), h = 1e-6;
    const central = step => {
      const p = child.residual(child.initial.map((v, i) => v + step * d[i])), m = child.residual(child.initial.map((v, i) => v - step * d[i]));
      return p.map((v, i) => (v - m[i]) / (2 * step));
    };
    const full = central(h), half = central(h / 2);
    errors.push(Math.max(...exact.map((v, i) => {
      const fd = (4 * half[i] - full[i]) / 3; return Math.abs(v - fd) / Math.max(1, Math.abs(v), Math.abs(fd));
    })));
  }
  assert.ok(errors.every(e => e < 5e-6), JSON.stringify(errors));
  const solved = solveCoupledStreamtubeBody(child, { tolerance: 1e-10, maxIterations: 16 });
  assert.equal(solved.converged, true, solved.reason); assert.equal(solved.mesh.quality.valid, true);
  assert.ok(Object.values(solved.families).every(v => v < 1e-10));
  solved.flow.nodes.forEach((nodes, g) => {
    const c = directChannelConservation({ nodes, sections: solved.flow.sections.map(row => row[g]), cells: solved.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(c[key].every(v => Math.abs(v) < 2e-9));
  });
  t.diagnostic(JSON.stringify({ unknowns: child.n, errors, iterations: solved.history.length - 1, families: solved.families }));
});

test('nonidentity automatic surface redistribution supplies an evaluable guess and reconverges', t => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), streamwiseMode: 'isentropic' };
  const parent = createCoupledStreamtubeBody(input, { transitionMode: 'automatic', edgeMatching: 'section-velocity' });
  const root = solveCoupledStreamtubeBody(parent, { tolerance: 1e-10, maxIterations: 12 });
  assert.equal(root.converged, true, root.reason);
  const fractions = structuredClone(parent.euler.fractions), before = root.x.slice(), phases = parent.bl.snapshotActive();
  for (const body of fractions) for (const side of ['upper', 'lower']) body[side][2] += .02 * (body[side][3] - body[side][2]);
  const r = redistributeCoupledSurfaceStations(input, parent, { initial: root.x, surfaceFractions: fractions });
  assert.deepEqual(root.x, before); assert.deepEqual(parent.bl.snapshotActive(), phases);
  assert.equal(r.system.admissible(r.system.initial), true);
  const solved = solveCoupledStreamtubeBody(r.system, { tolerance: 1e-10, maxIterations: 12 });
  assert.equal(solved.converged, true, solved.reason); assert.equal(solved.mesh.quality.valid, true);
  assert.ok(Object.values(solved.families).every(v => v < 1e-10));
  assert.equal(solved.boundaryLayer.surfaces.length, 4); assert.equal(solved.boundaryLayer.wakes.length, 2);
  t.diagnostic(JSON.stringify({ iterations: solved.history.length - 1, families: solved.families }));
});
