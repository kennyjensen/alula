import test from 'node:test';
import assert from 'node:assert/strict';
import { coupledDownstreamControl } from '../scripts/validation/coupled-downstream-control.js';
import { truncateCoupledDownstreamDomain } from '../src/euler/tests/streamtube-coupled-truncation.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

test('exact prefix truncation preserves multielement states and solves the changed outlet with all four BLs and two wakes', t => {
  const { input, source, mapped, seed, result } = coupledDownstreamControl(), before = source.evaluate(source.initial);
  // BL values are copied exactly; physical wake-bank nodes are reconstructed
  // by the new geometry chart. Allow only a few floating-point roundoff units
  // at the coordinate scale, including cancellation near y=0. The captured
  // pre-transition source has the same 2.17e-19 reconstruction difference.
  const coordinateScale = before.outer.nodes.reduce((m, grid) => grid.reduce((m, row) => row.reduce(
    (m, p) => Math.max(m, Math.abs(p.x), Math.abs(p.y)), m), m), source.euler.conditions.lengthScale);
  assert.ok(mapped.diagnostics.nodeError <= 8 * Number.EPSILON * coordinateScale);
  assert.equal(mapped.diagnostics.blError, 0);
  assert.deepEqual(seed.outer.allocation, before.outer.allocation);
  assert.deepEqual(mapped.system.bl.surfaces.map(b => b.tripParameter), source.bl.surfaces.map(b => b.tripParameter));
  assert.deepEqual(mapped.system.initial.slice(0, mapped.system.euler.layout.densityCount), source.initial.slice(0, mapped.system.euler.layout.densityCount));
  assert.ok(mapped.system.n < source.n);
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.boundaryLayer.surfaces.length, 4); assert.equal(result.boundaryLayer.wakes.length, 2);
  const geometry = directStreamtubeVolumeGeometry(result.flow.nodes);
  assert.equal(geometry.valid, true); assert.deepEqual(geometry.concavePrimal, []);
  result.flow.nodes.forEach((nodes, g) => {
    const c = directChannelConservation({ nodes, sections: result.flow.sections.map(row => row[g]), cells: result.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) for (const k of key === 'internalCancellation' ? [0, 1, 2, 3] : [0, 3])
      assert.ok(Math.abs(c[key][k]) < 2e-9, `${g}/${key}/${k}`);
  });
  const s = mapped.system, x = result.x, d = x.map((v, i) => (i < s.ne ? .001 : Math.max(.01, Math.abs(v))) * Math.sin(.43 * i + .2));
  const exact = sparseProduct(s.jacobian(x), d), h = 2e-6;
  const plus = s.residual(x.map((v, i) => v + h * d[i])), minus = s.residual(x.map((v, i) => v - h * d[i]));
  const error = Math.max(...exact.map((v, i) => {
    const fd = (plus[i] - minus[i]) / (2 * h); return Math.abs(fd - v) / Math.max(1, Math.abs(fd), Math.abs(v));
  })); assert.ok(error < 5e-6);
  const oldInput = structuredClone(input), oldInitial = source.initial.slice();
  const identity = truncateCoupledDownstreamDomain(input, source, { endIndex: source.euler.layout.nx });
  assert.deepEqual(identity.system.evaluate(identity.system.initial).families, before.families);
  for (const endIndex of [null, 1.5, source.euler.layout.nx + 1, Math.max(...source.euler.layout.bodies.map(b => b.trailingIndex))])
    assert.throws(() => truncateCoupledDownstreamDomain(input, source, { endIndex }), /intervals/);
  assert.deepEqual(input, oldInput); assert.deepEqual(source.initial, oldInitial);
  t.diagnostic(JSON.stringify({ parentUnknowns: source.n, unknowns: s.n, families: result.families, iterations: result.history.length - 1, jacobianError: error }));
});
