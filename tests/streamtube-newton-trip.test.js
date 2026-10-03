import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

test('ordinary Newton crosses the retained default material-trip event with all Euler/BL/wake equations active', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('fixtures/default-newton-trip-crossing.json', import.meta.url)));
  const build = () => createCoupledStreamtubeBody(fixture.input, { ...fixture.options,
    initialEuler: fixture.initialEuler, initialBL: Float64Array.from(fixture.initialBL) });
  const system = build(), x = system.initial.slice(), initial = system.evaluate(x);
  for (const key in fixture.families) assert.ok(Math.abs(initial.families[key] - fixture.families[key]) < 1e-10);
  const disabled = solveCoupledStreamtubeBody(build(), { maxIterations: 1, transitionEvents: false });
  const trips = system.bl.surfaces.map(s => [s.tripArc, s.tripParameter]);
  const r = solveCoupledStreamtubeBody(system, { maxIterations: 1 });
  assert.equal(r.stepMethod, 'newton'); assert.equal(r.transitionEvents, true);
  assert.equal(r.history.length, 2); assert.equal(r.history[1].activeChange, true);
  assert.equal(r.history[1].meritComparable, false);
  assert.ok(r.history[1].step >= .0625); assert.ok(r.history[1].step > 16 * disabled.history[1].step);
  assert.equal(r.linearDiagnostics.solves, 1); assert.equal(r.converged, false);
  assert.equal(system.bl.surfaces.length, 4); assert.equal(system.bl.wakes.length, 2);
  assert.equal(system.admissible(r.x), true); assert.equal(directStreamtubeVolumeGeometry(r.flow.nodes).valid, true);
  assert.deepEqual(system.bl.surfaces.map(s => [s.tripArc, s.tripParameter]), trips);
  assert.deepEqual(system.bl.activeTargets(r.x.subarray(0, system.ne)).filter(t => t.from !== t.to), []);
  const v = system.evaluate(r.x), j = system.jacobian(r.x), h = 2e-6;
  const d = r.x.map((q, i) => (i < system.ne ? .001 : Math.max(.01, Math.abs(q))) * Math.sin(i * .43 + .2));
  const exact = sparseProduct(j, d), plus = system.residual(r.x.map((q, i) => q + h * d[i])), minus = system.residual(r.x.map((q, i) => q - h * d[i]));
  const error = Math.max(...exact.map((q, i) => { const fd = (plus[i] - minus[i]) / (2 * h); return Math.abs(fd - q) / Math.max(1, Math.abs(fd), Math.abs(q)); }));
  assert.ok(error < 5e-6, `post-event Jacobian error ${error}`);
  t.diagnostic(JSON.stringify({ unknowns: system.n, disabledStep: disabled.history[1].step,
    acceptedStep: r.history[1].step, changes: r.history[1].changes, families: v.families, jacobianError: error }));
});
