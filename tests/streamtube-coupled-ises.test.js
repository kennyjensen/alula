import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { proposeCoupledDensityNewton } from '../src/euler/streamtube-density-newton.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody, coupledStreamtubeResult, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

test('Giles viscous caps apply one scalar to density, BL and auxiliary updates', () => {
  const system = { n: 5, ne: 1, bl: { stations: [{ id: 0 }], thicknesses: x => x }, euler: {
    layout: { n: 1, densityCount: 1, positions: [], globals: { stagnation: [] } }, curves: [],
    decode: () => ({ nodes: [] }), setDisplacement: () => {},
  } };
  const x = Float64Array.of(Math.log(1.2), .1, 1, 3, 2), d = Float64Array.of(.4, .02, -2, -1, 1);
  const r = proposeCoupledDensityNewton(system, x, d);
  assert.equal(r.step, .3); assert.equal(r.limiter.variable, 'theta');
  assert.ok(Math.abs(Math.exp(r.x[0]) - 1.344) < 1e-14);
  for (let i = 1; i < x.length; i++) assert.ok(Math.abs(r.x[i] - x[i] - .3 * d[i]) < 1e-15);
  const y = Float64Array.of(0, .1, 1, 1.5, 2);
  const h = proposeCoupledDensityNewton(system, y, Float64Array.of(0, 0, .8, -.1, 0));
  assert.equal(h.limiter.kind, 'bl-shape'); assert.ok(Math.abs(h.x[3] - h.x[2] - .125) < 1e-14);
  // Giles's H guard is conditional: do not clip a still-positive gap to
  // an unconditional 25% floor that is absent from the source listing.
  const near = proposeCoupledDensityNewton(system, y, Float64Array.of(0, 0, .4, -.04, 0));
  assert.equal(near.step, 1); assert.ok(Math.abs(near.x[3] - near.x[2] - .06) < 1e-14);
  const density = proposeCoupledDensityNewton(system, x, Float64Array.of(-2, .1, 0, 0, 0));
  assert.equal(density.step, .25); assert.ok(Math.abs(Math.exp(density.x[0]) - .6) < 1e-14);
  assert.throws(() => proposeCoupledDensityNewton(system, x, d, { maximumStep: 2 }), /maximum/);
});

test('a concave but volume-admissible coupled endpoint cannot pass the final gate', () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-ises-wake-fold.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
  assert.equal(system.admissible(system.initial), false);
  assert.equal(system.admissible(system.initial, { requireConvex: false }), true);
  const v = system.evaluate(system.initial);
  for (const key in f.families) assert.ok(Math.abs(v.families[key] - f.families[key]) < 1e-10);
  assert.equal(directStreamtubeVolumeGeometry(v.outer.nodes).valid, true);
  // Deliberately make the residual gate pass to isolate the final quality
  // gate. No production tolerance is changed by this regression.
  const r = coupledStreamtubeResult(system, system.initial, { tolerance: 1 });
  assert.ok(Math.max(...r.residual.map(Math.abs)) < 1); assert.equal(r.converged, false);
  assert.equal(r.mesh.quality.valid, false); assert.equal(r.mesh.initialization.flowSolved, false);
  assert.equal(r.status, 'unconverged'); assert.equal(r.cl, null);
  assert.throws(() => solveCoupledStreamtubeBody(system, { maxIterations: 0, tolerance: 1 }), /initial state/);
});

test('a rejected full coupled ISES proposal retains the last published state and material-trip phase', () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/default-newton-trip-crossing.json', import.meta.url)));
  const options = { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL), iterationGeometry: 'ises-sampled' };
  const initial = solveCoupledStreamtubeIses(f.input, { ...options, maxIterations: 0 });
  const snapshots = [], r = solveCoupledStreamtubeIses(f.input, { ...options, maxIterations: 1,
    stepAcceptance: 'listing', onMesh: m => snapshots.push(m.nodes) });
  assert.ok(r.lastRejectedStep); assert.equal(r.history.length, 1); assert.equal(snapshots.length, 1);
  assert.equal(r.linearDiagnostics.solves, 1); assert.deepEqual(r.x, initial.x); assert.deepEqual(r.flow.nodes, initial.flow.nodes);
  assert.deepEqual(r.families, initial.families); assert.deepEqual(snapshots[0], r.flow.nodes);
  const ne = r.x.length - 4 * r.boundaryLayer.stations.length;
  const system = createCoupledStreamtubeBody(r.solverInput, { ...r.coupledOptions, initialEuler: r.flow, initialBL: r.x.slice(ne) });
  assert.deepEqual(system.bl.activeTargets(system.initial.subarray(0, ne)).filter(t => t.from !== t.to), []);
  const v = system.evaluate(system.initial);
  for (const key in r.families) assert.ok(Math.abs(v.families[key] - r.families[key]) < 1e-10);
});

test('coupled ISES maintenance closes all four BLs and both wakes with conservation and restart Jacobian checks', t => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), before = structuredClone(input), snapshots = [];
  const r = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', maxIterations: 12, tolerance: 1e-10,
    stepAcceptance: 'admissible', onMesh: m => snapshots.push({ nodes: m.nodes, iteration: m.iteration.iteration }) });
  assert.deepEqual(input, before); assert.equal(r.converged, true, r.reason);
  assert.equal(r.initialRedistribution.accepted, true); assert.ok(r.initialRedistribution.passages.every(p => p.pairs === 5));
  assert.equal(r.boundaryLayer.surfaces.length, 4); assert.equal(r.boundaryLayer.wakes.length, 2);
  assert.equal(r.linearDiagnostics.solves, r.history.length - 1);
  assert.ok(r.linearDiagnostics.maxRelativeResidual < 1e-10);
  assert.deepEqual(snapshots.map(m => m.iteration), r.history.map(h => h.iteration));
  assert.deepEqual(snapshots.at(-1).nodes, r.flow.nodes);
  assert.equal(r.mesh.quality.valid, true); assert.equal(r.mesh.initialization.flowSolved, true);
  for (let g = 0; g < r.flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.flow.nodes[g], sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(Math.max(...c[key].map(Math.abs)) < 2e-9);
  }
  const ne = r.x.length - 4 * r.boundaryLayer.stations.length;
  const system = createCoupledStreamtubeBody(r.solverInput, { ...r.coupledOptions, initialEuler: r.flow, initialBL: r.x.slice(ne) });
  const x = system.initial, v = system.evaluate(x);
  assert.ok(Math.max(...v.residual.map(Math.abs)) < 1e-10);
  for (let i = 0; i < x.length; i++) assert.ok(Math.abs(x[i] - r.x[i]) < 1e-12);
  assert.deepEqual(system.bl.activeTargets(x.subarray(0, ne)).filter(t => t.from !== t.to), []);
  const jacobian = system.jacobian(x), h = 2e-6, errors = [];
  for (const family of ['euler', 'bl', 'both']) {
    const d = x.map((q, i) => (family === 'euler' && i >= ne || family === 'bl' && i < ne) ? 0 : (i < ne ? .001 : Math.max(.01, Math.abs(q))) * Math.sin(i * .43 + .2));
    const exact = sparseProduct(jacobian, d), p = system.residual(x.map((q, i) => q + h * d[i])), m = system.residual(x.map((q, i) => q - h * d[i]));
    const error = Math.max(...exact.map((q, i) => { const fd = (p[i] - m[i]) / (2 * h); return Math.abs(fd - q) / Math.max(1, Math.abs(fd), Math.abs(q)); }));
    assert.ok(error < 5e-6, `${family} ${error}`); errors.push({ family, error });
  }
  system.evaluate(x);
  for (const w of r.boundaryLayer.wakes) for (const id of w.ids.slice(1)) {
    const s = r.boundaryLayer.stations[id], a = r.flow.nodes[w.body][s.i].at(-1), b = r.flow.nodes[w.body + 1][s.i][0];
    assert.ok(Math.abs(Math.hypot(a.x - b.x, a.y - b.y) - system.euler.conditions.lengthScale * s.deltaStar) < 2e-12);
  }
  t.diagnostic(JSON.stringify({ unknowns: r.x.length, iterations: r.history.length - 1, families: r.families, errors }));
});
