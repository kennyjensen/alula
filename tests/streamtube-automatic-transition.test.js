import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody, coupledStreamtubeTripEvents } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { redistributeCoupledSurfaceStations } from '../src/euler/tests/streamtube-coupled-redistribution.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { evaluateTransitionInterval } from '../src/viscous/transition-interval.js';
import { evaluateLeadingTransitionInterval } from '../src/viscous/leading-transition-interval.js';

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const max = values => Math.max(...Array.from(values, Math.abs));
const fixture = () => intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
const options = { edgeMatching: 'section-velocity', transitionMode: 'automatic', reynolds: 1e6, ncrit: 9 };
const make = () => createCoupledStreamtubeBody(fixture(), options);

test('terminal event closes TE junctions and preserves downstream wake shear and physical primitives', () => {
  const s = createCoupledStreamtubeBody(fixture(), { ...options, ncrit: 100 });
  const candidate = s.initial.slice(), bl = candidate.subarray(s.ne), euler = candidate.subarray(0, s.ne);
  for (const branch of s.bl.surfaces) for (const id of branch.ids.slice(1)) bl[4 * id] = .3;
  s.bl.restoreActive(s.bl.surfaces.map(() => 1));
  const targets = s.bl.activeTargets(euler, bl);
  assert.ok(targets.every((t, k) => t.kind === 'trailing-edge' && t.to === s.bl.surfaces[k].ids.length - 1));
  const before = candidate.slice(), beforeGeometry = s.euler.decode(euler).nodes;
  const event = s.bl.updateActive(bl, euler), value = s.evaluate(candidate);
  assert.equal(event.changed, true);
  for (let i = 0; i < candidate.length; i++) if (i < s.ne || (i - s.ne) % 4 !== 0)
    assert.equal(candidate[i], before[i]);
  assert.deepEqual(s.euler.decode(euler).nodes, beforeGeometry);
  for (const branch of s.bl.surfaces) {
    const id = branch.ids.at(-1), p = value.layers.states[id];
    assert.equal(bl[4 * id], s.bl.kernel.station(p, 'turbulent').transitionShear);
    assert.ok(Math.abs(value.layers.residual[4 * id]) < 1e-12);
    const change = event.changes.find(c => c.body === branch.body && c.side === branch.side);
    assert.equal(change.converted.at(-1).terminalShearInitialization, true);
    assert.equal(change.terminalWakeShearInitialization.method, 'theta-weighted-TE-shear');
    assert.equal(change.downstreamWakeShearPreserved, true);
    assert.equal(change.wakeShearTransportInitialization, undefined);
  }
  for (const wake of s.bl.wakes) {
    assert.ok(Math.abs(value.layers.residual[4 * wake.ids[0]]) < 1e-12);
    for (const id of wake.ids.slice(1)) assert.equal(bl[4 * id], before[s.ne + 4 * id]);
  }
  const stable = candidate.slice();
  assert.deepEqual(s.bl.updateActive(bl, euler), { changed: false, changes: [] });
  assert.deepEqual(candidate, stable);
});

test('terminal transfer touches only affected surfaces and junctions and rolls back failed preparation', () => {
  const prepare = () => {
    const s = createCoupledStreamtubeBody(fixture(), { ...options, ncrit: 100 });
    const x = s.initial.slice(), bl = x.subarray(s.ne), euler = x.subarray(0, s.ne);
    const phase = s.bl.surfaces.map(branch => branch.ids.length - 1); phase[0] = 1;
    s.bl.restoreActive(phase);
    for (const id of s.bl.surfaces[0].ids.slice(1)) bl[4 * id] = .3;
    return { s, x, bl, euler };
  };
  const f = prepare(), before = f.x.slice();
  const event = f.s.bl.updateActive(f.bl, f.euler);
  assert.equal(event.changes.length, 1);
  const affected = f.s.bl.surfaces[0].body;
  for (const station of f.s.bl.stations) if (station.body !== affected || station.kind === 'surface' && station.side !== 'upper')
    for (let k = 0; k < 4; k++) assert.equal(f.bl[4 * station.id + k], before[f.s.ne + 4 * station.id + k]);
  for (const wake of f.s.bl.wakes) for (const id of wake.ids.slice(1))
    assert.equal(f.bl[4 * id], before[f.s.ne + 4 * id]);

  const failed = prepare(), original = failed.x.slice(), phase = failed.s.bl.snapshotActive();
  const station = failed.s.bl.kernel.station;
  failed.s.bl.kernel.station = (state, regime) => {
    if (regime === 'turbulent') throw new Error('deliberate terminal shear preparation failure');
    return station(state, regime);
  };
  assert.throws(() => failed.s.bl.updateActive(failed.bl, failed.euler), /terminal shear preparation failure/);
  assert.deepEqual(failed.x, original); assert.deepEqual(failed.s.bl.snapshotActive(), phase);
});

test('a natural interval event preserves downstream surface and wake shear', () => {
  const s = make(), x = s.initial.slice(), bl = x.subarray(s.ne), euler = x.subarray(0, s.ne);
  const targets = s.bl.activeTargets(euler, bl), phase = targets.map(t => t.to);
  const branchIndex = targets.findIndex((t, k) => t.kind === 'natural' && t.to > 0 && t.to < s.bl.surfaces[k].ids.length - 1);
  assert.ok(branchIndex >= 0, 'The fixture must contain a resolved natural interval with a downstream turbulent tail.');
  const branch = s.bl.surfaces[branchIndex]; phase[branchIndex]++;
  s.bl.restoreActive(phase);
  for (const id of branch.ids.slice(targets[branchIndex].to)) bl[4 * id] = .2;
  const before = x.slice(), event = s.bl.updateActive(bl, euler), value = s.evaluate(x);
  assert.equal(event.changes.length, 1);
  assert.equal(event.changes[0].downstreamShearPreserved, true);
  assert.equal(event.changes[0].surfaceShearTransportInitialization, undefined);
  const wake = s.bl.wakes.find(w => w.body === branch.body);
  const conversion = event.changes[0].converted.find(c => c.shearSelection);
  assert.equal(conversion.shearSelection.method, 'mixed-and-following-shear-minimum');
  const alternative = x.slice();
  alternative[s.ne + 4 * conversion.id] = conversion.shearSelection.initial.aux;
  const squared = residual => residual.reduce((sum, v) => sum + v * v, 0);
  assert.ok(squared(value.residual) <= squared(s.evaluate(alternative).residual),
    'Converted shear must improve the complete residual over closing only its mixed shear row.');
  if (conversion.shearSelection.selected === 'minimum') {
    assert.equal(conversion.transitionShearInitialization, undefined);
    assert.equal(conversion.shearInitialization, undefined);
  }
  for (const id of [...branch.ids.slice(Math.max(phase[branchIndex], targets[branchIndex].to) + 1), ...wake.ids])
    assert.equal(bl[4 * id], before[s.ne + 4 * id], `Existing downstream shear at BL station ${id}`);
  for (let i = 0; i < x.length; i++) if (i < s.ne || (i - s.ne) % 4 !== 0) assert.equal(x[i], before[i]);
});

test('downstream natural-event preparation selects the lower complete coupled residual', () => {
  const s = make(), x = s.initial.slice(), bl = x.subarray(s.ne), euler = x.subarray(0, s.ne);
  const targets = s.bl.activeTargets(euler, bl), phase = targets.map(t => t.to);
  const k = targets.findIndex((t, i) => t.kind === 'natural' && t.to > 1 && t.to + 1 < s.bl.surfaces[i].ids.length);
  assert.ok(k >= 0);
  const branch = s.bl.surfaces[k]; phase[k]--;
  s.bl.restoreActive(phase);
  for (const id of branch.ids.slice(phase[k])) bl[4 * id] = .2;
  const before = x.slice(), event = s.bl.updateActive(bl, euler);
  const conversion = event.changes.find(c => c.body === branch.body && c.side === branch.side)
    .converted.find(c => c.shearSelection);
  assert.ok(conversion, 'Coupled preparation must compare the existing and initialized shear.');
  const selection = conversion.shearSelection;
  const alternative = x.slice();
  alternative[s.ne + 4 * conversion.id] = selection[selection.selected === 'existing' ? 'initialized' : 'existing'].aux;
  const squared = state => s.evaluate(state).residual.reduce((sum, v) => sum + v * v, 0);
  assert.ok(squared(x) <= squared(alternative), 'Selection must reduce the complete residual, including downstream transport.');
  for (let i = 0; i < x.length; i++) if (i < s.ne || (i - s.ne) % 4 !== 0) assert.equal(x[i], before[i]);
});

test('fixed-trip controlled state, residual and complete sparse Jacobian retain their captured values', () => {
  const b = JSON.parse(fs.readFileSync(new URL('../docs/automatic-transition-integration/before.json', import.meta.url)));
  const s = createCoupledStreamtubeBody(b.input, b.options), v = s.evaluate(s.initial), j = s.jacobian(s.initial);
  // JSON does not retain the sign of zero; compare the same serialization.
  assert.deepEqual(serialize(s.initial), b.state); assert.deepEqual(v.outer.nodes, b.nodes);
  assert.deepEqual(v.families, b.families); assert.deepEqual(serialize(v.residual), b.residual);
  for (const key of ['rowPtr', 'colIndex', 'values']) assert.deepEqual(serialize(j[key]), b.jacobian[key]);
});

test('all automatic two-element Jacobian columns include moving natural roots and both coupling directions', t => {
  const s = make(), x = s.initial, { n, ne } = s, j = s.jacobian(x, { sparse: false });
  const transitions = s.evaluate(x).layers.transitions;
  assert.equal(transitions.length, 4); assert.ok(transitions.every(r => r.kind === 'natural'));
  const errors = [0, 0, 0, 0], strength = [0, 0];
  for (let col = 0; col < n; col++) {
    const h = 2e-7 * Math.max(1, Math.abs(x[col])), p = x.slice(), m = x.slice(); p[col] += h; m[col] -= h;
    const a = s.residual(p), b = s.residual(m);
    for (let row = 0; row < n; row++) {
      const exact = j[row * n + col], fd = (a[row] - b[row]) / (2 * h);
      const block = 2 * Number(row >= ne) + Number(col >= ne);
      errors[block] = Math.max(errors[block], Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd)));
      if (row < ne && col >= ne) strength[0] += Math.abs(exact);
      if (row >= ne && (row - ne) % 4 === 3 && col < ne) strength[1] += Math.abs(exact);
    }
  }
  assert.ok(errors.every(e => e < 5e-6), JSON.stringify(errors)); assert.ok(strength.every(v => v > 1e-8));
  t.diagnostic(JSON.stringify({ n, errors, strength }));
});

test('simultaneous two-element Newton closes four natural-transition BLs, both wakes and conservative flow', t => {
  const s = make(), before = s.evaluate(s.initial), r = solveCoupledStreamtubeBody(s, { maxIterations: 12, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason); assert.ok(max(r.residual) < 1e-10); assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.boundaryLayer.surfaces.length, 4); assert.equal(r.boundaryLayer.wakes.length, 2);
  assert.ok(r.boundaryLayer.transitions.every(q => q.kind === 'natural'));
  assert.ok(r.boundaryLayer.transitions.some((q, i) => Math.abs(q.s - before.layers.transitions[i].s) > 1e-5));
  for (let g = 0; g < r.flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.flow.nodes[g], sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(max(c[key]) < 2e-9, `${g}/${key}`);
  }
  for (const w of r.boundaryLayer.wakes) {
    const [a, b] = r.boundaryLayer.surfaces.filter(q => q.body === w.body).map(q => r.boundaryLayer.stations[q.ids.at(-1)]);
    const z = r.boundaryLayer.stations[w.ids[0]];
    assert.ok(Math.abs(z.theta - a.theta - b.theta) < 1e-12);
    assert.ok(Math.abs(z.deltaStar - a.deltaStar - b.deltaStar) < 1e-12);
    assert.ok(Math.abs(z.aux * z.theta - a.aux * a.theta - b.aux * b.theta) < 1e-12);
  }
  const restart = { ...options, initialEuler: r.flow, initialBL: r.x.slice(s.ne), transitionState: r.boundaryLayer.transitionState };
  const again = createCoupledStreamtubeBody(fixture(), restart);
  assert.ok(max(again.residual(again.initial)) < 1e-10); assert.deepEqual(again.bl.snapshotActive(), r.boundaryLayer.transitionState);
  const missing = { ...restart }; delete missing.transitionState;
  assert.throws(() => createCoupledStreamtubeBody(fixture(), missing), /requires its transition-interval map/);
  t.diagnostic(JSON.stringify({ iterations: r.history.length - 1, families: r.families, transitions: r.boundaryLayer.transitions }));
});

test('automatic ISES accepted checkpoints preserve transition phases through maintenance and resume', t => {
  const controls = { ...options, tolerance: 1e-10, maxIterations: 16, stepAcceptance: 'admissible' };
  const full = solveCoupledStreamtubeIses(fixture(), controls);
  assert.equal(full.converged, true, full.reason);
  let saved;
  assert.throws(() => solveCoupledStreamtubeIses(fixture(), { ...controls, onCheckpoint: (c, details) => {
    if (details.history.at(-1).iteration === 3) { saved = serialize(c); throw new Error('saved automatic interruption'); }
  } }), /saved automatic interruption/);
  assert.equal(saved.restart.options.transitionMode, 'automatic'); assert.equal(saved.restart.options.transitionState.length, 4);
  const before = structuredClone(saved), zero = solveCoupledStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 0 });
  assert.deepEqual(zero.families, saved.families); assert.deepEqual(zero.boundaryLayer.transitionState, saved.restart.options.transitionState);
  assert.equal(zero.initialRedistribution.resumed, true);
  const resumed = solveCoupledStreamtubeIses(undefined, { ...controls, resume: saved, maxIterations: 13 });
  assert.equal(resumed.converged, true, resumed.reason); assert.deepEqual(saved, before);
  assert.ok(max(resumed.x.map((v, i) => v - full.x[i])) < 1e-11);
  assert.equal(resumed.boundaryLayer.wakes.length, 2); assert.deepEqual(resumed.boundaryLayer.transitionState, full.boundaryLayer.transitionState);
  t.diagnostic(JSON.stringify({ iterations: full.history.length - 1, families: full.families }));
});

test('normal refinement and identity surface transfer retain automatic state meanings', () => {
  const input = { ...fixture(), streamwiseMode: 'isentropic' }, s = createCoupledStreamtubeBody(input, options);
  const identity = redistributeCoupledSurfaceStations(input, s, { surfaceFractions: s.euler.fractions });
  assert.equal(identity.options.transitionMode, 'automatic'); assert.deepEqual(identity.system.bl.snapshotActive(), s.bl.snapshotActive());
  const refined = refineCoupledStreamtubeBody(input, s, { streamwiseFactor: 1, normalFactor: 2 });
  assert.equal(refined.options.transitionMode, 'automatic'); assert.deepEqual(refined.options.transitionState, s.bl.snapshotActive());
  assert.equal(refined.system.bl.surfaces.length, 4); assert.equal(refined.system.bl.wakes.length, 2);
  assert.deepEqual(Array.from(refined.initialBL), Array.from(s.initial.slice(s.ne)));
});

test('a natural interval crossing near a coupled root commits only auxiliary conversions and Newton crosses back', t => {
  const s = make(), root = solveCoupledStreamtubeBody(s, { maxIterations: 12, tolerance: 1e-10 });
  assert.equal(root.converged, true);
  const before = root.x.slice(), events = coupledStreamtubeTripEvents(s), phase = events.snapshot();
  const baseline = serialize(s.residual(before)), candidate = before.slice(), id = s.bl.surfaces[0].ids[1];
  // A one-percent local displacement change moves this retained natural root
  // just into the next interval. This changes a solved flow variable, not Ncrit.
  candidate[s.ne + 4 * id + 2] *= .99;
  const physical = candidate.slice(), event = events.prepare(candidate, before);
  assert.equal(event.changed, true); assert.equal(event.changes.length, 1);
  assert.deepEqual([event.changes[0].from, event.changes[0].to, event.changes[0].kind], [1, 2, 'natural']);
  for (let i = 0; i < candidate.length; i++) if (i < s.ne || (i - s.ne) % 4 !== 0) assert.equal(candidate[i], physical[i]);
  assert.equal(s.admissible(candidate), true);
  events.restore(phase); assert.deepEqual(serialize(s.residual(before)), baseline);
  const moved = physical.slice(); events.prepare(moved, before);
  const r = solveCoupledStreamtubeBody(s, { initial: moved, maxIterations: 16, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason); assert.deepEqual(r.boundaryLayer.transitionState, phase);
  assert.ok(r.history.some(h => h.activeChange && h.changes.some(c => c.from === 2 && c.to === 1)));
  t.diagnostic(JSON.stringify({ iterations: r.history.length - 1, events: r.history.filter(h => h.activeChange), families: r.families }));
});

test('trips before the first resolved nodes couple all four leading intervals and their grid derivatives', () => {
  const s = createCoupledStreamtubeBody(fixture(), { ...options, tripFractions: [[.1, .1], [.15, .15]] });
  const x = s.initial, j = s.jacobian(x, { sparse: false });
  const columns = [...s.euler.layout.globals.stagnation, ...s.bl.surfaces.flatMap(b => [0, 1, 2, 3].map(k => s.ne + 4 * b.ids[0] + k))];
  for (const col of columns) {
    const h = 2e-7 * Math.max(1, Math.abs(x[col])), p = x.slice(), m = x.slice(); p[col] += h; m[col] -= h;
    const a = s.residual(p), b = s.residual(m);
    for (let row = 0; row < s.n; row++) {
      const exact = j[row * s.n + col], fd = (a[row] - b[row]) / (2 * h);
      assert.ok(Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd)) < 5e-6, `${row}/${col}: ${exact}/${fd}`);
    }
  }
  const r = solveCoupledStreamtubeBody(s, { maxIterations: 12, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason); assert.ok(r.boundaryLayer.transitionState.every(j => j === 0));
  assert.ok(r.boundaryLayer.transitions.every(q => q.kind === 'forced'));
});

test('a mixed natural/laminar-to-TE two-element solve closes both viscous wakes', t => {
  const s = createCoupledStreamtubeBody(fixture(), { ...options, reynolds: 3e5, ncrit: 14 });
  const r = solveCoupledStreamtubeBody(s, { maxIterations: 12, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason); assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.boundaryLayer.transitions.filter(q => q.kind === 'trailing-edge').length, 3);
  for (const q of r.boundaryLayer.transitions.filter(q => q.kind === 'trailing-edge')) {
    const surface = s.bl.surfaces.find(b => b.body === q.body && b.side === q.side);
    assert.equal(q.id, surface.ids.at(-1)); assert.ok(q.amplification < 14);
    assert.ok(surface.ids.slice(0, -1).every(id => ['similarity', 'laminar'].includes(s.bl.stations[id].regime)));
    assert.ok(r.boundaryLayer.stations[q.id].aux > 0 && r.boundaryLayer.stations[q.id].aux < .3);
  }
  assert.equal(r.boundaryLayer.wakes.length, 2);
  t.diagnostic(JSON.stringify({ iterations: r.history.length - 1, families: r.families, transitions: r.boundaryLayer.transitions }));
});

test('terminal transition keeps laminar thickness equations and supplies shear, never N, to a wake', () => {
  const k = createIntegralKernel({ reynolds: 1e5, mach: .2, ncrit: 9, exactJacobian: true });
  const upstream = { s: .8, ue: 1, theta: .002, deltaStar: .005, aux: .7 };
  const downstream = { s: 1, ue: .98, theta: .0022, deltaStar: .0057, aux: .03 };
  downstream.aux = k.station(downstream, 'turbulent').transitionShear;
  const r = evaluateTransitionInterval(k, { upstream, downstream, tripS: downstream.s }, { jacobian: true });
  assert.equal(r.transition.forced, true); assert.equal(r.transition.s, downstream.s);
  const n = k.transitionCheck({ upstream, downstream }).amplification;
  const laminar = k.interval({ upstream, downstream: { ...downstream, aux: n }, regime: 'laminar' });
  for (const row of [1, 2]) assert.ok(Math.abs(r.residual[row] - laminar.residual[row]) < 1e-12);
  assert.ok(Math.abs(r.residual[0]) < 1e-12); assert.notEqual(downstream.aux, n);
});

test('terminal transition agrees with executed original Fortran at the exact endpoint', () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/terminal-transition.json', import.meta.url)));
  for (const [path, hash] of Object.entries(f.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${path}`, import.meta.url))).digest('hex'), hash, path);
  for (const c of f.cases) {
    const k = createIntegralKernel({ ...c.parameters, exactJacobian: true }), r = evaluateTransitionInterval(k, c);
    assert.equal(r.transition.s, c.tripS); assert.equal(r.transition.forced, true);
    for (let row = 0; row < 3; row++) assert.ok(Math.abs(r.residual[row] - c.expected.transition[row]) < 3e-10);
    assert.ok(Math.abs(k.station(c.downstream, 'turbulent').transitionShear - c.expected.shear) < 3e-10);
    for (const row of [1, 2]) assert.ok(Math.abs(r.residual[row] - c.expected.laminar[row]) < 3e-10);
  }
});
test('leading transition derivatives include the virtual similarity state and its moving trip', () => {
  const k = createIntegralKernel({ reynolds: 1e6, mach: .2, ncrit: 9, exactJacobian: true });
  const input = { downstream: { s: .02, ue: .4, theta: .0001, deltaStar: .00022, aux: .02 }, tripS: .0001 };
  const r = evaluateLeadingTransitionInterval(k, input, { jacobian: true });
  assert.equal(r.transition.forced, true); assert.equal(r.transition.s, input.tripS);
  assert.ok(Math.abs(r.upstreamState.ue / r.upstreamState.s - input.downstream.ue / input.downstream.s) < 1e-12);
  for (const [col, key] of ['aux', 'theta', 'deltaStar', 'ue', 's', 'tripS'].entries()) {
    const v = col === 5 ? input.tripS : input.downstream[key], h = 2e-5 * v;
    const values = [-2, -1, 1, 2].map(offset => evaluateLeadingTransitionInterval(k, col === 5
      ? { ...input, tripS: v + offset * h } : { ...input, downstream: { ...input.downstream, [key]: v + offset * h } }).residual);
    for (let row = 0; row < 3; row++) {
      const fd = (values[0][row] - 8 * values[1][row] + 8 * values[2][row] - values[3][row]) / (12 * h);
      const exact = col === 5 ? r.partials.trip[row] : r.partials.downstream[row][col];
      assert.ok(Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd)) < 5e-6, `${key}/${row}: ${exact}/${fd}`);
    }
  }
});
