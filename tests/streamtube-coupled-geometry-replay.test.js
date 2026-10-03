// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { pathToFileURL } from 'node:url';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { transferStreamtubeGeometry } from '../src/euler/streamtube-geometry.js';
import { streamtubeMotionDirections } from '../src/euler/streamtube-geometry.js';
import { sparseProduct } from '../src/numerics/sparse.js';

// A private source copy can qualify this change before touching a running
// browser's numerical closure. Normal test runs exercise the live modules.
const root = process.env.MSES_GEOMETRY_REPLAY_RUNTIME;
const moduleURL = name => root ? pathToFileURL(`${root}/src/euler/${name}.js`)
  : new URL(`../src/euler/${name}.js`, import.meta.url);
const { createStreamtubeBodySystem } = await import(moduleURL('streamtube-body'));
const { createCoupledStreamtubeBody } = await import(moduleURL('streamtube-coupled'));
const { solveCoupledStreamtubeIses } = await import(moduleURL('streamtube-coupled-ises'));
const json = a => JSON.parse(JSON.stringify(a, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));

function fixture() {
  const f = twoActiveFiniteBaseWakes(), input = { ...f.input, wakeDisplacementMotion: 'te-center',
    stagnationMotion: 'walls-only', normalStencil: 'body-stations' };
  const euler = createStreamtubeBodySystem({ ...input, displacement: f.system.euler.displacement });
  const state = transferStreamtubeGeometry(f.system.euler, f.x.subarray(0, f.system.ne), euler);
  const decoded = euler.decode(state);
  const options = { reynolds: 1e6, ncrit: 9, edgeMatching: 'section-velocity',
    geometryReplay: 'preserve-undisplaced' };
  const initialEuler = { x: state, nodes: decoded.nodes, undisplacedNodes: decoded.undisplacedNodes };
  const initialBL = f.x.slice(f.system.ne);
  const system = createCoupledStreamtubeBody(input, { ...options, initialEuler, initialBL });
  return { input, options, initialEuler, initialBL, system };
}

test('physical NCALC direction refresh at a chart origin preserves every coordinate and rejects nonzero positions atomically', () => {
  const f = fixture(), s = f.system.euler, x = f.system.initial.slice(0, f.system.ne);
  const displacement = structuredClone(s.displacement);
  displacement.surfaces[0].upper[displacement.surfaces[0].upper.length - 1] *= 1.4;
  s.setDisplacement(displacement);
  const before = s.decode(x), packed = x.slice(), chart = s.geometryChart();
  const normals = streamtubeMotionDirections(s.layout, before.nodes, 'body-stations');
  s.refreshGeometryDirections(x);
  assert.deepEqual(x, packed); assert.deepEqual(s.decode(x), before);
  for (const p of s.geometryChart()) {
    assert.deepEqual(p.normal, normals.get(p.column));
    assert.deepEqual(p.offset, chart.find(q => q.column === p.column).offset);
  }
  const moved = x.slice(); moved[s.layout.positions[0].column] = Number.MIN_VALUE;
  const retained = s.geometryChart();
  assert.throws(() => s.refreshGeometryDirections(moved), /zero free-position/);
  assert.deepEqual(s.geometryChart(), retained); assert.deepEqual(s.decode(x), before);
});

test('the explicit raw-chart replay preserves full state, physical geometry, residual and phase across JSON restarts', () => {
  const f = fixture(); let system = f.system;
  const original = system.evaluate(system.initial), packed = system.initial.slice();
  for (let k = 0; k < 3; k++) {
    const value = system.evaluate(system.initial);
    const restart = json({ input: f.input, options: f.options,
      initialEuler: { x: system.initial.slice(0, system.ne), nodes: value.outer.nodes,
        undisplacedNodes: value.outer.undisplacedNodes }, initialBL: system.initial.slice(system.ne) });
    system = createCoupledStreamtubeBody(restart.input, { ...restart.options,
      initialEuler: restart.initialEuler, initialBL: restart.initialBL });
    const replay = system.evaluate(system.initial);
    assert.deepEqual(system.initial, packed);
    assert.deepEqual(replay.outer.nodes, original.outer.nodes);
    assert.deepEqual(replay.outer.undisplacedNodes, original.outer.undisplacedNodes);
    assert.deepEqual(replay.residual, original.residual);
    assert.deepEqual(replay.outer.allocation, original.outer.allocation);
    assert.deepEqual(system.bl.snapshotActive(), f.system.bl.snapshotActive());
    assert.equal(system.conditions.geometryReplay, 'preserve-undisplaced');
  }
  for (const geometryReplay of [null, true, 'unknown']) assert.throws(() =>
    createCoupledStreamtubeBody(f.input, { ...f.options, geometryReplay, initialEuler: f.initialEuler,
      initialBL: f.initialBL }), /geometry replay policy/);
  assert.throws(() => createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: {
    x: f.initialEuler.x, nodes: f.initialEuler.nodes }, initialBL: f.initialBL }), /supplied undisplaced/);
  assert.throws(() => createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler }), /supplied undisplaced/);
});

test('raw-chart replay retains full coupled state and TE-displacement derivatives', () => {
  const f = fixture(), s = f.system, x = s.initial, matrix = s.jacobian(x), h = 1e-6;
  for (const group of ['Euler', 'BL', 'TE']) {
    const d = x.map((v, i) => {
      if (group === 'Euler') return i < s.ne ? .001 * Math.sin(i * .47 + .2) : 0;
      if (group === 'BL') return i >= s.ne ? Math.max(.001, Math.abs(v)) * .01 * Math.sin(i * .31 + .4) : 0;
      const station = s.bl.stations[Math.floor((i - s.ne) / 4)];
      return station?.kind === 'surface' && (i - s.ne) % 4 === 2
        && station.i === s.euler.layout.bodies[station.body].trailingIndex ? .001 : 0;
    });
    const derivative = sparseProduct(matrix, d);
    const plus = s.residual(x.map((v, i) => v + h * d[i]));
    const minus = s.residual(x.map((v, i) => v - h * d[i]));
    const error = Math.max(...derivative.map((v, i) => {
      const finite = (plus[i] - minus[i]) / (2 * h);
      return Math.abs(v - finite) / Math.max(1, Math.abs(v), Math.abs(finite));
    }));
    assert.ok(error < 5e-6, `${group} directional derivative ${error}`);
  }
});

test('ISES persists the explicit replay policy and refuses incompatible resumes', () => {
  const f = fixture(), controls = { maxIterations: 0, tolerance: 1e-10,
    iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing' };
  const result = solveCoupledStreamtubeIses(f.input, { ...f.options, initialEuler: f.initialEuler,
    initialBL: f.initialBL, ...controls });
  assert.ok(result.checkpoint, result.reason);
  assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(result.checkpoint.restart.options.geometryReplay, 'preserve-undisplaced');
  const checkpoint = json(result.checkpoint);
  const again = solveCoupledStreamtubeIses(undefined, { resume: checkpoint, ...controls });
  assert.deepEqual(again.x, result.x); assert.deepEqual(again.flow.nodes, result.flow.nodes);
  assert.deepEqual(again.residual, result.residual);
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { resume: checkpoint,
    ...controls, geometryReplay: 'legacy' }), /geometry replay controls/);
  delete checkpoint.restart.options.geometryReplay;
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { resume: checkpoint,
    ...controls, geometryReplay: 'preserve-undisplaced' }), /geometry replay controls/);
});
