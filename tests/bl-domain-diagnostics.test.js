import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { selectSurfaceTransition } from '../src/viscous/transition-selection.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createStreamtubeBoundaryLayers } from '../src/euler/streamtube-boundary-layers.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const root = new URL('../', import.meta.url);
const before = JSON.parse(fs.readFileSync(new URL('docs/bl-domain-diagnostics/before.json', root)));
async function archivedModule(path) {
  const entry = before.files.find(row => row.path === path);
  const bytes = fs.readFileSync(new URL(entry.archive, root));
  assert.equal(createHash('sha256').update(bytes).digest('hex'), entry.sha256);
  const source = bytes.toString().replace(/from '(\.[^']+)'/g,
    (_, name) => `from '${new URL(name, new URL(path, root)).href}'`);
  return import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
}
const archived = await archivedModule('src/viscous/integral.js');
const archivedSelection = await archivedModule('src/viscous/transition-selection.js');
const state = { s: .2, theta: .0002, deltaStar: .0006, ue: 1, aux: 0 };
const failure = fn => {
  let caught;
  assert.throws(fn, error => { caught = error; return true; });
  return caught;
};

test('raw Hk rejection identifies the shape domain and retains exact gap subtraction and bound', () => {
  const kernel = createIntegralKernel({ mach: .4 });
  for (const wakeGap of [0, .03]) {
    const input = { ...state, wakeGap, deltaStar: wakeGap + 1.01 * state.theta };
    const saved = structuredClone(input), error = failure(() => kernel.station(input, wakeGap ? 'wake' : 'laminar'));
    assert.match(error.message, /^Inadmissible compressible BL edge state\. raw kinematic shape factor Hk must exceed 1/);
    assert.equal(error.code, 'BL_EDGE_STATE_DOMAIN');
    assert.equal(error.diagnostics.condition, 'raw-hk');
    assert.deepEqual(error.diagnostics.failedChecks, ['raw-hk']);
    const h = (input.deltaStar - wakeGap) / input.theta, m2 = .4 ** 2;
    assert.ok(Math.abs(error.diagnostics.rawHk - (h - .29 * m2) / (1 + .113 * m2)) < 2e-15);
    assert.ok(error.diagnostics.enthalpyRatio > 0 && error.diagnostics.density > 0 && error.diagnostics.viscosity > 0);
    assert.deepEqual(error.station, input); assert.deepEqual(input, saved);
    assert.throws(() => archived.createIntegralKernel({ mach: .4 }).station(input), /^Error: Inadmissible compressible BL edge state\./);
  }
  const zeroMach = createIntegralKernel();
  assert.throws(() => zeroMach.station({ ...state, deltaStar: state.theta }), /Inadmissible integral BL station/);
  assert.ok(zeroMach.station({ ...state, deltaStar: state.theta * (1 + 1e-10) }).rawHk > 1);
});

test('thermal, viscosity and Reynolds-number failures distinguish their actual failed quantities', () => {
  for (const [condition, parameters, input, text] of [
    ['thermal-energy', { mach: .4 }, { ...state, ue: 8 }, /thermal-energy ratio/],
    ['viscosity', { reynolds: Number.MIN_VALUE }, state, /edge viscosity/],
    ['re-theta', { reynolds: 1e308 }, { ...state, theta: 1e10, deltaStar: 3e10 }, /Reynolds number/]
  ]) {
    const error = failure(() => createIntegralKernel(parameters).station(input));
    assert.equal(error.diagnostics.condition, condition); assert.match(error.message, text);
    assert.throws(() => archived.createIntegralKernel(parameters).station(input), /Inadmissible compressible BL edge state/);
  }
});

test('all accepted native fixture blocks and derivatives remain exactly equal to the archived adapter', () => {
  const fixture = JSON.parse(fs.readFileSync(new URL('tests/fixtures/fortran/integral.json', root)));
  for (const { input: c } of fixture.cases) for (const exactJacobian of [false, true]) {
    const current = createIntegralKernel({ ...c.parameters, exactJacobian });
    const old = archived.createIntegralKernel({ ...c.parameters, exactJacobian });
    const evaluate = kernel => {
      if (c.regime !== 'te') return kernel.interval({ upstream: c.states[0], downstream: c.states[1], regime: c.regime, tripS: c.tripS });
      const half = { ...c.states[0], ...c.matched, theta: c.matched.theta / 2, deltaStar: c.matched.deltaStar / 2 };
      const [upper, lower] = c.surfaceStates ?? [half, half];
      return kernel.trailingEdge(upper, lower, c.states[1], c.gap ?? 0);
    };
    assert.deepEqual(evaluate(current), evaluate(old));
  }
});

test('selector validates a hidden tail with its exact index and preserves accepted selection output', () => {
  const states = [state, { ...state, s: .3 }, { ...state, s: .4, deltaStar: state.theta * 1.01 }];
  const error = failure(() => selectSurfaceTransition(createIntegralKernel({ mach: .4 }), states, { tripS: .1 }));
  assert.equal(error.stationIndex, 2); assert.equal(error.diagnostics.surfaceStationIndex, 2);
  assert.equal(error.diagnostics.condition, 'raw-hk');
  const valid = states.map(s => ({ ...s, deltaStar: state.deltaStar }));
  assert.deepEqual(selectSurfaceTransition(createIntegralKernel({ mach: .4 }), valid, { tripS: .1 }),
    archivedSelection.selectSurfaceTransition(archived.createIntegralKernel({ mach: .4 }), valid, { tripS: .1 }));
});

test('Euler automatic selection adds original element, side and station without changing phases or packed state', () => {
  // A prescribed tiny grid and packed BL data: no mesh generator, BL initializer or Newton solve.
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 1, contourPanels: 40, mach: .4 });
  input.bodies[0].element = 6;
  const euler = createStreamtubeBodySystem(input), bl = createStreamtubeBoundaryLayers(euler, euler.initial, { transitionMode: 'automatic' });
  const x = new Float64Array(bl.stations.length * 4);
  for (const station of bl.stations) x.set([0, state.theta / bl.scale, state.deltaStar / bl.scale, 1], 4 * station.id);
  const surface = bl.surfaces[0], id = surface.ids.at(-1);
  x[4 * id + 2] = state.theta * 1.01 / bl.scale;
  const saved = x.slice(), phases = bl.snapshotActive();
  const error = failure(() => bl.activeTargets(euler.initial, x));
  assert.match(error.message, /Element 7, upper, station \d+ \(BL id \d+, grid station \d+\)/);
  assert.deepEqual(error.diagnostics.stationLocation, { element: 6, body: 0, side: 'upper',
    surfaceStationIndex: surface.ids.length - 1, id, gridStation: bl.stations[id].i,
    regime: bl.stations[id].regime, indexBase: 0 });
  assert.equal(error.diagnostics.condition, 'raw-hk'); assert.deepEqual(x, saved); assert.deepEqual(bl.snapshotActive(), phases);
  // This is the exact message field forwarded by Worker and retained in ISES rejection reasons.
  assert.match(JSON.parse(JSON.stringify({ message: error.message, diagnostics: error.diagnostics })).message, /Hk must exceed 1.*Element 7/);
});
