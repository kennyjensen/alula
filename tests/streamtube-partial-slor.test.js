// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { relaxStreamtubeInitialGrid } from '../src/euler/streamtube-elliptic-initializer.js';
import { prepareStreamtubeMesh } from '../src/euler/streamtube-mesh-preview.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { requireStreamtubeInitialGridDomain } from '../src/euler/streamtube-grid-repair.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { slorNumericalError } from '../src/geometry/slor-termination.js';
import { smoothEllipticStreamtubeGrid, createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
import { smoothPairedBoundaryGrid } from '../src/geometry/paired-boundary-slor.js';
const xi = [0, .5, 1], eta = [0, .25, .65, 1], massFlows = [.25, .4, .35];
const fixture = curvature => [-.1, 0, .1].map(x => eta.map(y => ({ x, y: .5 * curvature * x * x + y })));
const prepared = nodes => ({ nodes: [nodes], initial: [], system: { layout: { primaryBody: 0, bodies: [] },
  decode: () => ({ allocation: { groups: [massFlows.map(massFlow => ({ massFlow }))] } }) } });
const controls = { seed: 'supplied', boundaryControl: 'wall-angle', streamwiseCoordinates: xi, maxSweeps: 200, tolerance: 1e-9 };

async function beforeInitializer() {
  const before = 'docs/solver-reliability/gui-defaults-current/30p30n/partial-slor-promotion/before';
  const cache = new Map();
  function moduleURL(relative) {
    if (cache.has(relative)) return cache.get(relative);
    const text = fs.readFileSync(path.join(before, relative), 'utf8').replace(/from (['"])(\.\.?\/[^'"]+)\1/g, (_, quote, specifier) => {
      const dependency = path.normalize(path.join(path.dirname(relative), specifier));
      const url = fs.existsSync(path.join(before, dependency)) ? moduleURL(dependency) : pathToFileURL(path.resolve(dependency)).href;
      return `from '${url}'`;
    });
    const url = `data:text/javascript;base64,${Buffer.from(text).toString('base64')}`;
    cache.set(relative, url); return url;
  }
  return (await import(moduleURL('src/euler/streamtube-elliptic-initializer.js'))).relaxStreamtubeInitialGrid;
}

test('successful angle and harmonic paths retain exact pre-promotion nodes and reports apart from WASM hypot roundoff', async () => {
  const old = await beforeInitializer();
  for (const curvature of [0, .8, 10]) {
    const input = prepared(fixture(curvature));
    const expected = old(input, controls), actual = relaxStreamtubeInitialGrid(input, { ...controls, allowPartialInitialGuess: true });
    assert.ok(expected.converged && actual.converged);
    assert.deepEqual(actual.nodes, expected.nodes);
    // AS and V8 hypot differ by a few ulps in quality diagnostics. All other
    // numeric fields, iteration counts, decisions and coordinates stay exact.
    const qualityRoundoff = value => JSON.parse(JSON.stringify(value, (key, v) =>
      key === 'minCornerSine' && Number.isFinite(v) ? Number(v.toFixed(13)) : v));
    assert.deepEqual(qualityRoundoff(actual), qualityRoundoff(expected));
  }
});

test('sweep exhaustion exposes a separate healthy initial guess without claiming the requested tolerance', () => {
  const input = prepared(fixture(.8));
  const strict = relaxStreamtubeInitialGrid(input, { ...controls, maxSweeps: 0 });
  const partial = relaxStreamtubeInitialGrid(input, { ...controls, maxSweeps: 0, allowPartialInitialGuess: true });
  assert.equal(strict.converged, false); assert.equal(strict.partialInitialGuess, undefined);
  assert.equal(partial.converged, false); assert.deepEqual(partial.regions, strict.regions);
  assert.deepEqual(partial.nodes, strict.nodes); assert.deepEqual(partial.partialInitialGuess.nodes, input.nodes);
  const selected = partial.partialInitialGuess.reports[0];
  assert.equal(selected.iteration, 0); assert.ok(selected.residual > controls.tolerance);
  assert.equal(selected.requestedFailure.requestedTolerance, controls.tolerance);
});

test('observer errors cannot trigger fallback, even if they impersonate a numerical termination', () => {
  for (const curvature of [.8, 10]) {
    const input = prepared(fixture(curvature)), publications = [];
    assert.throws(() => relaxStreamtubeInitialGrid(input, { ...controls, allowPartialInitialGuess: true, onSweep: h => {
      publications.push(h); throw slorNumericalError('No admissible decreasing Newton step for eta rows 1, 2.', 'no-admissible-decreasing-line-step');
    } }), error => error.code === 'slor-observer-failed');
    assert.equal(publications.length, 1);
  }
});

test('progress snapshots remain isolated from smoothing and retained recovery seeds', () => {
  for (const curvature of [.8, 10]) for (const maxSweeps of [1, 200]) {
    const input = prepared(fixture(curvature)), before = structuredClone(input.nodes);
    const options = { ...controls, maxSweeps, allowPartialInitialGuess: true };
    const expected = relaxStreamtubeInitialGrid(input, options), snapshots = [];
    const actual = relaxStreamtubeInitialGrid(input, { ...options, onSweep: (history, nodes) => {
      snapshots.push({ history, nodes });
      history.residual = NaN; history.invalidCells = Infinity;
      nodes[0][0].x = NaN; nodes[1][1].y = Infinity;
    } });
    assert.ok(snapshots.length > 0);
    assert.deepEqual(actual, expected);
    assert.deepEqual(input.nodes, before);
    for (const snapshot of snapshots) snapshot.nodes.length = 0;
    assert.deepEqual(actual, expected, 'retained snapshots cannot alias returned or recovery nodes');
  }
});

test('an observer failure is terminal across all passages; nonfinite measures are not numerical exhaustion', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ bodySegments: 4, tubes: 3, contourPanels: 40, elements: 2 }));
  let calls = 0;
  assert.throws(() => relaxStreamtubeInitialGrid({ system }, { seed: 'supplied', boundaryControl: 'wall-angle',
    allowPartialInitialGuess: true, onSweep: () => { calls++; throw new Error('observer stopped'); } }), /observer stopped/);
  assert.equal(calls, 1);
  const nodes = system.decode(system.initial).nodes, phases = [];
  assert.throws(() => prepareStreamtubeMesh({ system, initial: system.initial, nodes, diagnostics: {} }, {
    ellipticSmoothing: { boundaryControl: 'wall-angle', maxSweeps: 0 },
    onMesh: (mesh, phase) => { phases.push(phase); if (phase === 'smoothing') throw new Error('mesh observer stopped'); },
  }), error => error.code === 'slor-observer-failed' && error.message === 'mesh observer stopped');
  assert.deepEqual(phases, ['initial', 'smoothing'], 'no further publication after observer cancellation');
  const grid = createEllipticStreamtubeGrid({ nodes: fixture(.8), massFlows, streamwiseCoordinates: xi, discretization: 'giles-1985' });
  for (const residual of [NaN, Infinity, -Infinity]) {
    const scalar = smoothEllipticStreamtubeGrid({ ...grid, residuals: () => ({ residual }) }, { maxSweeps: 0 });
    assert.equal(scalar.converged, false); assert.equal(scalar.termination.origin, 'exception');
    const paired = smoothPairedBoundaryGrid({ ...grid, coordinateEquations: { lineGrouping: 'boundary-pairs', lineSearch: 'armijo' },
      residuals: () => ({ residual, rows: [{ x: 0, y: 0 }] }) }, { maxSweeps: 0 });
    assert.equal(paired.converged, false); assert.equal(paired.termination.origin, 'exception');
  }
});

test('mesh preparation applies physical/chart gates to a partial initial guess and keeps the status incomplete', () => {
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 3, contourPanels: 40 });
  const system = createStreamtubeBodySystem(input), initial = system.initial, nodes = system.decode(initial).nodes;
  system.evaluate = () => { throw new Error('Partial mesh acceptance must not evaluate gas or flow.'); };
  const source = { input, system, initial, nodes, diagnostics: {}, guideField: {
    admissibleNode: () => true, diagnosticsForNodes: () => ({}),
  } };
  const result = prepareStreamtubeMesh(source, { ellipticSmoothing: { boundaryControl: 'wall-angle', maxSweeps: 0, tolerance: 1e-9 } });
  const report = result.diagnostics.gridSmoothing;
  assert.equal(report.initialGuessAccepted, true); assert.equal(report.converged, false);
  assert.ok(report.geometryAcceptance.valid && result.mesh.quality.valid);
  assert.ok(directStreamtubeVolumeGeometry(result.nodes).valid);
  assert.deepEqual(system.decode(result.initial).nodes, result.nodes);
  const deniedSystem = createStreamtubeBodySystem(input); let adopted = false;
  deniedSystem.adoptGeometry = () => { adopted = true; throw new Error('must not adopt a solid-crossing candidate'); };
  const denied = prepareStreamtubeMesh({ ...source, system: deniedSystem,
    guideField: { ...source.guideField, admissibleNode: () => false } }, {
    ellipticSmoothing: { boundaryControl: 'wall-angle', maxSweeps: 0, tolerance: 1e-9 },
  });
  assert.equal(denied.diagnostics.gridSmoothing.initialGuessAccepted, undefined);
  assert.equal(denied.diagnostics.gridSmoothing.retainedOriginal, true); assert.equal(adopted, false);
});

test('partial grid geometry checks reject overlaps and retain shared edges and positive staggered volumes', () => {
  const rectangle = offset => [0, 1, 2, 3].map(x => [0, 1, 2].map(y => ({ x, y: y + offset })));
  for (const nodes of [[rectangle(0)], [rectangle(0), rectangle(2)]]) {
    assert.ok(requireStreamtubeInitialGridDomain(nodes).valid);
    assert.ok(directStreamtubeVolumeGeometry(nodes).valid);
  }
  assert.throws(() => requireStreamtubeInitialGridDomain([rectangle(0), rectangle(0)]), /overlapping/);
  assert.throws(() => requireStreamtubeInitialGridDomain([rectangle(0), rectangle(.5)]), /crossed|overlapping/);
  const folded = rectangle(0); folded[1][1].y = 3;
  assert.throws(() => requireStreamtubeInitialGridDomain([folded]), /folded/);
});
