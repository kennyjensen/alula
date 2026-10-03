// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { relaxStreamtubeInitialGrid } from '../src/euler/streamtube-elliptic-initializer.js';
import { reconstructOrthogonalBoundary } from '../src/geometry/orthogonal-boundary-control.js';
import { createEllipticStreamtubeGrid } from '../src/geometry/elliptic-streamtube-grid.js';
const xi = [0, .5, 1], eta = [0, .25, .65, 1], massFlows = [.25, .4, .35];
const fixture = curvature => [-.1, 0, .1].map(x => eta.map(y => ({ x, y: .5 * curvature * x * x + y })));
const prepared = nodes => ({ nodes: [nodes], initial: [], system: { layout: { primaryBody: 0, bodies: [] },
  decode: () => ({ allocation: { groups: [massFlows.map(massFlow => ({ massFlow }))] } }) } });
const controls = { seed: 'supplied', boundaryControl: 'wall-angle', streamwiseCoordinates: xi, maxSweeps: 200, tolerance: 1e-9 };

test('an infeasible additional angle target gets a certified harmonic SLOR solve with unchanged fixed boundaries', () => {
  const nodes = fixture(10), before = structuredClone(nodes), events = [];
  const result = relaxStreamtubeInitialGrid(prepared(nodes), { ...controls, onSweep: (h, n) => events.push({ h, n }) });
  assert.ok(result.converged);
  const report = result.regions[0], original = report.harmonicFallback.originalAttempt;
  assert.equal(original.code, 'orthogonal-normal-branch'); assert.equal(original.stage, 'boundary-control-initialization');
  assert.ok(original.diagnostics.discriminantRatio < 0); assert.equal(original.diagnostics.station, 1);
  assert.ok(events.length > 1 && events.every(e => e.h.harmonicFallback && e.h.invalidCells === 0));
  assert.equal(report.spacingFallback, undefined); assert.equal(report.harmonicFallback.retainedAdditionalPoissonSpacingControl, false);
  assert.equal(report.coordinateEquations.xi, 'Laplace'); assert.equal(report.coordinateEquations.eta, 'Laplace');
  const plain = createEllipticStreamtubeGrid({ nodes, massFlows, streamwiseCoordinates: xi, discretization: 'giles-1985' });
  assert.ok(plain.residuals(result.nodes[0]).residual <= controls.tolerance);
  for (let i = 0; i < nodes.length; i++) for (let j = 0; j < nodes[i].length; j++) if (!i || i === 2 || !j || j === 3)
    assert.deepEqual(result.nodes[0][i][j], before[i][j]);
  assert.deepEqual(nodes, before);
});

test('the fallback preserves fold/nonfinite rejection and is absent when the original angle solve works', () => {
  const flat = fixture(0), ordinary = relaxStreamtubeInitialGrid(prepared(flat), controls);
  assert.ok(ordinary.converged); assert.equal(ordinary.regions[0].harmonicFallback, undefined);
  assert.deepEqual(ordinary.nodes[0], flat);
  const nodes = fixture(10); nodes[1][1].x = -.5;
  let publications = 0;
  const rejected = relaxStreamtubeInitialGrid(prepared(nodes), { ...controls, onSweep: () => publications++ });
  assert.equal(rejected.converged, false); assert.equal(rejected.regions[0].quality.valid, false); assert.equal(publications, 0);
  assert.throws(() => reconstructOrthogonalBoundary({ point: { x: 0, y: 0 }, tangent: { x: 1, y: 0 },
    secondTangent: { x: 0, y: 1e308 }, firstInterior: { x: 0, y: 1e308 }, secondInterior: { x: 0, y: 1e308 },
    firstDistance: .1, secondDistance: .2 }), e => e.code === undefined);
});
