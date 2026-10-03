// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { findUniqueGuideXBracket, sampleUniqueGuideX } from '../src/geometry/guide-x-bracket.js';
import { velocityAt } from '../src/inviscid/linear-vortex.js';
import { tracePotentialCurve } from '../src/inviscid/potential-curve.js';

const point = (x, potential) => ({ x, y: potential, potential });

test('request-local brackets preserve order and deduplicate ordinary shared vertices', () => {
  const path = [point(-4, 0), point(-2, 1), point(0, 2), point(2, 3), point(1, 4)];
  const before = structuredClone(path);
  assert.deepEqual(findUniqueGuideXBracket(path, -3), { vertex: null, lower: 0, upper: 1, orientation: 1 });
  assert.deepEqual(findUniqueGuideXBracket(path, -2), { vertex: 1, lower: 1, upper: 1, orientation: 0 });
  assert.throws(() => findUniqueGuideXBracket(path, 1.5), { code: 'GUIDE_X_AMBIGUOUS' });
  assert.throws(() => findUniqueGuideXBracket(path, 1), { code: 'GUIDE_X_AMBIGUOUS' });
  assert.throws(() => findUniqueGuideXBracket(path, 3), { code: 'GUIDE_X_UNRESOLVED' });
  assert.deepEqual(path, before);
});

test('vertical x-plane intervals and repeated crossings are rejected exactly', () => {
  const path = [point(-2, 0), point(0, 1), point(0, 2), point(2, 3), point(-2, 4)];
  assert.throws(() => findUniqueGuideXBracket(path, 0), /vertical guide interval/);
  assert.throws(() => findUniqueGuideXBracket(path, -1), /multiple resolved guide intersections/);
  assert.throws(() => findUniqueGuideXBracket([point(0, 1), point(1, 1)], .5), /strictly increasing potential/);
  assert.throws(() => findUniqueGuideXBracket([point(0, 1), point(NaN, 2)], .5), /finite points/);
});

test('nonlinear inverse uses either physical-x orientation and returns the actual curve', () => {
  // r(phi)=(phi^2, phi), with q=r'/|r'|^2, makes phi exactly the
  // velocity-potential parameter. Positive and negative phi branches have
  // opposite x orientations but the same closed-form inverse magnitude.
  for (const sign of [-1, 1]) {
    const values = sign < 0 ? [-3, -2, -1] : [1, 2, 3];
    const path = values.map(phi => ({ x: phi * phi, y: phi, potential: phi }));
    const sample = (_, phi) => ({ x: phi * phi, y: phi, potential: phi });
    const velocity = p => ({ u: 2 * p.y / (4 * p.y * p.y + 1), v: 1 / (4 * p.y * p.y + 1) });
    const result = sampleUniqueGuideX(path, 3, { sample, velocity, tolerance: 1e-13 });
    assert.ok(Math.abs(result.potential - sign * Math.sqrt(3)) < 2e-14);
    assert.ok(Math.abs(result.x - 3) < 1e-13);
    assert.equal(result.y, result.potential);
  }
});

test('a hidden continuous turn at a sampled extremum is not relabeled unique', () => {
  // x(phi)=1-(phi-.2)^2: x(0)=.96 is reached again at phi=.4,
  // although the sampled polyline [phi=-1,0,1] has one .96 vertex.
  const path = [-1, 0, 1].map(phi => ({ x: 1 - (phi - .2) ** 2, y: phi, potential: phi }));
  const sample = (_, phi) => ({ x: 1 - (phi - .2) ** 2, y: phi, potential: phi });
  const velocity = p => { const dx = -2 * (p.y - .2); return { u: dx / (dx * dx + 1), v: 1 / (dx * dx + 1) }; };
  assert.equal(findUniqueGuideXBracket(path, path[1].x).vertex, 1);
  assert.throws(() => sampleUniqueGuideX(path, path[1].x, { sample, velocity, tolerance: 1e-12 }),
    /unresolved x turning interval/);
});

test('safeguarded inverse enforces its residual and bounded iteration budget', () => {
  const path = [point(1, 1), point(4, 2)];
  const sample = (_, phi) => ({ x: phi * phi, y: phi, potential: phi });
  const velocity = p => ({ u: 2 * p.y / (4 * p.y * p.y + 1), v: 1 / (4 * p.y * p.y + 1) });
  assert.throws(() => sampleUniqueGuideX(path, 2, { sample, velocity, tolerance: 1e-15, maxIterations: 1 }),
    /did not converge/);
  const result = sampleUniqueGuideX(path, 2, { sample, velocity, tolerance: 1e-13 });
  assert.ok(Math.abs(result.x - 2) <= 1e-13);
});

test('frozen RAE 64/11 remote request is unique; actual near-LE ambiguities stay rejected', () => {
  const file = 'docs/rae2822/refined-guide-64x11.json', bytes = fs.readFileSync(file);
  const report = JSON.parse(bytes), c = report.failure.guideCapture;
  assert.equal(report.reproduced, true);
  assert.equal(createHash('sha256').update(bytes).digest('hex'), 'ba7dd1d8743291b8d9856f9e9566abd955e8c2c80f30f4068ec781a91d767fb5');
  const { path } = c;
  const bracket = findUniqueGuideXBracket(path, c.requestedX);
  assert.equal(bracket.orientation, 1);
  const velocity = p => velocityAt(p, c.field);
  const sample = (_, potential) => {
    const a = path[bracket.lower];
    const result = tracePotentialCurve({ seed: a, initialPotential: a.potential, endPotential: potential,
      velocity, tolerance: c.tolerance, maxSpatialStep: .1, maxStep: .05 });
    assert.equal(result.converged, true);
    return result.points.at(-1);
  };
  const result = sampleUniqueGuideX(path, c.requestedX, { sample, velocity, tolerance: c.tolerance });
  assert.ok(Math.abs(result.x - c.requestedX) <= c.tolerance);
  assert.throws(() => findUniqueGuideXBracket(path, .0034), { code: 'GUIDE_X_AMBIGUOUS' });
  assert.throws(() => sampleUniqueGuideX(path, path[116].x, { sample, velocity, tolerance: c.tolerance }),
    /unresolved x turning interval/);
  const initialLeading = c.suppliedInput.bodies[0].leadingIndex;
  // The matcher pins the inlet to the actual traced plane event; its
  // finite-tolerance endpoint is not resampled at the exact topology x.
  assert.ok(Math.abs(c.suppliedInput.cutPaths[0][0].x - path[0].x) <= c.tolerance);
  const inletXs = c.suppliedInput.cutPaths[0].slice(1, initialLeading).map(p => p.x);
  for (const x of inletXs) assert.doesNotThrow(() => findUniqueGuideXBracket(path, x));
  console.log(JSON.stringify({ frozenRaeSuppliedInletRequests: inletXs.length,
    actualRequestedXError: result.x - c.requestedX, permittedTolerance: c.tolerance }));
});
