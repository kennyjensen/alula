// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { pathToFileURL } from 'node:url';
import { resolve } from 'node:path';
import { findUniqueGuideXBracket, sampleUniqueGuideX } from '../src/geometry/guide-x-bracket.js';
import { pointInside, segmentsTouch } from '../src/geometry/airfoil.js';
import { velocityAt } from '../src/inviscid/linear-vortex.js';
import { tracePotentialCurve } from '../src/inviscid/potential-curve.js';
import { potentialDifference } from '../src/inviscid/streamfunction.js';
import { createStagnationGuideConnector } from '../src/inviscid/stagnation-guide-connector.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/rae-mit-max-resolution-guide.json', import.meta.url)));
const close = (a, b, tolerance) => assert.ok(Math.abs(a - b) <= tolerance,
  `Absolute error ${Math.abs(a - b)} exceeds ${tolerance}`);
const fluid = p => !fixture.bodies.some(body => pointInside(p, body));
const segment = (a, b) => !fixture.field.panels.some(p => segmentsTouch(a, b, p.a, p.b, 1e-12));
const anchorDirection = { x: -fixture.anchorDerivative.y, y: fixture.anchorDerivative.x };

test('frozen MIT maximum-resolution request has two intersections and a unique forward source-panel join', () => {
  for (const binding of [fixture.provenance.capture, fixture.provenance.publicInput, fixture.provenance.oldInitializer])
    assert.equal(createHash('sha256').update(fs.readFileSync(binding.path)).digest('hex'), binding.sha256);
  assert.deepEqual(JSON.parse(fs.readFileSync(fixture.provenance.publicInput.path)), fixture.publicInput);
  assert.equal(fixture.publicInput.gridIntervals, 128);
  assert.equal(fixture.publicInput.gridTubes, 11);
  assert.equal(fixture.field.panels.length, 128);
  assert.equal(fixture.bodies[0].length, 129);
  assert.throws(() => findUniqueGuideXBracket(fixture.path, fixture.requestedX), error => {
    assert.equal(error.code, 'GUIDE_X_AMBIGUOUS');
    assert.deepEqual(error.diagnostics.intervals, [[114, 115], [121, 122]]);
    assert.deepEqual(error.diagnostics.vertices, []);
    return true;
  });
  assert.deepEqual(findUniqueGuideXBracket(fixture.path, fixture.joinX),
    { vertex: null, lower: 108, upper: 109, orientation: 1 });
  const maximumX = Math.max(...fixture.path.map(p => p.x));
  assert.ok(fixture.path.at(-1).x < fixture.requestedX && fixture.requestedX < maximumX);
  assert.ok(maximumX < fixture.anchor.x - fixture.tolerance);
  const panel = fixture.field.panels[fixture.sheetStagnation.segment];
  assert.equal(fixture.joinX, Math.min(panel.a.x, panel.b.x));
  assert.ok(fixture.joinX + fixture.tolerance < fixture.requestedX);
});

function checkConnector(c) {
  assert.ok(c.diagnostics.certificateLeaves > 0);
  assert.ok(c.diagnostics.minimumPotentialDerivative > 0);
  assert.ok(Math.abs(c.diagnostics.joinPotentialError) <= 8 * fixture.tolerance);
  assert.match(c.diagnostics.interpretation, /not an exact constant-streamfunction/);
  assert.ok(c.diagnostics.maximumStreamfunctionDefect > 1e-6,
    'Keep the declared geometric approximation visible.');
  const controls = c.diagnostics.controls, q0 = velocityAt(c.join, fixture.field);
  close((controls[1].y - controls[0].y) / (controls[1].x - controls[0].x), q0.v / q0.u, 2e-14);
  close((controls[3].y - controls[2].y) / (controls[3].x - controls[2].x),
    anchorDirection.y / anchorDirection.x, 2e-14);
  assert.deepEqual(c.anchor, { ...fixture.anchor, potential: fixture.stagnationPotential });
  assert.equal(fluid(c.join), true);
  let previous = c.join;
  for (const x of [...fixture.nearLeadingEdgeRequests, fixture.anchor.x]) {
    const p = c.atX(x);
    close(p.x, x, 2e-18);
    assert.equal(fluid(p), true);
    assert.equal(segment(previous, p), true);
    assert.ok(p.potential > previous.potential && p.derivative > 0);
    // These chords are checked before using the panel potential primitive.
    assert.equal(segment(p, c.anchor), true);
    close(c.anchor.potential - p.potential, potentialDifference(p, c.anchor, fixture.field), 3e-14);
    previous = p;
  }
  for (const x of [fixture.nearLeadingEdgeRequests[0], fixture.requestedX,
    fixture.nearLeadingEdgeRequests.at(-1), fixture.anchor.x]) {
    const p = c.atX(x);
    close(c.atPotential(p.potential).x, p.x, fixture.tolerance);
  }
}

test('MIT local reintegration and unchanged connector fill every refined near-LE request with exterior ordered potential', t => {
  const before = structuredClone(fixture), started = performance.now();
  const velocity = p => velocityAt(p, fixture.field);
  const bracket = findUniqueGuideXBracket(fixture.path, fixture.joinX);
  const sample = (_, potential) => {
    const a = fixture.path[bracket.lower];
    const result = tracePotentialCurve({ seed: a, initialPotential: a.potential, endPotential: potential,
      velocity, tolerance: fixture.tolerance, maxSpatialStep: .1, maxStep: .05, minStep: 1e-12,
      admissible: fluid, admissibleSegment: segment });
    assert.equal(result.converged, true);
    return result.points.at(-1);
  };
  const join = sampleUniqueGuideX(fixture.path, fixture.joinX, { sample, velocity, tolerance: fixture.tolerance });
  close(join.x, fixture.joinX, fixture.tolerance);
  checkConnector(createStagnationGuideConnector({ join, anchor: fixture.anchor, anchorDirection,
    stagnationPotential: fixture.stagnationPotential, field: fixture.field, tolerance: fixture.tolerance,
    streamfunctionLevel: fixture.streamfunctionLevel, admissible: fluid, admissibleSegment: segment }));
  assert.deepEqual(fixture, before);
  t.diagnostic(JSON.stringify({ localKernelSeconds: (performance.now() - started) / 1000,
    refinedRequests: fixture.nearLeadingEdgeRequests.length, fullMeshBuilds: 0, flowSolves: 0 }));
});

function replaceOnce(source, needle, replacement) {
  assert.equal(source.split(needle).length, 2, 'The validation-only insertion point must be unique.');
  return source.replace(needle, replacement);
}

function moduleURL(source, original, replacements = {}) {
  const location = pathToFileURL(resolve(original));
  return 'data:text/javascript;base64,' + Buffer.from(source.replace(/\bfrom\s+(['"])(\.[^'"]+)\1/g,
    (_, quote, name) => `from ${quote}${replacements[name] ?? new URL(name, location).href}${quote}`)).toString('base64');
}

async function stoppedPublicPrefix(initializerSource) {
  const original = 'src/euler/streamtube-body-initializer.js';
  // Only record an existing exception, and stop immediately after the
  // unmodified production connector is certified and installed. No inputs,
  // numerical branches, field construction, or matcher arithmetic change.
  let source = replaceOnce(initializerSource,
    '      const profile = profiles[body], anchor = profile?.curve.evaluate(profile.stag).point;',
    `      const profile = profiles[body], anchor = profile?.curve.evaluate(profile.stag).point;
      error.testGuideRequest = { x, path: structuredClone(path) };`);
  source = replaceOnce(source, '          ...connector.diagnostics });',
    `          ...connector.diagnostics });
        throw Object.assign(new Error('Validation stop after connector publication.'), {
          code: 'TEST_CONNECTOR_PUBLISHED', connector, originalCause: error.cause,
          requestedX: error.requestedX, originalPath: path, connected: profile.upstream,
          published: diagnostics.geometricStagnationConnectors, panelBoundaryCondition });`);
  source = replaceOnce(source, '    ({ guides, outer } = matched);',
    `    throw new Error('Validation stop: matcher completed without publishing the expected connector.');
    ({ guides, outer } = matched);`);
  const initializer = moduleURL(source, original), adapter = 'src/euler/streamtube-result.js';
  const { prepareStreamtubeAssembly } = await import(moduleURL(fs.readFileSync(adapter, 'utf8'), adapter,
    { './streamtube-body-initializer.js': initializer }));
  const input = structuredClone(fixture.publicInput), before = structuredClone(input);
  let failure, meshPublications = 0;
  try {
    prepareStreamtubeAssembly(input, { meshOnly: true, onMesh: () => {
      meshPublications++;
      throw new Error('Unexpected mesh publication in bounded prefix test.');
    } });
  } catch (error) { failure = error; }
  assert.equal(meshPublications, 0);
  assert.deepEqual(input, before);
  assert.ok(failure);
  return failure;
}

test('actual public maximum-resolution routing fails before recovery in archived source and publishes the certified connector now', async t => {
  const started = performance.now();
  const old = await stoppedPublicPrefix(fs.readFileSync(fixture.provenance.oldInitializer.path, 'utf8'));
  assert.equal(old.code, 'streamtube-guide-trace');
  assert.match(old.message, /left the admissible fluid region/);
  assert.equal(old.testGuideRequest.x, fixture.requestedX);
  assert.deepEqual(old.testGuideRequest.path, fixture.path);

  const current = await stoppedPublicPrefix(fs.readFileSync('src/euler/streamtube-body-initializer.js', 'utf8'));
  assert.equal(current.code, 'TEST_CONNECTOR_PUBLISHED');
  assert.equal(current.panelBoundaryCondition, 'streamfunction');
  assert.equal(current.originalCause.code, 'GUIDE_X_AMBIGUOUS');
  assert.equal(current.requestedX, fixture.requestedX);
  assert.deepEqual(current.originalPath, fixture.path);
  assert.equal(current.published.length, 1);
  assert.equal(current.published[0].sourcePanel, fixture.sheetStagnation.segment);
  checkConnector(current.connector);
  const prefix = fixture.path.filter(p => p.potential < current.connector.join.potential);
  assert.deepEqual(current.connected.slice(0, prefix.length), prefix);
  assert.deepEqual(current.connected.slice(prefix.length), current.connector.points);
  assert.ok(current.connected.every((p, k, points) => !k
    || p.x > points[k - 1].x && p.potential > points[k - 1].potential));
  t.diagnostic(JSON.stringify({ publicPrefixSeconds: (performance.now() - started) / 1000,
    publicPrefixRequests: 2, meshPublications: 0, smoothingCalls: 0, flowSolves: 0 }));
});
