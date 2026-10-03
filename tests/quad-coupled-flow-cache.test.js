// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createQuadCoupledFlowCache } from '../src/ui/quad-coupled-flow-cache.js';

function result(mach = .2) {
  const families = { euler: 1e-12, boundaryLayer: 2e-12, edgeMatching: 3e-12 };
  const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
  return { model: 'research-streamtube-euler-bl', converged: true, mesh: { quality: { valid: true, minCornerSine: .2 }, vertices: nodes },
    sourceCase: { mach, alpha: 4, reynolds: 1e6, referenceChord: 1, flowModel: 'streamtube-grid', quadBoundaryLayers: true,
      gridIntervals: 16, gridTubes: 9, gridEllipticSmoothing: true, transitionMode: 'automatic', materialTrips: [[1, 1]],
      elements: [{ name: 'Main', points: [{ x: 0, y: 0 }, { x: 1, y: 0 }, { x: 0, y: .1 }] }] },
    checkpoint: { version: 1, families: { ...families }, continuation: { fractions: [[0, 1]], preferredOrdering: 'amd' },
      restart: { input: { mach, streamwiseMode: 'isentropic' }, options: { edgeMatching: 'section-velocity' },
        initialEuler: { x: new Float64Array([.1, .2, .3]), nodes, undisplacedNodes: structuredClone(nodes) },
        initialBL: new Float64Array([.1, .2, .3, .4, .5, .6, .7, .8]) } },
    families, referenceChord: 1, referenceReynolds: 1e6, solverLength: .98, kernelReynolds: 980000,
    solverSettings: { tolerance: 1e-10, maxIterations: 40, edgeMatching: 'section-velocity' },
    initialization: { euler: { gridSmoothing: { attempted: true, converged: true, regions: [{ history: Array(100).fill(1) }] } },
      attempts: [{ history: Array(100).fill(1) }], history: Array(1000).fill(1) },
    flow: { sections: Array(1000).fill(1) }, residual: Array(1000).fill(0), x: Array(1000).fill(1),
    boundaryLayer: { stations: Array(1000).fill(1) }, history: Array(1000).fill(1),
    refinement: { kind: 'paired-nested', level: 1, parent: { nx: 139, tubes: [12, 15, 12] },
      transfer: { history: Array(100).fill(1), nodes }, initialBL: Array(1000).fill(1) } };
}

test('same-Mach and Mach-only queries return one complete detached checkpoint with compact metadata', () => {
  const cache = createQuadCoupledFlowCache(), good = result(), before = structuredClone(good);
  assert.equal(cache.forCase(good.sourceCase), null); assert.equal(cache.remember(good), true);
  for (const mach of [.2, .3, .4]) {
    const parent = cache.forCase({ ...good.sourceCase, mach });
    assert.ok(parent); assert.equal(parent.sourceCase.mach, .2);
    assert.deepEqual(parent.checkpoint, good.checkpoint);
    assert.deepEqual(parent.mesh, { quality: good.mesh.quality });
    assert.deepEqual(parent.initialization.euler.gridSmoothing, { attempted: true, converged: true });
    assert.equal(parent.refinement.level, 1); assert.deepEqual(parent.refinement.parent.tubes, [12, 15, 12]);
    for (const key of ['flow', 'x', 'residual', 'boundaryLayer', 'history']) assert.equal(Object.hasOwn(parent, key), false);
    assert.equal(Object.hasOwn(parent.initialization, 'attempts'), false);
    assert.equal(Object.hasOwn(parent.refinement, 'initialBL'), false);
  }
  assert.deepEqual(good, before);
  assert.ok(JSON.stringify(cache.forCase(good.sourceCase)).length < JSON.stringify(good).length / 3);
});

test('cached private station policies preserve the exact checkpoint without adding a case control', () => {
  for (const policy of [undefined, 'station', 'station-auto']) {
    const cache = createQuadCoupledFlowCache(), good = result();
    if (policy !== undefined) good.checkpoint.continuation.linearOrdering = policy;
    if (policy === 'station-auto') good.checkpoint.continuation.stationFallback = true;
    const before = structuredClone(good);
    assert.equal(cache.remember(good), true);
    for (const mach of [.2, .3]) {
      const parent = cache.forCase({ ...good.sourceCase, mach });
      assert.deepEqual(parent.checkpoint, before.checkpoint);
      assert.equal(Object.hasOwn(parent.sourceCase, 'linearOrdering'), false);
      assert.equal(Object.hasOwn(parent.sourceCase, 'stationFallback'), false);
      parent.checkpoint.continuation.linearOrdering = 'changed';
      delete parent.checkpoint.continuation.stationFallback;
    }
    assert.deepEqual(cache.forCase(good.sourceCase).checkpoint, before.checkpoint);
    assert.deepEqual(good, before);
  }
});

test('remember/get detach nested controls, packed states, geometry and maintenance metadata', () => {
  const cache = createQuadCoupledFlowCache(), good = result(), original = structuredClone(good);
  assert.equal(cache.remember(good), true);
  good.sourceCase.elements[0].points[0].x = 9; good.checkpoint.restart.initialBL[0] = 9;
  good.checkpoint.continuation.fractions[0][0] = 9; good.solverSettings.tolerance = 1;
  let parent = cache.forCase(original.sourceCase);
  assert.deepEqual(parent.checkpoint, original.checkpoint); assert.equal(parent.solverSettings.tolerance, 1e-10);
  parent.sourceCase.alpha = 9; parent.checkpoint.restart.initialEuler.nodes[0][0][0].x = 8;
  parent.checkpoint.restart.initialEuler.x[0] = 8; parent.initialization.euler.gridSmoothing.converged = false;
  parent = cache.forCase(original.sourceCase);
  assert.deepEqual(parent.checkpoint, original.checkpoint); assert.equal(parent.sourceCase.alpha, 4);
  assert.equal(parent.initialization.euler.gridSmoothing.converged, true);
});

test('any other case edit is stale, restored controls regain eligibility, and property order does not matter', () => {
  const cache = createQuadCoupledFlowCache(), good = result(); cache.remember(good);
  const changes = [c => c.alpha++, c => c.reynolds++, c => c.gridIntervals++, c => c.gridTubes++,
    c => c.gridEllipticSmoothing = false, c => c.materialTrips[0][1] = .5,
    c => c.elements[0].points[1].x += .001, c => c.elements[0].name = 'Changed',
    c => c.quadBoundaryLayers = false, c => c.flowModel = 'subcritical', c => c.transitionMode = 'fixed-trip',
    c => c.addedControl = true, c => delete c.referenceChord];
  for (const change of changes) {
    const input = structuredClone(good.sourceCase); change(input); assert.equal(cache.forCase(input), null);
    assert.ok(cache.forCase(good.sourceCase));
  }
  const reordered = Object.fromEntries(Object.entries(good.sourceCase).reverse());
  assert.ok(cache.forCase(reordered));
  assert.equal(cache.forCase({ ...good.sourceCase, addedControl: undefined }), null);
  for (const mach of [0, -1, 1, NaN, Infinity, undefined]) assert.equal(cache.forCase({ ...good.sourceCase, mach }), null);
});

test('one good result replaces the prior seed; failed target or iterate cannot displace it', () => {
  const cache = createQuadCoupledFlowCache(), first = result(.2), second = result(.3);
  assert.equal(cache.remember(first), true); assert.equal(cache.remember(second), true);
  assert.equal(cache.forCase(first.sourceCase).sourceCase.mach, .3);
  const failed = result(.4); failed.converged = false; failed.stateConverged = true;
  failed.continuation = { reachedTarget: false, currentMach: .3, targetMach: .4 };
  assert.equal(cache.remember(failed), false);
  assert.equal(cache.forCase(first.sourceCase).sourceCase.mach, .3);
  cache.clear(); assert.equal(cache.forCase(second.sourceCase), null);
});

test('cheap completeness, gas-mode, residual, normalization and source-Mach guards preserve the accepted cache', () => {
  const cache = createQuadCoupledFlowCache(), good = result(); cache.remember(good);
  const mutations = [r => r.model = 'panel', r => r.mesh.quality.valid = false, r => r.sourceCase.mach = .3,
    r => r.sourceCase.quadBoundaryLayers = false, r => r.checkpoint.version = 2,
    r => delete r.checkpoint.continuation, r => r.checkpoint.restart.initialBL = [],
    r => r.checkpoint.restart.initialBL = [1, 2, 3], r => r.checkpoint.restart.initialEuler.x[0] = NaN,
    r => delete r.checkpoint.restart.initialEuler.undisplacedNodes,
    r => r.checkpoint.restart.initialEuler.nodes[0][0][0].x = Infinity,
    r => r.checkpoint.restart.input.streamwiseMode = 'momentum', r => r.checkpoint.restart.options.edgeMatching = 'other',
    r => r.checkpoint.restart.input.streamwiseMode = 'hybrid', r => r.families.euler = 2e-12,
    r => r.checkpoint.families.boundaryLayer = r.families.boundaryLayer = 1e-3,
    r => r.checkpoint.families.extra = 0, r => r.solverSettings.tolerance = 0,
    r => r.solverLength = 0, r => r.referenceReynolds = Infinity];
  for (const mutate of mutations) {
    const bad = result(); mutate(bad); assert.equal(cache.remember(bad), false);
    assert.deepEqual(cache.forCase(good.sourceCase).checkpoint, good.checkpoint);
  }
  const hybrid = result(.4); hybrid.checkpoint.restart.input.streamwiseMode = 'hybrid';
  hybrid.checkpoint.restart.options.blThermodynamics = 'historical-common-isentrope';
  assert.equal(cache.remember(hybrid), true); assert.equal(cache.forCase(good.sourceCase).sourceCase.mach, .4);
});

test('excluded flow payloads are never inspected or cloned and the cache imports no numerical module', () => {
  const cache = createQuadCoupledFlowCache(), good = result();
  for (const key of ['flow', 'x', 'residual', 'boundaryLayer', 'history']) Object.defineProperty(good, key, {
    get() { throw new Error('Large flow payload must not be read.'); }, enumerable: true });
  assert.equal(cache.remember(good), true); assert.ok(cache.forCase(good.sourceCase));
  const source = fs.readFileSync(new URL('../src/ui/quad-coupled-flow-cache.js', import.meta.url), 'utf8');
  assert.equal(/\bimport\b\s*(?:\(|[\w{*])/.test(source), false);
});

test('an intermediate Ncrit root cannot be mislabeled or cached as the requested condition', () => {
  const cache = createQuadCoupledFlowCache(), good = result();
  good.sourceCase.ncrit = 9; good.checkpoint.restart.options.ncrit = 9;
  assert.equal(cache.remember(good), true);
  const mislabeled = result(); mislabeled.sourceCase.ncrit = 9; mislabeled.checkpoint.restart.options.ncrit = 8;
  assert.equal(cache.remember(mislabeled), false);
  const intermediate = structuredClone(mislabeled); intermediate.sourceCase.ncrit = 8;
  intermediate.ncritContinuation = { actualNcrit: 8, targetNcrit: 9, reachedTarget: false };
  assert.equal(cache.remember(intermediate), false);
  assert.equal(cache.forCase(good.sourceCase).checkpoint.restart.options.ncrit, 9);
  delete intermediate.ncritContinuation;
  assert.equal(cache.remember(intermediate), true);
  assert.equal(cache.forCase(good.sourceCase), null);
  assert.equal(cache.forCase(intermediate.sourceCase).checkpoint.restart.options.ncrit, 8);
});

test('retained grid roots and mismatched source resolution cannot replace a valid cached parent', () => {
  const cache = createQuadCoupledFlowCache(), good = result();
  assert.equal(cache.remember(good), true);
  for (const field of ['gridSequence', 'automaticRefinement']) {
    const bad = result(.3);
    bad[field] = { kind: 'coarse-to-fine', reachedTarget: false, actualGridIntervals: 16, requestedGridIntervals: 128 };
    assert.equal(bad.converged, true, 'Exercise an incorrectly retained converged flag.');
    assert.equal(cache.remember(bad), false);
    bad[field].reachedTarget = true;
    assert.equal(cache.remember(bad), false, 'Resolution evidence still contradicts the requested grid.');
    bad[field].requestedGridIntervals = 16; bad.sourceCase.gridIntervals = 128;
    assert.equal(cache.remember(bad), false, 'Source-case count must describe the actual retained grid.');
    assert.deepEqual(cache.forCase(good.sourceCase).checkpoint, good.checkpoint);
  }
  const complete = result(.3);
  complete.gridSequence = { kind: 'coarse-to-fine', reachedTarget: true, actualGridIntervals: 16,
    requestedGridIntervals: 16, levels: [{ intervals: 8, converged: true }] };
  assert.equal(cache.remember(complete), true);
});
