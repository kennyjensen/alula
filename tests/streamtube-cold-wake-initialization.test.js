import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { createHash } from 'node:crypto';
import { coldStreamtubeWakeCorrespondenceEligible, initializeCoupledStreamtubeBody, panelCoupledEdgeGuess, scaleStreamtubeBLThicknesses } from '../src/euler/streamtube-coupled-initializer.js';
import { naca4Standard } from '../src/geometry/airfoil.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('panel BL edge guesses preserve surface and wake speeds when physical units change', () => {
  const points = naca4Standard('0012', 160), bodies = [{ element: 0 }];
  points.push({ ...points[0] });
  const bl = { surfaces: [{ body: 0, side: 'upper', ids: [0, 1] }, { body: 0, side: 'lower', ids: [2, 3] }],
    wakes: [{ ids: [4, 5, 6] }] };
  const geometry = { coordinates: [.2, 1.01, .3, 1.02, 1.02, 1.12, 1.62].map(s => ({ s })) };
  const stations = geometry.coordinates.map((_, id) => ({ id, body: 0,
    kind: id < 4 ? 'surface' : 'wake', side: id < 2 ? 'upper' : 'lower' }));
  const speeds = [1, 3].map(scale => {
    const input = { elements: [{ points: points.map(p => ({ x: p.x * scale, y: p.y * scale })),
      trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 160 } }],
      referenceChord: scale, alpha: 4, mach: .01, reynolds: 1e6, ncrit: 9, materialTrips: [[1, 1]] };
    const before = structuredClone(input), guess = panelCoupledEdgeGuess(input, bodies, .9958 * scale);
    assert.equal(guess.accepted, true); assert(guess.residual < 1e-8);
    assert.deepEqual(input, before);
    return stations.map(station => guess.initialEdgeVelocity(station, bl, geometry));
  });
  speeds[0].forEach((speed, i) => {
    assert(speed > 0 && Number.isFinite(speed));
    assert(Math.abs(speed - speeds[1][i]) < 1e-9, `station ${i} changes with physical units`);
  });
});

test('finite-base NACA 64x24 at Mach .01 initializes its full cold BL without thinning or altering the precursor', () => {
  // Exact app geometry/settings from the reported failure, after its Euler
  // precursor has converged. No Newton solve or generated substitute mesh.
  const { input, options } = JSON.parse(gunzipSync(fs.readFileSync(
    new URL('./fixtures/naca64x24-low-mach-euler.json.gz', import.meta.url))));
  const before = structuredClone({ input, options });
  const seed = createCoupledStreamtubeBody(input, { ...options, blInitialization: 'mrchue' });
  const prepared = initializeCoupledStreamtubeBody(input, options, { maximumBacktracks: 0 });
  const { system, initialization, mesh } = prepared;
  assert.equal(input.mach, .01);
  assert.equal(input.alpha, 4);
  assert.equal(system.bl.hasFiniteBase, true);
  assert.equal(initialization.thicknessFactor, 1);
  assert.equal(initialization.history.length, 1);
  const adjustment = initialization.history[0].wakeCorrespondence;
  assert.equal(adjustment.beforeQuality.valid, false);
  assert.equal(adjustment.accepted, true);
  assert.equal(adjustment.geometry.equationsChanged, false);
  assert(adjustment.geometry.maximumGapChange < 1e-14);
  assert(adjustment.geometry.maximumCenterChange < 1e-14);
  assert.deepEqual(system.initial.slice(system.ne), seed.initial.slice(seed.ne));
  assert.deepEqual(system.bl.snapshotActive(), seed.bl.snapshotActive());
  assert.equal(mesh.quality.valid, true);
  assert(mesh.quality.minCornerSine > .2);
  assert.equal(initialization.flowSolved, false);
  assert(system.evaluate(system.initial).families.boundaryLayer < .04);
  assert.deepEqual({ input, options }, before);
});

test('panel BL thickness transfer removes the lower aft-body kink before the coupled solve', () => {
  const { input, options } = JSON.parse(gunzipSync(fs.readFileSync(
    new URL('./fixtures/naca64x24-low-mach-euler.json.gz', import.meta.url))));
  const points = naca4Standard('0012', 160); points.push({ ...points[0] });
  const physical = { elements: [{ points, trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 160 } }],
    referenceChord: 1, alpha: 4, mach: .01, reynolds: 1e6, ncrit: 9, materialTrips: [[1, 1]] };
  const lengthScale = createCoupledStreamtubeBody(input, options).euler.conditions.lengthScale;
  const guess = panelCoupledEdgeGuess(physical, input.bodies, lengthScale);
  assert.equal(guess.accepted, true);
  const settings = { ...options, edgeMatching: 'pressure', initialEdgeVelocity: guess.initialEdgeVelocity };
  delete settings.blThermodynamics;
  const baseline = initializeCoupledStreamtubeBody(input, settings);
  const before = baseline.system.initial.slice(), phases = baseline.system.bl.snapshotActive();
  const mapped = guess.initialBoundaryLayer(baseline.system);
  // A geometric startup backtrack must not leave a thin stagnation prefix
  // attached to the full panel profile when that full profile is retried.
  const thinState = baseline.system.initial.slice();
  for (let i = baseline.system.ne; i < thinState.length; i += 4) {
    thinState[i + 1] *= .125; thinState[i + 2] *= .125;
  }
  const restored = guess.initialBoundaryLayer({ ...baseline.system, initial: thinState }, { thicknessFactor: .125 });
  mapped.initialBL.forEach((v, i) => assert(Math.abs(v - restored.initialBL[i]) < 1e-10 * Math.max(1, Math.abs(v)),
    `Mapped BL slot ${i} depends on the discarded startup thickness factor`));
  assert.deepEqual(baseline.system.initial, before);
  assert.deepEqual(baseline.system.bl.snapshotActive(), phases);
  const scaled = { ...physical, referenceChord: 3,
    elements: physical.elements.map(e => ({ ...e, points: e.points.map(p => ({ x: 3 * p.x, y: 3 * p.y })) })) };
  const mappedScaled = panelCoupledEdgeGuess(scaled, input.bodies, 3 * lengthScale).initialBoundaryLayer(baseline.system);
  assert.deepEqual(mappedScaled.transitionState, mapped.transitionState);
  mapped.initialBL.forEach((v, i) => assert(Math.abs(v - mappedScaled.initialBL[i]) < 1e-7 * Math.max(1, Math.abs(v)),
    `Mapped BL slot ${i} changes with physical units`));
  const prepared = initializeCoupledStreamtubeBody(input, { ...settings, ...mapped }, { maximumBacktracks: 0 });
  assert.equal(prepared.initialization.thicknessFactor, 1);
  assert.equal(prepared.mesh.quality.valid, true);
  assert(prepared.mesh.quality.minCornerSine > .2);
  const turn = system => {
    const rows = system.evaluate(system.initial).outer.nodes[0], j = rows[0].length - 2;
    let maximum = 0;
    for (let i = 1; i < rows.length - 1; i++) {
      const a = rows[i - 1][j], b = rows[i][j], c = rows[i + 1][j];
      if (b.x < .95 || b.x > .99) continue;
      const ux = b.x - a.x, uy = b.y - a.y, vx = c.x - b.x, vy = c.y - b.y;
      maximum = Math.max(maximum, Math.abs(Math.atan2(ux * vy - uy * vx, ux * vx + uy * vy)) * 180 / Math.PI);
    }
    return maximum;
  };
  assert(turn(baseline.system) > 20, 'reproduce the inverse/direct thickness discontinuity');
  assert(turn(prepared.system) < 1, 'mapped panel thickness must remove the initial near-wall kink');
  assert(prepared.system.evaluate(prepared.system.initial).residual.every(Number.isFinite));
});

test('saved fine-grid cell locations distinguish immovable surface defects from first and second wake intervals', () => {
  const layout = { nx: 193, elements: 1, tubes: [14, 14], bodies: [{ trailingIndex: 129 }], independentWakeBanks: true };
  const eligible = invalidCells => coldStreamtubeWakeCorrespondenceEligible({ layout, quality: { valid: false, invalidCells } });
  // Frozen 32/11 source: i126 upper surface is upstream of the TE at129.
  assert.equal(eligible([4466, 4467, 4468, 4469, 4470, 4471, 4508]), false);
  // Half-thickness defects span several tubes and TWO wake intervals.
  assert.equal(eligible([1833, 4508, 4509, 4510, 4522, 4526]), true);
  assert.equal(eligible([0, 4508]), false);
  for (const ids of [[], [-1], [5404], [4508.5], [NaN]]) assert.equal(eligible(ids), false);
  assert.equal(coldStreamtubeWakeCorrespondenceEligible({ layout, quality: { valid: true, invalidCells: [4508] } }), false);
  assert.equal(coldStreamtubeWakeCorrespondenceEligible({ layout: { ...layout, independentWakeBanks: false }, quality: { valid: false, invalidCells: [4508] } }), false);
});

test('multielement eligibility uses each passage and its adjoining wake start', () => {
  const layout = { nx: 12, elements: 2, tubes: [2, 3, 4], bodies: [{ trailingIndex: 5 }, { trailingIndex: 9 }], independentWakeBanks: true };
  const eligible = ids => coldStreamtubeWakeCorrespondenceEligible({ layout, quality: { valid: false, invalidCells: ids } });
  assert.equal(eligible([10, 24 + 5 * 3, 60 + 9 * 4]), true);
  // The upper passage cannot be repaired using the other body's earlier TE.
  assert.equal(eligible([60 + 8 * 4]), false);
  assert.equal(eligible([24 + 4 * 3]), false);
});

test('an already valid real coupled seed is identical to the archived initializer, including all residual rows', async () => {
  const moduleURL = new URL('../src/euler/streamtube-coupled-initializer.js', import.meta.url);
  const archived = fs.readFileSync(new URL('../docs/rae2822/cold-wake-initialization/streamtube-coupled-initializer.before.js.txt', import.meta.url), 'utf8')
    .replace(/from '(\.\/[^']+)'/g, (_, path) => `from '${new URL(path.replace('./streamtube-displacement-initializer.js', './streamtube-displacement.js'), moduleURL).href}'`);
  const old = await import(`data:text/javascript;base64,${Buffer.from(archived).toString('base64')}`);
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), wakeGeometry: 'independent-banks' };
  const base = createCoupledStreamtubeBody(input), x = base.initial.slice(0, base.ne);
  const options = { reynolds: 1e6, ncrit: 9, initialBL: base.initial.slice(base.ne), initialEuler: { x, ...base.euler.decode(x) } };
  const before = structuredClone({ input, options });
  const expected = old.initializeCoupledStreamtubeBody(input, options, { maximumBacktracks: 0 });
  const actual = initializeCoupledStreamtubeBody(input, options, { maximumBacktracks: 0 });
  assert.equal(actual.initialization.history.length, 1);
  assert.equal(actual.initialization.history[0].wakeCorrespondence, undefined);
  assert.deepEqual(actual.initialization, expected.initialization);
  assert.deepEqual(actual.mesh, expected.mesh);
  assert.deepEqual(actual.system.initial, expected.system.initial);
  assert.deepEqual(actual.system.evaluate(actual.system.initial).residual, expected.system.evaluate(expected.system.initial).residual);
  assert.deepEqual({ input, options }, before);
});

const serial = value => JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v);
const digest = value => createHash('sha256').update(serial(value)).digest('hex');
const frozen = JSON.parse(fs.readFileSync(new URL('./fixtures/cold-wake-selection.json', import.meta.url)));

// This older receipt predates native transition/wake shear reseeding. Preserve
// its physical theta/delta/Ue/grid identity; retain its auxiliary/residual hashes
// as historical evidence and compare ALL current fields against the pre-swap
// initializer with precisely the same numerical dependencies.
const historicalPath = new URL(`../${frozen.provenance.sources.coldExpected.path}`, import.meta.url);
const historicalText = fs.readFileSync(historicalPath, 'utf8');
assert.equal(createHash('sha256').update(historicalText).digest('hex'), frozen.provenance.sources.coldExpected.sha256);
const historicalCold = JSON.parse(historicalText);
assert.equal(digest(historicalCold.restart.initialBL), frozen.cold.expected.packedBL);
assert.equal(digest(historicalCold.residual), frozen.cold.expected.residual);
assert.equal(digest(historicalCold.initialization.history), frozen.cold.expected.history);
assert.equal(digest(historicalCold.mesh), frozen.cold.expected.mesh);

let beforeWakeFirst;
async function initializerBeforeWakeFirst() {
  if (!beforeWakeFirst) {
    const moduleURL = new URL('../src/euler/streamtube-coupled-initializer.js', import.meta.url);
    const original = fs.readFileSync(new URL('../docs/solver-reliability/te-wake-seed-improvement/initializer.before.js.txt', import.meta.url), 'utf8');
    assert.equal(createHash('sha256').update(original).digest('hex'),
      '8e15c8e8fc1d605bd7bba7c0f5f1992423d4afb9dd7d28baf33091a1a83c2c47');
    const bound = original.replace(/from '(\.\/[^']+)'/g, (_, path) => `from '${new URL(path.replace('./streamtube-displacement-initializer.js', './streamtube-displacement.js'), moduleURL).href}'`);
    beforeWakeFirst = import(`data:text/javascript;base64,${Buffer.from(bound).toString('base64')}`);
  }
  return beforeWakeFirst;
}

function verifyFrozenSeed(prepared, expected, baseline) {
  const { system, initialization, mesh } = prepared, value = system.evaluate(system.initial);
  assert.equal(initialization.thicknessFactor, expected.thicknessFactor);
  assert.deepEqual(system.bl.snapshotActive(), expected.phase);
  assert.equal(mesh.quality.valid, true);
  assert.equal(initialization.flowSolved, false);
  assert.equal(initialization.equationsChanged, false);
  const frozenFields = { eulerState: system.initial.slice(0, system.ne),
    physicalNodes: value.outer.nodes, undisplacedNodes: value.outer.undisplacedNodes };
  for (const [name, data] of Object.entries(frozenFields))
    assert.equal(digest(data), expected[name], `Frozen ${name} must remain exactly unchanged.`);
  const packed = system.initial.slice(system.ne);
  assert.equal(packed.length, historicalCold.restart.initialBL.length);
  for (let k = 0; k < packed.length; k++) if (k % 4 !== 0)
    assert.equal(packed[k], historicalCold.restart.initialBL[k], `Historical theta/delta/Ue slot ${k} changed.`);
  // No tolerance: the ordering change must preserve every current auxiliary,
  // density, coordinate, phase and residual row in the accepted fallback.
  const baselineValue = baseline.system.evaluate(baseline.system.initial);
  assert.deepEqual(system.initial, baseline.system.initial);
  assert.deepEqual(system.bl.snapshotActive(), baseline.system.bl.snapshotActive());
  assert.deepEqual(mesh, baseline.mesh);
  assert.deepEqual(value.outer.nodes, baselineValue.outer.nodes);
  assert.deepEqual(value.outer.undisplacedNodes, baselineValue.outer.undisplacedNodes);
  assert.deepEqual(value.residual, baselineValue.residual);
  assert.deepEqual(value.families, baselineValue.families);
  assert.deepEqual(value.families, expected.families);
  const history = initialization.history;
  assert.equal(history.length, 3);
  assert.equal(initialization.originalBLInitialization.boundaryLayer.some(row =>
    row.method === 'mrchue' && row.localConvergenceWarnings?.length > 0), false);
  assert.deepEqual(history, baseline.initialization.history);
  return value;
}

for (const missing of ['omitted', 'null']) test(`unsupplied ${missing} cold BL with a completed surface march preserves historical physical fields and the current original fallback sequence`, async () => {
  const { input, options, expected } = structuredClone(frozen.cold);
  assert.equal(options.transitionMode, 'automatic');
  assert.equal(Object.hasOwn(options, 'initialBL'), false);
  if (missing === 'null') options.initialBL = null;
  const before = digest({ input, options });
  // Uses only the frozen converged Euler state. These are initialization
  // calls, not Euler or coupled Newton solves. The old full-thickness seed
  // became convex after alignment but had local M=.998 and a poor BL basin.
  const old = await initializerBeforeWakeFirst();
  const baseline = old.initializeCoupledStreamtubeBody(structuredClone(input), structuredClone(options));
  const prepared = initializeCoupledStreamtubeBody(input, options);
  const value = verifyFrozenSeed(prepared, expected, baseline);
  assert.deepEqual(prepared.initialization.history.map(h => [h.thicknessFactor, h.wakeInitialization ?? null, h.accepted]),
    [[1, null, false], [1, 'iset-linear-shape', false], [.5, null, true]]);
  assert.ok(prepared.initialization.history.every(h => h.wakeCorrespondence === undefined));
  assert.ok(value.outer.diagnostics.maxMach < .32);
  assert.equal(digest({ input, options }), before);
});

test('supplied mapped fine-grid BL retains physical fields while pairing both wake banks forward at full thickness', () => {
  const { input, options, expected } = structuredClone(frozen.fine), before = digest({ input, options });
  assert.equal(options.transitionMode, 'automatic');
  assert.ok(options.initialBL.length > 0);
  const prepared = initializeCoupledStreamtubeBody(input, options, { maximumBacktracks: 0 });
  const { system, initialization, mesh } = prepared, value = system.evaluate(system.initial);
  assert.equal(initialization.thicknessFactor, expected.thicknessFactor);
  assert.deepEqual(system.bl.snapshotActive(), expected.phase);
  assert.equal(mesh.quality.valid, true);
  assert.equal(initialization.flowSolved, false);
  assert.equal(initialization.equationsChanged, false);
  // The previous fixture's exact physical-node/residual/history hashes bind
  // the archived constant-local-tau rule. The frame-consistent rule changes
  // those coordinates, while retaining the density/globals and full BL state.
  // Omitted/null cold guesses above preserve historical physical fields and
  // all current pre-swap state/residual rows exactly, including auxiliary shear.
  assert.equal(digest(system.initial.slice(0, system.ne)), expected.eulerState);
  assert.equal(digest(system.initial.slice(system.ne)), expected.packedBL);
  // The saved checkpoint encodes derivative Maps as entry arrays.
  const allocation = JSON.parse(JSON.stringify(value.outer.allocation, (_, v) => v instanceof Map ? [...v] : v));
  assert.deepEqual(allocation, options.initialEuler.allocation);
  assert.ok(Array.from(value.residual).every(Number.isFinite));
  assert.equal(prepared.initialization.history.length, 1);
  const adjustment = prepared.initialization.history[0].wakeCorrespondence;
  assert.equal(adjustment.beforeQuality.valid, false);
  assert.equal(adjustment.afterQuality.valid, true);
  assert.equal(adjustment.accepted, true);
  assert.equal(adjustment.thicknessChanged, false);
  assert.equal(adjustment.geometry.coordinateRule, 'equal-center-edge-projected-bank-advance');
  const { nodes } = value.outer, { layout } = system.euler;
  const center = (b, i) => ({ x: .5 * (nodes[b][i].at(-1).x + nodes[b + 1][i][0].x),
    y: .5 * (nodes[b][i].at(-1).y + nodes[b + 1][i][0].y) });
  const thicknesses = system.bl.thicknesses(system.initial.slice(system.ne));
  for (let b = 0; b < layout.elements; b++) for (let i = layout.bodies[b].trailingIndex + 1; i <= layout.nx; i++) {
    const left = center(b, i - 1), c = center(b, i), right = center(b, Math.min(i + 1, layout.nx));
    const dx = c.x - left.x, dy = c.y - left.y, advance = Math.hypot(dx, dy);
    assert.ok(advance > 0);
    for (const [g, j] of [[b, layout.tubes[b]], [b + 1, 0]]) {
      const p = nodes[g][i][j], q = nodes[g][i - 1][j];
      const forward = ((p.x - q.x) * dx + (p.y - q.y) * dy) / advance;
      assert.ok(forward > 0);
      assert.ok(Math.abs(forward - advance) < 1e-13);
    }
    // Independent normal-gap check uses the actual centered physical chord.
    const tx = right.x - left.x, ty = right.y - left.y, length = Math.hypot(tx, ty);
    const lower = nodes[b][i].at(-1), upper = nodes[b + 1][i][0];
    const gap = ((upper.y - lower.y) * tx - (upper.x - lower.x) * ty) / length;
    assert.ok(Math.abs(gap - thicknesses.wakes[b][i - layout.bodies[b].trailingIndex - 1]) < 1e-13);
  }
  // Every physical BL field retains the mapped value; a geometry repair
  // must not silently thin theta/delta or change the prescribed edge speed.
  const packed = prepared.system.initial.slice(prepared.system.ne);
  for (let i = 0; i < packed.length; i++) if (i % 4 !== 0)
    assert.equal(packed[i], options.initialBL[i], `Mapped physical BL slot ${i} changed.`);
  assert.equal(digest({ input, options }), before);
});

test('one unresolved panel transition retains other elements and restores both fallback surfaces and their wake', () => {
  const packet = JSON.parse(gunzipSync(fs.readFileSync(
    new URL('./fixtures/flap-panel-profile-transfer.json.gz', import.meta.url))));
  const beforePacket = structuredClone(packet), f = packet.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options,
    initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
  const original = system.initial.slice(), phase = system.bl.snapshotActive();
  const guess = panelCoupledEdgeGuess(packet.input, f.input.bodies, system.euler.conditions.lengthScale);
  assert.equal(guess.accepted, true);
  const mapped = guess.initialBoundaryLayer(system);
  assert.deepEqual(mapped.profileTransfer.mappedBodies, [0]);
  assert.deepEqual(mapped.profileTransfer.fallbackBodies, [1]);
  assert.deepEqual(mapped.profileTransfer.rejectedSurfaces,
    [{ body: 1, side: 'upper', reason: 'unresolved transferred transition' }]);
  let mappedChanges = 0, surfaceCount = 0, wakeCount = 0;
  for (const station of system.bl.stations) {
    const k = 4 * station.id, prior = original.slice(system.ne + k, system.ne + k + 4);
    if (station.body === 1) {
      assert.deepEqual(mapped.initialBL.slice(k, k + 4), prior, `fallback station ${station.id}`);
      if (station.kind === 'surface') surfaceCount++; else wakeCount++;
    } else if (mapped.initialBL[k + 1] !== prior[1]) mappedChanges++;
  }
  assert(mappedChanges > 10 && surfaceCount > 10 && wakeCount > 10);
  system.bl.surfaces.forEach((surface, k) => {
    if (surface.body === 1) assert.equal(mapped.transitionState[k], phase[k]);
  });
  // Geometric thinning must be undone exactly once, without scaling the
  // finite-base gap or leaving a thin prefix attached to the mapped body.
  const factor = .125, thin = original.slice();
  thin.set(scaleStreamtubeBLThicknesses(system.bl, original.slice(system.ne), factor), system.ne);
  const restored = guess.initialBoundaryLayer({ ...system, initial: thin }, { thicknessFactor: factor });
  assert.deepEqual(restored.transitionState, mapped.transitionState);
  mapped.initialBL.forEach((v, i) => assert(Math.abs(v - restored.initialBL[i]) < 1e-10 * Math.max(1, Math.abs(v)),
    `transferred slot ${i} depends on discarded thinning`));
  assert.deepEqual(system.initial, original);
  assert.deepEqual(system.bl.snapshotActive(), phase);
  assert.deepEqual(packet, beforePacket);
});
