// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { finiteBaseBodyFixture } from './fixtures/finite-base-body.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { recoverFiniteBaseWakeTangency, finiteBaseWakeRecoveryPlan, recoverFiniteBaseWakeInterior } from '../src/euler/streamtube-finite-base-wake-initializer.js';
import { pointInside } from '../src/geometry/airfoil.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { streamtubeWakeGap } from '../src/euler/streamtube-wake-geometry.js';
import { refineStreamtubeBody } from '../src/euler/streamtube-refinement.js';

const foldedWake = () => {
  // Exact wake blocks from the reported finite-TE NACA 0012 64x24,
  // alpha 4, Mach .2, automatic spacing, 32 inlet/outlet intervals.
  const f = JSON.parse(fs.readFileSync(new URL('./fixtures/naca64x24-folded-wake.json', import.meta.url)));
  return { ...f, system: { inviscidBaseWake: true, layout: { tubes: [24, 24] } } };
};

test('finite-base downstream bank crossings retain the TE crossline and physical wake boundaries', () => {
  const f = foldedWake(), before = structuredClone(f), touched = [];
  const r = recoverFiniteBaseWakeInterior(f, { admissibleNode: p => { touched.push(p); return true; } });
  assert.equal(r.report.originalQuality.invalidCells.length, 10);
  assert.equal(r.report.accepted, true);
  assert.ok(r.report.quality.valid);
  assert.ok(directStreamtubeVolumeGeometry(r.nodes).valid);
  assert.deepEqual(r.report.groups.map(g => g.group), [0, 1]);
  for (let g = 0; g < 2; g++) {
    assert.deepEqual(r.nodes[g][0], f.nodes[g][0]);
    for (let i = 0; i < f.nodes[g].length; i++) {
      assert.deepEqual(r.nodes[g][i][0], f.nodes[g][i][0]);
      assert.deepEqual(r.nodes[g][i].at(-1), f.nodes[g][i].at(-1));
    }
  }
  assert.notDeepEqual(r.nodes[0].at(-1)[23], f.nodes[0].at(-1)[23]);
  assert.equal(touched.length, 2 * (f.nodes[0].length - 1) * 23);
  assert.deepEqual(f, before);
  assert.equal(recoverFiniteBaseWakeInterior({ ...f, nodes: r.nodes }), null);
});

test('wake interior reconstruction preserves failure guards and does not repair unrelated folds', () => {
  const f = foldedWake(), before = structuredClone(f);
  assert.throws(() => recoverFiniteBaseWakeInterior(f), /original solid-region guard/);
  const rejected = recoverFiniteBaseWakeInterior(f, { admissibleNode: () => false });
  assert.equal(rejected.report.accepted, false);
  assert.equal(rejected.nodes, f.nodes);
  for (const patch of [{ displacement: {} }, { wakeGeometry: 'independent-banks' }])
    assert.equal(recoverFiniteBaseWakeInterior({ ...f, input: { ...f.input, ...patch } }), null);
  assert.equal(recoverFiniteBaseWakeInterior({ ...f, system: { ...f.system, inviscidBaseWake: false } }), null);
  const damaged = structuredClone(f.nodes);
  damaged[0][1][12] = { ...damaged[0][0][12] };
  assert.equal(recoverFiniteBaseWakeInterior({ ...f, nodes: damaged }), null);
  assert.deepEqual(f, before);
});

function fixture() {
  const input = finiteBaseBodyFixture({ bodySegments: 4, tubes: 2 }), b = input.bodies[0], te = b.trailingIndex;
  b.points[0].x += .002; b.points.at(-1).x = b.points[0].x;
  for (const line of [input.outerLower, input.outerUpper, ...input.cutPaths]) { line[te].x = 1.001; line[te + 1].x = 1.00102; }
  const system = createStreamtubeBodySystem(input), nodes = system.decode(system.initial).nodes;
  // Interior stations advance; the purely normal bank alone loses the
  // upper solid TE's tangential component and reverses its first edge.
  nodes[1][te + 1][1].x = nodes[1][te][1].x + .00002;
  return { input, system, nodes, admissibleNode: p => !pointInside(p, b.points) };
}

test('a typed first-finite-base-wake fold retains normal width, solid endpoints and an exact restorable chart', () => {
  const f = fixture(), before = { input: structuredClone(f.input), nodes: structuredClone(f.nodes), chart: f.system.geometryChart() };
  const result = recoverFiniteBaseWakeTangency(f, { admissibleNode: f.admissibleNode });
  assert.deepEqual(result.report.bodies, [0]); assert.equal(result.input.bodies[0].wakeTangentialReference, 'material-te');
  assert.ok(result.report.quality.valid && result.report.quality.minCornerSine > .05);
  assert.equal(result.report.replayError, 0);
  assert.deepEqual(f.input, before.input); assert.deepEqual(f.nodes, before.nodes); assert.deepEqual(f.system.geometryChart(), before.chart);
  const geometry = directStreamtubeVolumeGeometry(result.nodes);
  assert.ok(geometry.valid); assert.equal(geometry.concavePrimal.length, 0);
  const te = f.input.bodies[0].trailingIndex;
  for (const [g, j] of [[0, 2], [1, 0]]) assert.deepEqual(result.nodes[g][te][j], f.nodes[g][te][j]);
  for (let i = te + 1; i < result.system.layout.nx; i++) {
    const bank = (g, j) => [i - 1, i, i + 1].map(k => result.nodes[g][k][j]);
    assert.ok(Math.abs(streamtubeWakeGap(bank(0, 2), bank(1, 0)).gap - result.system.displacement.wakes[0][i - te - 1]) < 1e-14);
    for (const k of ['x', 'y']) assert.equal(.5 * (result.nodes[0][i][2][k] + result.nodes[1][i][0][k]),
      .5 * (f.nodes[0][i][2][k] + f.nodes[1][i][0][k]));
  }
  const restored = createStreamtubeBodySystem(JSON.parse(JSON.stringify(result.input)));
  assert.deepEqual(restored.decode(restored.adoptGeometry(restored.initial, result.nodes)).nodes, result.nodes);
  const refined = refineStreamtubeBody(result.input, result.system, { initial: result.initial, streamwiseFactor: 1,
    normalSubdivisions: [[2, 1], [1, 2]], initializeFlow: false });
  assert.equal(refined.input.bodies[0].wakeTangentialReference, 'material-te');
  assert.ok(directStreamtubeVolumeGeometry(refined.initialEuler.nodes).valid);
});

test('healthy meshes, other folds, independent wakes and prescribed BL data do not activate this initializer', () => {
  const f = fixture(), input = finiteBaseBodyFixture({ bodySegments: 4, tubes: 2 }), healthy = createStreamtubeBodySystem(input);
  assert.equal(recoverFiniteBaseWakeTangency({ input, system: healthy, nodes: healthy.decode(healthy.initial).nodes }), null);
  const damaged = structuredClone(f.nodes); damaged[0][1][1] = { ...damaged[0][0][1] };
  assert.equal(finiteBaseWakeRecoveryPlan(f.input, f.system, damaged), null);
  assert.equal(finiteBaseWakeRecoveryPlan({ ...f.input, wakeGeometry: 'independent-banks' }, f.system, f.nodes), null);
  assert.equal(finiteBaseWakeRecoveryPlan({ ...f.input, displacement: {} }, f.system, f.nodes), null);
  const marked = structuredClone(f.input); marked.bodies[0].wakeTangentialReference = 'material-te';
  assert.equal(finiteBaseWakeRecoveryPlan(marked, f.system, f.nodes), null);
  assert.throws(() => recoverFiniteBaseWakeTangency(f), /original solid-region guard/);
  assert.throws(() => recoverFiniteBaseWakeTangency(f, { admissibleNode: () => false }), e => e.code === 'finite-base-wake-tangency-rejected');
});

for (const [name, mixed] of [['30p-mixed-wake-folds', true], ['three-wake-interior-folds', false]])
  test(`${name}: wake reconstruction composes with tangency and preserves neighboring solid banks`, async () => {
    const { gunzipSync } = await import('node:zlib');
    const f = JSON.parse(gunzipSync(fs.readFileSync(new URL(`./fixtures/${name}.json.gz`, import.meta.url))));
    const before = structuredClone(f), system = createStreamtubeBodySystem(f.input);
    const admissibleNode = p => !f.input.bodies.some(b => pointInside(p, b.points));
    const r = (mixed ? recoverFiniteBaseWakeTangency : recoverFiniteBaseWakeInterior)(
      { ...f, system }, { admissibleNode });
    assert.equal(r.report.accepted, true);
    assert.equal(r.report.quality.valid, true);
    assert.equal(directStreamtubeVolumeGeometry(r.nodes).valid, true);
    if (mixed) assert.equal(r.report.interiorRecovery.accepted, true);
    for (let g = 0; g < f.nodes.length; g++) for (let i = 0; i < f.nodes[g].length; i++)
      for (let j = 0; j < f.nodes[g][i].length; j++) {
        if (system.layout.nodes[g][i][j].kind === 'wall') assert.deepEqual(r.nodes[g][i][j], f.nodes[g][i][j]);
        else assert.ok(admissibleNode(r.nodes[g][i][j]));
      }
    assert.deepEqual(f, before);
  });
