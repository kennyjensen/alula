// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { refineStreamtubeBody } from '../src/euler/streamtube-refinement.js';
import { streamtubeWakeGap } from '../src/euler/streamtube-wake-geometry.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';
import { finiteBaseBodyFixture } from './fixtures/finite-base-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const counts = [[1, 2, 3], [3, 2, 1]];
const controls = { streamwiseFactor: 1, normalSubdivisions: counts, initializeFlow: false };
const distance = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);

for (const normalInterpolation of ['linear', 'streamfunction-quadratic'])
  test(`finite-base ${normalInterpolation} normal refinement retains physical banks and parent mass`, () => {
    const input = finiteBaseBodyFixture(), originalInput = structuredClone(input);
    const source = createStreamtubeBodySystem(input), original = source.decode(source.initial), chart = source.geometryChart();
    const r = refineStreamtubeBody(input, { ...source, evaluate: () => assert.fail('Preview must not evaluate flow') },
      { ...controls, normalInterpolation });
    const next = r.system.decode(r.initial), te = input.bodies[0].trailingIndex;
    assert.equal(r.system.inviscidBaseWake, true);
    assert.deepEqual(r.system.displacement, source.displacement);
    assert.deepEqual(r.input.bodies[0].points, input.bodies[0].points);
    assert.deepEqual(r.input.bodies[0].trailingEdge, input.bodies[0].trailingEdge);
    assert.deepEqual(next.captured, original.captured);
    assert.equal(r.system.conditions.lengthScale, source.conditions.lengthScale);
    assert.equal(r.system.conditions.massScale, source.conditions.massScale);
    assert.equal(r.diagnostics.finiteBaseWake.retainedStreamwiseStations, true);
    assert.equal(r.diagnostics.flowInitialized, false);
    assert.equal(r.diagnostics.residual, undefined);
    for (let g = 0; g < counts.length; g++) {
      let child = 0;
      for (let j = 0; j < counts[g].length; j++) {
        const sum = next.allocation.groups[g].slice(child, child + counts[g][j]).reduce((s, tube) => s + tube.massFlow, 0);
        assert.ok(Math.abs(sum - original.allocation.groups[g][j].massFlow) < 1e-14);
        for (let i = 0; i <= source.layout.nx; i++)
          assert.ok(distance(next.nodes[g][i][child], original.nodes[g][i][j]) < 2e-14);
        child += counts[g][j];
      }
      for (let i = 0; i <= source.layout.nx; i++)
        assert.ok(distance(next.nodes[g][i].at(-1), original.nodes[g][i].at(-1)) < 2e-14);
    }
    assert.ok(distance(next.nodes[0][te].at(-1), next.nodes[1][te][0]) > 0);
    assert.deepEqual(next.nodes[0][te].at(-1), input.bodies[0].points[input.bodies[0].trailingEdge.lowerIndex]);
    assert.deepEqual(next.nodes[1][te][0], input.bodies[0].points[input.bodies[0].trailingEdge.upperIndex]);
    for (let i = te + 1; i <= source.layout.nx; i++) {
      const ids = [i - 1, i, Math.min(source.layout.nx, i + 1)];
      const gap = streamtubeWakeGap(ids.map(k => next.nodes[0][k].at(-1)), ids.map(k => next.nodes[1][k][0]));
      assert.ok(Math.abs(gap.gap - source.baseGeometry[0].width) < 2e-14);
    }
    assert.equal(r.diagnostics.quality.valid, true);
    assert.equal(directStreamtubeVolumeGeometry(next.nodes).valid, true);
    assert.deepEqual(source.decode(source.initial).nodes, original.nodes);
    assert.deepEqual(source.geometryChart(), chart);
    assert.deepEqual(input, originalInput);
  });

test('finite-base initialized refinement keeps exactly the preview geometry and admissible subsonic gas', () => {
  const input = finiteBaseBodyFixture(), source = createStreamtubeBodySystem(input);
  const preview = refineStreamtubeBody(input, source, controls);
  const initialized = refineStreamtubeBody(input, source, { ...controls, initializeFlow: true });
  assert.deepEqual(initialized.initialEuler.nodes, preview.initialEuler.nodes);
  assert.ok(initialized.diagnostics.maxMach < 1);
  assert.ok(Number.isFinite(initialized.diagnostics.residual));
  assert.equal(initialized.diagnostics.flowInitialized, true);
});

test('finite-base refinement rejects BL displacement, altered gap, streamwise insertion and oversized output', () => {
  const input = finiteBaseBodyFixture(), source = createStreamtubeBodySystem(input);
  assert.throws(() => refineStreamtubeBody(input, source, { ...controls, streamwiseFactor: 2 }), /normal subdivisions only/);
  assert.throws(() => refineStreamtubeBody(input, source, { ...controls, maxNodes: 1 }), /node budget/);
  const original = structuredClone(source.displacement), surface = structuredClone(original);
  surface.surfaces[0].upper[1] = 1e-5; source.setDisplacement(surface);
  assert.throws(() => refineStreamtubeBody(input, source, controls), /unchanged zero-wall, constant-width/);
  const wake = structuredClone(original); wake.wakes[0][0] *= 1.001; source.setDisplacement(wake);
  assert.throws(() => refineStreamtubeBody(input, source, controls), /unchanged zero-wall, constant-width/);
  const explicitInput = { ...input, displacement: original }, explicitSource = createStreamtubeBodySystem(explicitInput);
  assert.equal(explicitSource.inviscidBaseWake, false);
  assert.throws(() => refineStreamtubeBody(explicitInput, explicitSource, controls), /Use coupled refinement/);
});

test('sharp geometry-only refinement remains identical to the archived pre-change implementation', async () => {
  const archive = new URL('../docs/nlr-finite-base/before/refinement/streamtube-refinement.js.txt', import.meta.url);
  const oldText = fs.readFileSync(archive, 'utf8').replace(/from '([^']+)'/g,
    (_, p) => `from '${new URL(p, new URL('../src/euler/streamtube-refinement.js', import.meta.url)).href}'`);
  const oldModule = await import('data:text/javascript;base64,' + Buffer.from(oldText).toString('base64'));
  const input = intrinsicBodyFixture({ bodySegments: 8, tubes: 3 });
  const current = refineStreamtubeBody(input, createStreamtubeBodySystem(input), controls);
  const old = oldModule.refineStreamtubeBody(input, createStreamtubeBodySystem(input), controls);
  for (const key of ['input', 'initial', 'initialEuler', 'diagnostics']) assert.deepEqual(current[key], old[key]);
});
