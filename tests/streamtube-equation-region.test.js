// SPDX-License-Identifier: GPL-2.0-or-later
// Small synthetic checks for a fixed, explicitly inherited equation region.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { normalizeStreamtubeEquationRegions, normalizeStreamtubeEquationSelection,
  validateStreamtubeEquationRegionTopology, streamtubeEquationAt } from '../src/euler/streamtube-equation-selection.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const topology = { nx: 24, tubes: [6, 3, 7], bodies: [{ leadingIndex: 2, trailingIndex: 9 }, { leadingIndex: 8, trailingIndex: 15 }] };
const region = () => ({ version: 1, parentCheckpointSha256: 'a'.repeat(64), topology: structuredClone(topology),
  regions: [{ body: 0, throughRow: 13, lowerTube: 2, upperTube: 10 }, { body: 1, throughRow: 19, lowerTube: 5, upperTube: 13 }] });

test('inherited region is immutable and explicit; incompatible controls and malformed provenance/bounds reject', () => {
  const raw = region(), original = structuredClone(raw), value = normalizeStreamtubeEquationRegions(raw);
  assert.deepEqual(value, original);
  for (const object of [value, value.topology, value.topology.tubes, value.topology.bodies,
    ...value.topology.bodies, value.regions, ...value.regions]) assert.ok(Object.isFrozen(object));
  raw.regions[0].throughRow = 14; raw.topology.tubes[0] = 30;
  assert.deepEqual(value, original);
  for (const ismom of [undefined, 1, 2, 4])
    assert.throws(() => normalizeStreamtubeEquationSelection({ epsilonP: .001, ismom, entropyRegions: region() }), /ISMOM3/);
  assert.deepEqual(normalizeStreamtubeEquationSelection({ epsilonP: .001, ismom: 3, entropyRegions: region() }).entropyRegions, original);
  for (const mutate of [r => r.version = 2, r => r.parentCheckpointSha256 = 'unbound', r => r.extra = true,
    r => r.topology.tubes[0] = 0, r => r.topology.bodies[0].trailingIndex = 24,
    r => r.regions.pop(), r => r.regions[0].body = 1, r => r.regions[0].throughRow = 24,
    r => r.regions[0].throughRow = 1, r => r.regions[0].lowerTube = 6,
    r => r.regions[0].lowerTube = -.5, r => r.regions[1].upperTube = 9,
    r => r.regions[1].upperTube = 17, r => r.regions[1].upperTube = Infinity,
    r => delete r.topology.tubes[0], r => delete r.topology.bodies[0], r => delete r.regions[0]]) {
    const bad = region(); mutate(bad); assert.throws(() => normalizeStreamtubeEquationRegions(bad), /ISMOM3/);
  }
});

test('fixed inherited multielement union uses inclusive row and half-open global tube bounds with strict grid binding', () => {
  const entropyRegions = normalizeStreamtubeEquationRegions(region()), hybrid = { ismom: 3, epsilonP: .001, entropyRegions };
  validateStreamtubeEquationRegionTopology(entropyRegions, topology);
  for (let i = 1; i < topology.nx; i++) for (let global = 0; global < 16; global++) {
    let group = 0, tube = global;
    while (tube >= topology.tubes[group]) tube -= topology.tubes[group++];
    const expected = (i <= 13 && global >= 2 && global < 10) || (i <= 19 && global >= 5 && global < 13);
    assert.equal(streamtubeEquationAt({ hybrid, ...topology, i, group, tube }), expected ? 'isentropic' : 'momentum');
  }
  for (const mutate of [t => t.nx++, t => t.tubes[0]++, t => t.bodies[0].leadingIndex++, t => t.bodies[1].trailingIndex++]) {
    const other = structuredClone(topology); mutate(other);
    assert.throws(() => validateStreamtubeEquationRegionTopology(entropyRegions, other), /topology does not match/);
    assert.throws(() => streamtubeEquationAt({ hybrid, ...other, i: 1, group: 0, tube: 0 }), /topology does not match/);
  }
});

test('inherited region selects the same exact residual/J rows on moving two-element geometry', t => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 3, mach: .5, alpha: .25 }),
    streamwiseMode: 'hybrid', hybrid: { epsilonP: .001, ismom: 1 },
    upwind: { mucon: 1, mcrit: .5, boundary: { kind: 'unfiltered-first-two' } } };
  const a = createStreamtubeBodySystem(input), layout = a.layout;
  const entropyRegions = { version: 1, parentCheckpointSha256: 'b'.repeat(64),
    topology: { nx: layout.nx, tubes: layout.tubes, bodies: layout.bodies.map(({ leadingIndex, trailingIndex }) => ({ leadingIndex, trailingIndex })) },
    regions: layout.bodies.map((body, b) => ({ body: b, throughRow: body.leadingIndex + 1,
      lowerTube: layout.groupOffsets[b + 1] - 1, upperTube: layout.groupOffsets[b + 1] + 1 })) };
  const b = createStreamtubeBodySystem({ ...input, hybrid: { ...input.hybrid, ismom: 2 } });
  const c = createStreamtubeBodySystem({ ...input, hybrid: { ...input.hybrid, ismom: 3, entropyRegions } });
  const x = a.initial.map((_, k) => (k < layout.densityCount ? 1e-3 : 1e-5) * Math.sin(k + .3));
  const values = [a, b, c].map(s => s.evaluate(x)), matrices = [a, b, c].map(s => s.jacobian(x, { sparse: true }));
  const row = (m, r) => Array.from({ length: m.rowPtr[r + 1] - m.rowPtr[r] }, (_, k) => {
    const at = m.rowPtr[r] + k; return [m.colIndex[at], m.values[at]];
  }).filter(([, v]) => v !== 0).sort((p, q) => p[0] - q[0]);
  let entropy = 0, momentum = 0;
  for (const r of layout.rows) {
    const global = layout.groupOffsets[r.group] + r.tube;
    const useEntropy = r.kind === 'streamwise' && entropyRegions.regions.some(e => r.i <= e.throughRow && global >= e.lowerTube && global < e.upperTube);
    const source = useEntropy ? 1 : 0;
    assert.equal(values[2].residual[r.index], values[source].residual[r.index]);
    assert.deepEqual(row(matrices[2], r.index), row(matrices[source], r.index));
    if (r.kind === 'streamwise') useEntropy ? entropy++ : momentum++;
  }
  assert.ok(entropy > 0 && momentum > 0);
  const d = x.map((_, k) => .02 * Math.cos(k + .4)), h = 7e-7;
  const samples = [-2, -1, 1, 2].map(m => c.residual(x.map((v, k) => v + m * h * d[k]))), jv = sparseProduct(matrices[2], d);
  let worst = 0;
  for (let k = 0; k < x.length; k++) {
    const fd = (samples[0][k] - 8 * samples[1][k] + 8 * samples[2][k] - samples[3][k]) / (12 * h);
    worst = Math.max(worst, Math.abs(fd - jv[k]) / Math.max(1, Math.abs(fd), Math.abs(jv[k])));
  }
  assert.ok(worst < 2e-7, String(worst));
  t.diagnostic(JSON.stringify({ unknowns: x.length, entropy, momentum, directionalError: worst }));
  assert.throws(() => createStreamtubeBodySystem({ ...input, hybrid: { ...input.hybrid, ismom: 3,
    entropyRegions: { ...entropyRegions, topology: { ...entropyRegions.topology, nx: layout.nx + 1 } } } }), /topology does not match/);
});
