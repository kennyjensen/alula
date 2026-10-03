// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { prepareStreamtubeTransportChain } from '../src/euler/streamtube-transport-chain.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
const banks = () => ({ lower: [0, 1, 2, 3, 4].map(x => ({ x, y: 0 })),
  upper: [0, 1, 2, 3, 4].map(x => ({ x, y: .5 })),
  densities: [1, 1, 1, 1], massFlow: .4, stagnationEnthalpy: 5, upwind });

test('crossed transport geometry preserves the original rejection and actual local trial points', () => {
  const input = { ...banks(), geometryDomain: 'positive-simple' }; input.upper[3].y = -.5;
  const before = structuredClone(input); let caught;
  assert.throws(() => prepareStreamtubeTransportChain(input), error => {
    caught = error;
    assert.equal(error.message, 'Crossed, reversed or degenerate quadrilateral polygon.');
    assert.equal(error.code, 'streamtube-transport-geometry');
    assert.deepEqual(error.diagnostics.stencil, { i: 2, bankStations: [1, 2, 3], geometryDomain: 'positive-simple',
      lower: before.lower.slice(1, 4), upper: before.upper.slice(1, 4) });
    return true;
  });
  assert.deepEqual(input, before);
  input.upper[3].y = 7;
  assert.equal(caught.diagnostics.stencil.upper[2].y, -.5);
});

test('real body rejection retains group, tube and the six rejected stencil points', () => {
  const system = createStreamtubeBodySystem({ ...intrinsicBodyFixture({ bodySegments: 4, tubes: 3, mach: .2 }),
    geometryDomain: 'positive-simple', streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5, ismom: 3 }, upwind });
  const original = system.initial.slice(), trial = original.slice();
  const column = system.layout.nodes[0][2][1].column;
  assert.notEqual(column, null); trial[column] += 100;
  const decoded = system.decode(trial);
  assert.throws(() => system.evaluate(trial), error => {
    assert.equal(error.code, 'streamtube-transport-geometry');
    assert.match(error.message, /Body transport chain group=\d+, tube=\d+:/);
    assert.equal(error.cause.code, error.code);
    const { group, tube, stencil } = error.diagnostics;
    assert.ok(Number.isInteger(group) && Number.isInteger(tube));
    assert.deepEqual(stencil.lower, decoded.nodes[group].slice(stencil.i - 1, stencil.i + 2).map(row => row[tube]));
    assert.deepEqual(stencil.upper, decoded.nodes[group].slice(stencil.i - 1, stencil.i + 2).map(row => row[tube + 1]));
    return true;
  });
  assert.deepEqual(system.initial, original);
});

test('successful transport-chain values remain identical to the archived implementation', async () => {
  const originalURL = new URL('../src/euler/streamtube-transport-chain.js', import.meta.url);
  const source = fs.readFileSync(new URL('../docs/transport-stencil-diagnostics/before/streamtube-transport-chain.js.txt', import.meta.url), 'utf8')
    .replace(/from '(\.[^']+)'/g, (_, path) => `from '${new URL(path, originalURL).href}'`);
  const old = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
  const input = banks();
  assert.deepEqual(prepareStreamtubeTransportChain(input), old.prepareStreamtubeTransportChain(input));
});
