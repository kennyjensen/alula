// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const archive = '../docs/surface-pchip-refinement/mach-continuation/before/';
const originalModuleURL = new URL('../src/euler/streamtube-body.js', import.meta.url);
const originalSource = fs.readFileSync(new URL(`${archive}streamtube-body.js.boundary-diagnostics.txt`, import.meta.url), 'utf8')
  .replace(/from '(\.[^']+)'/g, (_, name) => `from '${new URL(name, originalModuleURL).href}'`);
const previous = await import(`data:text/javascript;base64,${Buffer.from(originalSource).toString('base64')}`);
const fixture = () => ({ ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2, mach: .2 }),
  streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5, ismom: 3 },
  upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } });
const close = (a, b) => assert.ok(Math.abs(a - b) < 3e-13 * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

test('real inlet and outlet guards preserve their message and report the actual rejected gas state', () => {
  const system = createStreamtubeBodySystem(fixture()), original = system.initial.slice();
  const flow = system.evaluate(original), { nx } = system.layout, h0 = system.conditions.h0, gamma = system.conditions.gamma;
  for (const i of [0, nx - 1]) {
    const group = 0, tube = 0, targetMachSquared = 1.21;
    const a = flow.nodes[group][i][tube], b = flow.nodes[group][i + 1][tube];
    const c = flow.nodes[group][i][tube + 1], d = flow.nodes[group][i + 1][tube + 1];
    const sx = .5 * (b.x - a.x + d.x - c.x), sy = .5 * (b.y - a.y + d.y - c.y);
    const wx = .5 * (c.x - a.x + d.x - b.x), wy = .5 * (c.y - a.y + d.y - b.y);
    const area = (sx * wy - sy * wx) / Math.hypot(sx, sy), mass = flow.allocation.groups[group][tube].massFlow;
    const targetQ = Math.sqrt(targetMachSquared * (gamma - 1) * h0 / (1 + .5 * (gamma - 1) * targetMachSquared));
    const state = original.slice(); state[system.layout.densityIndex(i, group, tube)] = Math.log(mass / (area * targetQ));
    let error;
    try { system.evaluate(state); } catch (value) { error = value; }
    assert.ok(error); assert.equal(error.message, `Body inlet/outlet section must remain subsonic: i=${i}, group=0, tube=0.`);
    assert.equal(error.code, 'streamtube-boundary-subsonic');
    const diagnostic = error.diagnostics, rho = Math.exp(state[system.layout.densityIndex(i, group, tube)]);
    const q = mass / (rho * area), enthalpy = h0 - .5 * q * q;
    const machSquared = q * q / ((gamma - 1) * enthalpy);
    assert.deepEqual({ i: diagnostic.i, group: diagnostic.group, tube: diagnostic.tube, boundary: diagnostic.boundary },
      { i, group, tube, boundary: i === 0 ? 'inlet' : 'outlet' });
    close(diagnostic.rho, rho); close(diagnostic.q, q); close(diagnostic.enthalpy, enthalpy);
    close(diagnostic.p, (gamma - 1) / gamma * rho * enthalpy);
    close(diagnostic.machSquared, machSquared); close(diagnostic.machSquared, targetMachSquared);
    assert.equal(diagnostic.machSquaredUpperBound, 1);
    assert.equal(diagnostic.subsonicMargin, 1 - diagnostic.machSquared); assert.ok(diagnostic.subsonicMargin < 0);
    assert.deepEqual(system.initial, original);
  }
});

test('ordinary tiny Euler evaluation matches the archived success path exactly', () => {
  const input = fixture(), current = createStreamtubeBodySystem(input), old = previous.createStreamtubeBodySystem(input);
  assert.deepEqual(current.initial, old.initial);
  const a = current.evaluate(current.initial), b = old.evaluate(old.initial);
  for (const key of ['residual', 'nodes', 'undisplacedNodes', 'sections', 'cells', 'allocation', 'captured',
    'stagnation', 'transportSpeeds', 'diagnostics']) assert.deepEqual(a[key], b[key], key);
});
