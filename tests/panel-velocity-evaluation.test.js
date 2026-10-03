// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { makePanel, sourceVelocity } from '../src/inviscid/panel.js';
import { velocityAt, vortexBasis } from '../src/inviscid/linear-vortex.js';

function reference(point, field) {
  let { u, v } = field;
  for (const p of field.panels) {
    const basis = vortexBasis(point, p);
    for (let k = 0; k < 2; k++) { u += field.gamma[p.node + k] * basis[k].u; v += field.gamma[p.node + k] * basis[k].v; }
  }
  for (const p of field.basePanels ?? []) {
    const source = sourceVelocity(point, p);
    u += p.sourceStrength * source.u - p.vortexStrength * source.v;
    v += p.sourceStrength * source.v + p.vortexStrength * source.u;
  }
  return { u, v };
}

test('allocation-reduced field agrees bit-for-bit with the analytic basis sum, including finite bases', () => {
  const panels = Array.from({ length: 40 }, (_, i) => ({ ...makePanel(
    { x: Math.cos(i), y: Math.sin(i) }, { x: Math.cos(i + .17), y: Math.sin(i + .17) }), node: i }));
  const field = { panels, gamma: Float64Array.from({ length: 41 }, (_, i) => Math.sin(i * .3)), u: .97, v: .12,
    basePanels: [{ ...makePanel({ x: 1, y: -.003 }, { x: 1, y: .003 }), sourceStrength: .014, vortexStrength: -.008 }] };
  const points = Array.from({ length: 150 }, (_, i) => ({ x: 1.01 * Math.cos(i * .11), y: 1.3 * Math.sin(i * .11) }));
  points.push({ x: -12, y: 30 }, { x: 1, y: 1e-12 }, { x: .99, y: -.003 });
  for (const point of points) assert.deepEqual(velocityAt(point, field), reference(point, field));
  field.gamma[3] += .5;
  for (const point of points) assert.deepEqual(velocityAt(point, field), reference(point, field));
  assert.throws(() => velocityAt(panels[0].a, field), /singular at a panel endpoint/);
});
