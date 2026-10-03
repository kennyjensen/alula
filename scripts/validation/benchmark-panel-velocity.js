// SPDX-License-Identifier: GPL-2.0-or-later
// node scripts/validation/benchmark-panel-velocity.js
// The retained analytic basis sum is an independent, allocation-heavy oracle.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { solveInviscid, velocityAt, vortexBasis } from '../../src/inviscid/linear-vortex.js';
import { sourceVelocity } from '../../src/inviscid/panel.js';
const { input } = JSON.parse(fs.readFileSync(new URL('../../docs/solver-reliability/nlr64x24-subsonic/browser.json', import.meta.url)));
const { field } = solveInviscid({ ...input, mach: 0 });
function basisSum(point, field) {
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
const points = Array.from({ length: 4000 }, (_, i) => ({ x: -1 + 4 * i / 4000, y: .2 + .13 * Math.sin(i) }));
for (const p of points) assert.deepEqual(velocityAt(p, field), basisSum(p, field));
const results = [];
for (const method of ['basis-sum', 'optimized', 'optimized', 'basis-sum']) {
  const start = performance.now(); let sum = 0;
  for (let k = 0; k < 3; k++) for (const p of points) { const v = (method === 'basis-sum' ? basisSum : velocityAt)(p, field); sum += v.u + v.v; }
  results.push({ method, milliseconds: performance.now() - start, sum });
}
console.log(JSON.stringify({ panels: field.panels.length, evaluationsPerSample: 12000, identical: true, results }, null, 2));
