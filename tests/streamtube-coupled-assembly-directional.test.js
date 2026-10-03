// SPDX-License-Identifier: GPL-2.0-or-later
// Independent whole-residual differences on a supplied349-state assembly.
// No mesh initialization, precursor, LU or Newton solve.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const fixture = JSON.parse(fs.readFileSync(new URL('../docs/coupled-current-profile-preparation/two-element-six-update/initial.json', import.meta.url))).checkpoint.restart;
const cases = [
  ['momentum pressure', 'momentum', 'pressure'],
  ['momentum section speed', 'momentum', 'section-velocity'],
  ['isentropic distance-weighted speed', 'isentropic', 'section-velocity-distance'],
  ...[1, 2, 3, 4].map(ismom => [`ISMOM${ismom} section speed`, 'hybrid', 'section-velocity', ismom]),
];
for (const [label, mode, edgeMatching, ismom] of cases) test(`${label}: tiny coupled Euler, BL and wake Jacobian products match independent residual differences`, t => {
  const f = structuredClone(fixture); f.input.streamwiseMode = mode; f.options.edgeMatching = edgeMatching;
  if (mode === 'hybrid') {
    f.input.hybrid = { ismom, epsilonP: 1e-5 };
    f.input.upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
    f.options.blThermodynamics = 'historical-common-isentrope';
  }
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
  assert.equal(system.n, 349); assert.equal(system.bl.surfaces.length, 4); assert.equal(system.bl.wakes.length, 2);
  const x = system.initial, original = x.slice(), phase = system.bl.snapshotActive(), value = system.evaluate(x), matrix = system.jacobian(x);
  const wake = new Set(system.bl.wakes.flatMap(w => w.ids.flatMap(id => [0, 1, 2, 3].map(k => system.ne + 4 * id + k))));
  const errors = [], base = x.map((v, k) => k < system.ne ? 1e-3 * Math.sin(.73 * k + .2)
    : .02 * Math.max(.01, Math.abs(v)) * Math.sin(.43 * k + .4));
  for (const group of ['euler', 'bl', 'wake', 'mixed']) {
    const direction = base.map((v, k) => group === 'euler' ? k < system.ne ? v : 0
      : group === 'bl' ? k >= system.ne ? v : 0 : group === 'wake' ? wake.has(k) ? v : 0 : v);
    const h = 2e-5, exact = sparseProduct(matrix, direction);
    const sampled = [-2, -1, 1, 2].map(m => system.residual(x.map((v, k) => v + m * h * direction[k])));
    let maximum = 0, worstRow = -1;
    for (let row = 0; row < system.n; row++) {
      const numerical = (sampled[0][row] - 8 * sampled[1][row] + 8 * sampled[2][row] - sampled[3][row]) / (12 * h);
      const error = Math.abs(numerical - exact[row]) / Math.max(1, Math.abs(numerical), Math.abs(exact[row]));
      if (error > maximum) { maximum = error; worstRow = row; }
    }
    assert.ok(maximum < 5e-6, `${group}: row${worstRow} relative error${maximum}`);
    if (group === 'bl' || group === 'wake') assert.ok(exact.subarray(0, system.ne).some(v => Math.abs(v) > 1e-10), `${group} must couple back to Euler`);
    errors.push({ group, maximum, worstRow });
  }
  assert.deepEqual(x, original); assert.deepEqual(system.bl.snapshotActive(), phase); assert.deepEqual(system.residual(x), value.residual);
  t.diagnostic(JSON.stringify({ unknowns: system.n, errors }));
});
