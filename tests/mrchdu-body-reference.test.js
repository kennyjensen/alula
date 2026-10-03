import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { relaxXfoilBody } from '../scripts/validation/mrchdu-body.js';

test('complete MRCHDU preparation matches native surfaces, transition crossings and merged wakes', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-body.json', import.meta.url)));
  for (const [path, hash] of Object.entries(fixture.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(path)).digest('hex'), hash);
  let stations = 0, changedIntervals = 0, maximumError = 0, maximumTEDefect = 0, maximumDefaultChange = 0, defaultBodies = 0;
  for (const c of fixture.cases) {
    const input = structuredClone(c.input), actual = relaxXfoilBody(input, c.parameters);
    assert.deepEqual(input, c.input); assert.equal(actual.flowSolved, false);
    assert.deepEqual(c.nativeWarnings, []); assert.deepEqual(actual.localConvergenceWarnings, []);
    actual.surfaces.forEach((s, side) => {
      assert.equal(s.transition, c.native.surfaces[side].transition);
      assert.equal(s.forced, c.native.surfaces[side].forced);
      assert.ok(Math.abs(s.s - c.native.surfaces[side].s) < 1e-12);
      if (s.transition !== input.phases[side]) changedIntervals++;
    });
    const profiles = [...actual.surfaces.map(s => s.states), actual.wake];
    const expected = [...c.native.surfaces.map(s => s.states), c.native.wake];
    profiles.forEach((p, region) => {
      assert.equal(p.length, expected[region].length);
      p.forEach((a, i) => {
        const b = expected[region][i]; assert.equal(a.s, b.s);
        for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
          const error = Math.abs(a[key] - b[key]) / (key === 'aux' ? Math.max(.01, Math.abs(b[key])) : Math.abs(b[key]));
          maximumError = Math.max(maximumError, error);
          assert.ok(error < 1e-10, `${c.name}/${region}/${i}/${key}: ${error}`);
        }
        stations++;
      });
    });
    const [a, b] = actual.surfaces.map(s => s.states.at(-1)), w = actual.wake[0];
    for (const defect of [w.theta - a.theta - b.theta, w.deltaStar - a.deltaStar - b.deltaStar,
      w.aux * w.theta - a.aux * a.theta - b.aux * b.theta]) {
      maximumTEDefect = Math.max(maximumTEDefect, Math.abs(defect)); assert.ok(Math.abs(defect) < 1e-12);
    }
    if (c.name.startsWith('default-two-element-root')) {
      defaultBodies++;
      const before = [...input.surfaces, input.wake];
      profiles.forEach((p, region) => p.forEach((a, i) => {
        for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
          const b = before[region][i][key], change = Math.abs(a[key] - b) / (key === 'aux' ? Math.max(.01, Math.abs(b)) : Math.abs(b));
          maximumDefaultChange = Math.max(maximumDefaultChange, change); assert.ok(change < 1e-9);
        }
      }));
    }
  }
  assert.equal(stations, 2590); assert.equal(changedIntervals, 2); assert.equal(defaultBodies, 2);
  t.diagnostic(JSON.stringify({ stations, changedIntervals, defaultBodies, maximumError, maximumTEDefect, maximumDefaultChange }));
});
