import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { relaxXfoilSurfaces } from '../scripts/validation/mrchdu-surface.js';

test('MRCHDU agrees with original Fortran on refined-grid roots, Newton rays and transition crossings', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/mrchdu-wall.json', import.meta.url)));
  for (const [path, hash] of Object.entries(fixture.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(path)).digest('hex'), hash);
  let maximumError = 0, stations = 0, changedIntervals = 0;
  for (let i = 0; i < fixture.cases.length; i += 2) {
    const pair = fixture.cases.slice(i, i + 2), profiles = pair.map(c => structuredClone(c.profile));
    assert.deepEqual(pair[0].parameters, pair[1].parameters);
    const result = relaxXfoilSurfaces(profiles, pair[0].parameters, pair.map(c => c.phase));
    assert.deepEqual(profiles, pair.map(c => c.profile));
    assert.equal(result.flowSolved, false); assert.deepEqual(result.localConvergenceWarnings, []);
    result.surfaces.forEach((surface, side) => {
      const c = pair[side], expected = c.native;
      assert.equal(surface.transition, expected.transition, c.name);
      assert.equal(surface.forced, expected.forced, c.name);
      assert.ok(Math.abs(surface.s - expected.s) < 1e-12);
      if (surface.transition !== c.phase) changedIntervals++;
      surface.states.forEach((p, j) => {
        for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
          const a = p[key], b = expected.states[j][key], scale = key === 'aux' ? Math.max(.01, Math.abs(b)) : Math.abs(b);
          const error = Math.abs(a - b) / scale; maximumError = Math.max(maximumError, error);
          assert.ok(error < 1e-10, `${c.name}/${j}/${key}: ${error}`);
        }
        stations++;
      });
    });
  }
  assert.equal(stations, 1182); assert.equal(changedIntervals, 2);
  t.diagnostic(JSON.stringify({ maximumError, stations, changedIntervals }));
});
