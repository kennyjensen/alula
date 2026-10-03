import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { initializeXfoilSurfaces } from '../src/viscous/xfoil-surface-initializer.js';

test('MRCHUE surface adapter agrees with executed unmodified Fortran, including flagged initializer failures', t => {
  const fixture = JSON.parse(fs.readFileSync(new URL('../docs/current-quad-automatic-mrchue-before.json', import.meta.url)));
  for (const [path, hash] of Object.entries(fixture.provenance.sha256))
    assert.equal(createHash('sha256').update(fs.readFileSync(path)).digest('hex'), hash);
  let maximumError = 0, stationCount = 0, warnings = 0;
  for (const c of fixture.cases) {
    const before = structuredClone(c.profile), r = initializeXfoilSurfaces([c.profile, c.profile], c.parameters);
    assert.deepEqual(c.profile, before); assert.equal(r.flowSolved, false);
    warnings += r.localConvergenceWarnings.length;
    for (const side of r.surfaces) {
      assert.equal(side.transition, c.native.transition, c.name); assert.equal(side.forced, c.native.forced, c.name);
      assert.ok(Math.abs(side.s - c.native.s) < 1e-12, c.name);
      side.states.forEach((p, i) => {
        for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
          const a = p[key], b = c.native.states[i][key], scale = key === 'aux' ? Math.max(.01, Math.abs(b)) : Math.abs(b);
          const error = Math.abs(a - b) / scale; maximumError = Math.max(maximumError, error);
          assert.ok(error < 1e-10, `${c.name}/${i}/${key}: ${error}`);
        }
        stationCount++;
      });
    }
  }
  assert.ok(warnings > 0, 'Coarse native initializer failures must remain explicit.');
  t.diagnostic(JSON.stringify({ maximumError, stationCount, warnings }));
});

test('different upper/lower profiles retain their independent native initialization', () => {
  const fixture = JSON.parse(fs.readFileSync(new URL('../docs/current-quad-automatic-mrchue-before.json', import.meta.url)));
  for (const body of [0, 1]) {
    const cases = fixture.cases.filter(c => c.name.startsWith(`default-cold/${body}/`));
    assert.equal(cases.length, 2);
    const r = initializeXfoilSurfaces(cases.map(c => c.profile), cases[0].parameters);
    r.surfaces.forEach((surface, side) => {
      const reference = cases[side].native;
      assert.equal(surface.transition, reference.transition); assert.equal(surface.states.length, reference.states.length);
      surface.states.forEach((p, i) => {
        for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
          const expected = reference.states[i][key];
          assert.ok(Math.abs(p[key] - expected) < 1e-10 * Math.max(Math.abs(expected), 1e-6));
        }
      });
    });
  }
});
