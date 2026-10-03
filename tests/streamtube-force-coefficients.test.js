import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeForceCoefficients } from '../src/euler/streamtube-forces.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('last-iterate coefficients equal independent wall-traction integration at the requested chord and moment origin', () => {
  const system = createStreamtubeBodySystem(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }));
  const flow = solveStreamtubeBody(system, { maxIterations: 0 });
  assert.equal(flow.converged, false);
  for (const referenceChord of [.7, 2.3]) for (const momentReference of [{ x: .25, y: 0 }, { x: -.4, y: .8 }]) {
    const r = streamtubeForceCoefficients(flow, system.conditions, { referenceChord, momentReference });
    // Direct pressure force and nose-up moment about the requested point:
    // integrate each constant-pressure half edge at its own centroid.
    let fx = 0, fy = 0, mz = 0;
    for (const surface of flow.surfaces) for (let k = 1; k < surface.points.length; k++) {
      const a = surface.points[k - 1], b = surface.points[k], sign = surface.side === 'upper' ? -1 : 1;
      for (const [cp, fraction] of [[a.cp, .25], [b.cp, .75]]) {
        const xforce = -sign * cp * (b.y - a.y) / 2, yforce = sign * cp * (b.x - a.x) / 2;
        fx += xforce; fy += yforce;
        mz -= (a.x + fraction * (b.x - a.x) - momentReference.x) * yforce
          - (a.y + fraction * (b.y - a.y) - momentReference.y) * xforce;
      }
    }
    const angle = system.conditions.alpha * Math.PI / 180;
    const expected = { cx: fx / referenceChord, cy: fy / referenceChord,
      cl: (fy * Math.cos(angle) - fx * Math.sin(angle)) / referenceChord,
      cd: (fx * Math.cos(angle) + fy * Math.sin(angle)) / referenceChord, cm: mz / referenceChord ** 2 };
    for (const k of Object.keys(expected)) {
      assert.ok(Math.abs(r[k] - expected[k]) < 2e-13, k);
      assert.equal(r[k], r.perBody.reduce((s, f) => s + f[k], 0));
    }
    assert.equal(r.dragKind, 'pressure'); assert.equal(r.physicalValidation, false);
  }
  assert.throws(() => streamtubeForceCoefficients(flow, system.conditions, { referenceChord: 0 }), /reference/);
  assert.throws(() => streamtubeForceCoefficients({ ...flow, displacement: {} }, system.conditions), /solid-surface/);
  assert.throws(() => streamtubeForceCoefficients({ diagnosticForces: [{ cx: NaN, cy: 0, cm: 0 }] }, system.conditions), /finite/);
});
