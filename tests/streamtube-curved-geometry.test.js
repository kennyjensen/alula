import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createCurvedStreamtubeRegions } from '../src/euler/streamtube-curved-geometry.js';
import { solveInviscid, velocityAt } from '../src/inviscid/linear-vortex.js';

let cached;
function specimen() {
  if (cached) return cached;
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/streamtube-solved-level-initialization.json', import.meta.url), 'utf8'));
  const input = saved.preparedInput, system = createStreamtubeBodySystem(input);
  // Reconstruct only the global panel field and body metadata. The stored
  // grid nodes are reused directly; no initializer or tracing is invoked.
  const panel = solveInviscid({ elements: input.bodies.map(body => ({ points: body.points })),
    alpha: input.alpha, boundaryCondition: saved.controls.panelBoundaryCondition });
  assert.equal(panel.status, 'solved');
  const prepared = { input, system, nodes: saved.original.nodes,
    guideField: { velocityAt: p => velocityAt(p, panel.field) }, diagnostics: saved.original.diagnostics };
  cached = { prepared, regions: createCurvedStreamtubeRegions(prepared) }; return cached;
}

test('saved default curved builder retains original nodes and exact C2 wall endpoint derivatives', t => {
  const start = performance.now(), { prepared, regions } = specimen();
  assert.equal(regions.length, 3);
  const { system, input, nodes } = prepared, decoded = system.decode(system.initial);
  for (let g = 0; g < regions.length; g++) {
    const region = regions[g]; assert.deepEqual(region.geometry.initial, nodes[g]);
    assert.deepEqual(region.massFlows, decoded.allocation.groups[g].map(t => t.massFlow));
    assert.equal(region.boundaryApproximation.physicsValidated, false);
    assert.ok(region.boundaryApproximation.maximumWallEndpointMismatch < 1e-12);
    for (const [body, side, j] of [[g - 1, 'upper', 0], [g, 'lower', nodes[g][0].length - 1]]) {
      if (body < 0 || body >= input.bodies.length) continue;
      const range = input.bodies[body], fractions = system.fractions[body][side], curve = system.curves[body], stag = decoded.stagnation[body];
      for (let i = range.leadingIndex; i < range.trailingIndex; i++) {
        const k = i - range.leadingIndex, factor = (fractions[k + 1] - fractions[k]) * (side === 'upper' ? -stag : curve.length - stag);
        for (const s of [0, 1]) {
          const expected = curve.branch(side, fractions[k + s], stag).derivative;
          const actual = region.geometry.at(i, j ? j - 1 : 0, s, j ? 1 : 0).ds;
          assert.ok(Math.hypot(actual.x - factor * expected.x, actual.y - factor * expected.y) < 2e-12);
        }
      }
    }
  }
  t.diagnostic(JSON.stringify({ seconds: (performance.now() - start) / 1000,
    boundaryApproximation: regions.map(r => r.boundaryApproximation) }));
});

test('shared dividing and wake curves coincide across neighboring fluid regions', () => {
  const { prepared, regions } = specimen();
  for (let body = 0; body < prepared.input.bodies.length; body++) {
    const range = prepared.input.bodies[body], lowerNt = prepared.nodes[body][0].length - 1;
    for (let i = 0; i < prepared.system.layout.nx; i++) {
      if (i >= range.leadingIndex && i < range.trailingIndex) continue;
      for (const s of [0, .21, .57, 1]) {
        const a = regions[body].geometry.at(i, lowerNt - 1, s, 1);
        const b = regions[body + 1].geometry.at(i, 0, s, 0);
        assert.ok(Math.hypot(a.point.x - b.point.x, a.point.y - b.point.y) < 2e-12);
        assert.ok(Math.hypot(a.ds.x - b.ds.x, a.ds.y - b.ds.y) < 2e-12);
      }
    }
  }
});

test('curved builder preserves field failures and rejects invalid prepared dimensions', () => {
  const { prepared } = specimen();
  assert.throws(() => createCurvedStreamtubeRegions({ ...prepared, nodes: [] }), /prepared body/);
  assert.throws(() => createCurvedStreamtubeRegions({ ...prepared, guideField: { velocityAt: () => ({ u: 0, v: 0 }) } }), /Unresolved/);
  assert.throws(() => createCurvedStreamtubeRegions({ ...prepared, guideField: { velocityAt: () => ({ u: NaN, v: 0 }) } }), /Nonfinite/);
});
