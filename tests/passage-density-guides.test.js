import test from 'node:test';
import assert from 'node:assert/strict';
import { matchPassageDensityGuides } from '../src/euler/passage-density-guides.js';
import { createPanelPotentialBlocks } from '../src/geometry/panel-potential-blocks.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { createContourArc } from '../src/geometry/contour-arc.js';
import { naca4 } from '../src/geometry/airfoil.js';

test('passage placement preserves exact wall geometry, shared cuts, circulation and fixed outer/hit samples', () => {
  const curve = createContourCurve(naca4('0012', 40)), arc = createContourArc(curve), stag = curve.length / 2, origin = arc.at(stag);
  const phase = s => 1 + Math.abs(s - stag) / stag + (s < stag ? .2 : 0) * Math.abs(s - stag) / stag;
  const profile = { curve, stag, phase, upstream: 'upstream', wake: 'wake' };
  const blocks = createPanelPotentialBlocks({ bodies: [{ leading: 1, trailing: { upper: 2.2, lower: 2 }, inlet: 0, outletIncrement: 1 }],
    outer: [{ inlet: 0, outlet: 3.5 }, { inlet: 0, outlet: 3.5 }] });
  const N = 32, ranks = Array.from({ length: 3 * N + 1 }, (_, i) => i / N);
  const le = curve.branch('upper', 0, stag).point, te = curve.branch('upper', 1, stag).point;
  const input = { bodies: [{ leadingIndex: N, trailingIndex: 2 * N, surfaceFractions: { upper: [], lower: [] } }] };
  const guides = [{ upper: [], lower: [] }], descriptors = {};
  for (const side of ['upper', 'lower']) {
    const at = f => (side === 'upper' ? -1 : 1) * (arc.at(curve.branch(side, f, stag).parameter) - origin);
    descriptors[side] = { at, rows: [{ rank: 1, value: 0 }, { rank: 2, value: at(1) }], fractionAtPosition: s => {
      let lo = 0, hi = 1; for (let k = 0; k < 56; k++) { const mid = .5 * (lo + hi); if (at(mid) < s) lo = mid; else hi = mid; } return .5 * (lo + hi);
    } };
    for (let i = 0; i < ranks.length; i++) {
      if (i < N) guides[0][side].push({ x: le.x - 1 + i / N, y: le.y, potential: i / N });
      else if (i > 2 * N) guides[0][side].push({ x: te.x + (i - 2 * N) / N, y: te.y, potential: (side === 'upper' ? 2.2 : 2) + (i - 2 * N) / N });
      else { const f = (i - N) / N, v = curve.branch(side, f, stag); input.bodies[0].surfaceFractions[side].push(f); guides[0][side].push({ ...v.point, potential: phase(v.parameter) }); }
    }
  }
  const outer = [-2, 2].map(y => ranks.map(r => ({ x: r - 1, y, potential: r }))), before = structuredClone({ guides, outer });
  const demand = Array.from({ length: 17 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / 16)));
  const fitted = matchPassageDensityGuides({ input, profiles: [profile], guides, outer, blocks, ranks, surfaceMaps: [descriptors],
    resolvedDemands: [{ upper: demand, lower: demand }], sampleCutX: (path, x) => ({ x, y: path === 'upstream' ? le.y : te.y,
      potential: path === 'upstream' ? x - le.x + 1 : x - te.x }) });
  assert.deepEqual(outer, before.outer);
  assert.notDeepEqual(guides, before.guides);
  for (const side of ['upper', 'lower']) {
    for (const i of [0, N, 2 * N, 3 * N]) assert.deepEqual(guides[0][side][i], before.guides[0][side][i]);
    for (let i = N; i <= 2 * N; i++) {
      const v = curve.branch(side, input.bodies[0].surfaceFractions[side][i - N], stag), p = guides[0][side][i];
      assert.equal(p.x, v.point.x); assert.equal(p.y, v.point.y); assert.equal(p.potential, phase(v.parameter));
    }
  }
  for (let i = 0; i < ranks.length; i++) if (i <= N || i >= 2 * N) {
    const a = guides[0].lower[i], b = guides[0].upper[i]; assert.equal(a.x, b.x); assert.equal(a.y, b.y);
    assert.ok(Math.abs(b.potential - a.potential - (i >= 2 * N ? .2 : 0)) < 1e-14);
  }
  assert.equal(fitted.physicalAcceptance, false);
  assert.equal(fitted.achievedDensity.length, 4);
  assert.ok(fitted.fitted.every(f => f.projection.primalResidual <= 1e-10));
});
