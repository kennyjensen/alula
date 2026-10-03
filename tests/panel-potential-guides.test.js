import test from 'node:test';
import assert from 'node:assert/strict';
import { matchPanelPotentialGuides } from '../src/euler/streamtube-potential-guides.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { createContourArc } from '../src/geometry/contour-arc.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

test('geometric block hits use facing potentials, stop at walls and transfer wake circulation', () => {
  // Prescribed monotone surface potential phi=phi_LE+(phi_TE-phi_LE)*f².
  // The known intersections below have closed-form surface fractions.
  const values = [{ leading: 4.8, upper: 5.2, lower: 5 }, { leading: 4.1, upper: 5.7, lower: 4.9 }];
  const profiles = values.map((v, b) => {
    const curve = createContourCurve(transform(naca4('0012', 24), { chord: b ? 1 : .3, y: b }));
    const stag = curve.length / 2;
    return { curve, stag, phiStag: v.leading,
      phase: s => v.leading + (v[s < stag ? 'upper' : 'lower'] - v.leading) * ((s - stag) / stag) ** 2,
      upstream: [{ x: -4, y: b, potential: 0 }], wake: [{ x: 6, y: b, potential: 3 }] };
  });
  const input = { bodies: [{}, {}], gridSpacing: { inlet: { intervals: 16 }, outlet: { intervals: 16 } } };
  const sample = (path, potential) => ({ x: potential, y: path[0].y, potential });
  const result = matchPanelPotentialGuides({ input, profiles,
    surfaceFractions: profiles.map(() => ({ upper: [0, .1, .4, .7, 1], lower: [0, .1, .4, .7, 1] })),
    outerPaths: [-1, 2].map(y => [{ x: 0, y, potential: 0 }, { x: 9, y, potential: 9 }]),
    sample, sampleCut: sample, growth: .25, maxIntervals: 1000, outerSpread: .85, surfaceMetric: 'arc' });
  const [flap, main] = input.bodies;
  assert.ok(main.leadingIndex < flap.leadingIndex && flap.leadingIndex < main.trailingIndex && main.trailingIndex < flap.trailingIndex);
  const atMainTE = result.guides[0].upper[main.trailingIndex];
  const expectedFlap = profiles[0].curve.branch('upper', .5, profiles[0].stag).point;
  assert.ok(Math.hypot(atMainTE.x - expectedFlap.x, atMainTE.y - expectedFlap.y) < 1e-13);
  const atFlapLE = result.guides[1].lower[flap.leadingIndex];
  const expectedMain = profiles[1].curve.branch('lower', Math.sqrt(.875), profiles[1].stag).point;
  assert.ok(Math.hypot(atFlapLE.x - expectedMain.x, atFlapLE.y - expectedMain.y) < 1e-13);
  // Crossing the main wake transfers the TE jump .8 into the upper region.
  assert.ok(Math.abs(result.outer[1][flap.trailingIndex].potential - 6) < 1e-13);
  assert.equal(result.outer[0][flap.trailingIndex].potential, 5);
  // Within the actual facing-wall block, the two reconciled arc maps have
  // the same normalized progress, including interior nodes. Shared rank
  // samples alone could not give this with the old independent cubics.
  const leading = flap.leadingIndex, trailing = main.trailingIndex;
  const arcRows = [[0, 'upper'], [1, 'lower']].map(([b, side]) => {
    const arc = createContourArc(profiles[b].curve), body = input.bodies[b];
    const values = result.guides[b][side].slice(leading, trailing + 1).map((_, j) => {
      const f = body.surfaceFractions[side][leading + j - body.leadingIndex];
      return arc.at(profiles[b].curve.branch(side, f, profiles[b].stag).parameter);
    });
    return values.map(s => (s - values[0]) / (values.at(-1) - values[0]));
  });
  arcRows[0].forEach((s, i) => assert.ok(Math.abs(s - arcRows[1][i]) < 1e-12));
  assert.equal(result.diagnostics.facingSpacing.length, 1);
  assert.equal(result.diagnostics.facingSpacing[0].exactModernMsetLaw, false);
  for (const [b, body] of input.bodies.entries()) for (let i = 0; i < result.diagnostics.x.length; i++) {
    if (i > body.leadingIndex && i < body.trailingIndex) continue;
    const a = result.guides[b].lower[i], c = result.guides[b].upper[i];
    assert.equal(a.x, c.x); assert.equal(a.y, c.y);
    assert.ok(Math.abs((c.potential - a.potential) - (i < body.trailingIndex ? 0 : values[b].upper - values[b].lower)) < 1e-13);
  }
  for (const event of result.diagnostics.panelBlocks.events) {
    const i = result.diagnostics.x.indexOf(event.rank);
    for (const [g, target] of event.passagePotentials.entries()) if (target !== null) {
      const a = g === 0 ? result.outer[0][i] : result.guides[g - 1].upper[i];
      const c = g === profiles.length ? result.outer[1][i] : result.guides[g].lower[i];
      assert.ok(Math.abs(a.potential - target) < 1e-13);
      assert.ok(Math.abs(c.potential - target) < 1e-13);
    }
  }
});
