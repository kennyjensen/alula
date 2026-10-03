import test from 'node:test';
import assert from 'node:assert/strict';
import { matchPanelPotentialGuides } from '../src/euler/streamtube-potential-guides.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { createMonotoneCubicMap } from '../src/numerics/monotone-cubic.js';
import { naca4, transform } from '../src/geometry/airfoil.js';

test('full physical outlines reconcile cuts and walls while preserving physical hits, curves and circulation', () => {
  const laws = new Map(), values = [{ leading: 4.8, upper: 5.2, lower: 5 }, { leading: 4.1, upper: 5.7, lower: 4.9 }];
  const profiles = values.map((v, b) => {
    const curve = createContourCurve(transform(naca4('0012', 24), { chord: b ? 1 : .3, x: b ? 0 : .94, y: b }));
    const stag = curve.length / 2, le = curve.branch('upper', 0, stag).point, te = curve.branch('upper', 1, stag).point;
    const upstream = [{ x: le.x - Math.sqrt(2 * v.leading), y: le.y, potential: 0 },
      { x: le.x - .002, y: le.y, potential: v.leading - .5 * .002 ** 2 }];
    const wake = [{ x: te.x + .002, y: te.y, potential: .002 }, { x: te.x + 3, y: te.y, potential: 3 }];
    laws.set(upstream, { x: phi => le.x - Math.sqrt(2 * (v.leading - phi)), phi: x => v.leading - .5 * (x - le.x) ** 2 });
    laws.set(wake, { x: phi => te.x + phi, phi: x => x - te.x });
    return { curve, stag, upstream, wake, phiStag: v.leading,
      phase: s => v.leading + (v[s < stag ? 'upper' : 'lower'] - v.leading) * ((s - stag) / stag) ** 2 };
  });
  const sample = (path, potential) => ({ x: laws.has(path) ? laws.get(path).x(potential) : potential, y: path[0].y, potential });
  let samples = 0;
  const sampleCutX = (path, x) => { samples++; return { x, y: path[0].y, potential: laws.get(path).phi(x) }; };
  const build = (cutStationSpacing, surfaceStationPlacement = 'common-rank', fractions = [0, .1, .4, .7, 1], passageCountPlanning = false) => {
    const input = { bodies: [{}, {}], gridSpacing: { inlet: { intervals: 16 }, outlet: { intervals: 16 } } };
    const result = matchPanelPotentialGuides({ input, profiles, cutStationSpacing, sampleCutX, surfaceStationPlacement, passageCountPlanning,
      surfaceFractions: profiles.map(() => ({ upper: fractions, lower: fractions })),
      outerPaths: [-1, 2].map(y => [{ x: 0, y, potential: 0 }, { x: 9, y, potential: 9 }]),
      sample, sampleCut: sample, growth: .25, maxIntervals: 1000, outerSpread: .85, surfaceMetric: 'arc' });
    return { input, ...result };
  };
  const legacy = build('potential'), physical = build('physical-x');
  const local = build('physical-x', 'local-density');
  const dense = Array.from({ length: 33 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / 32)));
  const localFine = build('physical-x', 'local-density', dense), joined = build('physical-x', 'local-density-joined', dense);
  const automatic = build('physical-x', 'passage-density', dense, true), planned = automatic.diagnostics.passageDensity;
  assert.equal(planned.countPlan.allOriginalResolutionSatisfied, true);
  assert.ok(planned.achievedDensity.every(d => d.referenceOnly || d.resolutionSatisfied));
  assert.ok(automatic.diagnostics.blocks[0].intervals >= localFine.diagnostics.blocks[0].intervals);
  assert.ok(automatic.diagnostics.blocks.at(-1).intervals >= localFine.diagnostics.blocks.at(-1).intervals);
  for (let b = 0; b < profiles.length; b++) for (const side of ['upper', 'lower']) {
    const body = automatic.input.bodies[b];
    for (const hit of automatic.diagnostics.panelBlocks.constraints[side][b]) {
      const p = automatic.guides[b][side][automatic.diagnostics.x.indexOf(hit.rank)], q = localFine.guides[b][side][localFine.diagnostics.x.indexOf(hit.rank)];
      assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 1e-14);
    }
    for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
      const p = profiles[b].curve.branch(side, body.surfaceFractions[side][i - body.leadingIndex], profiles[b].stag).point;
      assert.equal(p.x, automatic.guides[b][side][i].x); assert.equal(p.y, automatic.guides[b][side][i].y);
    }
    for (let i = 0; i < automatic.diagnostics.x.length; i++) if (i <= body.leadingIndex || i >= body.trailingIndex) {
      const p = automatic.guides[b].upper[i], q = automatic.guides[b].lower[i];
      assert.equal(p.x, q.x); assert.equal(p.y, q.y);
      assert.ok(Math.abs(p.potential - q.potential - (i >= body.trailingIndex ? values[b].upper - values[b].lower : 0)) < 1e-13);
    }
  }
  assert.deepEqual(joined.diagnostics.x, localFine.diagnostics.x);
  assert.deepEqual(joined.outer, localFine.outer);
  assert.equal(joined.diagnostics.endpointSpacing.length, 2);
  for (let b = 0; b < profiles.length; b++) {
    const body = joined.input.bodies[b];
    for (const side of ['upper', 'lower']) {
      for (const hit of joined.diagnostics.panelBlocks.constraints[side][b]) {
        const i = joined.diagnostics.x.indexOf(hit.rank);
        assert.ok(Math.hypot(joined.guides[b][side][i].x - localFine.guides[b][side][i].x,
          joined.guides[b][side][i].y - localFine.guides[b][side][i].y) < 1e-14);
      }
      for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
        const expected = profiles[b].curve.branch(side, body.surfaceFractions[side][i - body.leadingIndex], profiles[b].stag).point;
        assert.equal(joined.guides[b][side][i].x, expected.x); assert.equal(joined.guides[b][side][i].y, expected.y);
      }
    }
    for (let i = 0; i < joined.diagnostics.x.length; i++) if (i <= body.leadingIndex || i >= body.trailingIndex) {
      const a = joined.guides[b].upper[i], c = joined.guides[b].lower[i];
      assert.equal(a.x, c.x); assert.equal(a.y, c.y);
      assert.ok(Math.abs(a.potential - c.potential - (i < body.trailingIndex ? 0 : values[b].upper - values[b].lower)) < 1e-13);
    }
  }
  assert.deepEqual(local.diagnostics.x, physical.diagnostics.x, 'placement comparison must retain all counts and ranks');
  assert.deepEqual(local.outer, physical.outer);
  for (let b = 0; b < profiles.length; b++) for (const side of ['upper', 'lower']) {
    const body = local.input.bodies[b];
    for (let i = 0; i < local.diagnostics.x.length; i++)
      if (i <= body.leadingIndex || i >= body.trailingIndex)
        assert.deepEqual(local.guides[b][side][i], physical.guides[b][side][i], 'shared cuts must stay unchanged');
    for (const hit of local.diagnostics.panelBlocks.constraints[side][b]) {
      const i = local.diagnostics.x.indexOf(hit.rank);
      assert.deepEqual(local.guides[b][side][i], physical.guides[b][side][i], 'actual physical hits must stay unchanged');
    }
    for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
      const p = profiles[b].curve.branch(side, body.surfaceFractions[side][i - body.leadingIndex], profiles[b].stag).point;
      assert.equal(local.guides[b][side][i].x, p.x); assert.equal(local.guides[b][side][i].y, p.y);
    }
  }
  assert.equal(local.diagnostics.surfaceDensity.length, 4);
  assert.ok(local.diagnostics.surfaceDensity.every(d => !d.facingReconciled && !d.joinDerivativeMatched));
  assert.notDeepEqual(local.input.bodies.map(b => b.surfaceFractions), physical.input.bodies.map(b => b.surfaceFractions));
  assert.ok(samples > 0);
  assert.deepEqual(physical.diagnostics.anchors, legacy.diagnostics.anchors);
  assert.deepEqual(physical.diagnostics.panelBlocks, legacy.diagnostics.panelBlocks);
  const index = (result, rank) => result.diagnostics.x.indexOf(rank);
  for (const [b, body] of physical.input.bodies.entries()) {
    for (const side of ['upper', 'lower']) {
      // Every wall sample stays on the original curve. Free station positions
      // can change; actual passage intersections cannot.
      for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
        const value = profiles[b].curve.branch(side, body.surfaceFractions[side][i - body.leadingIndex], profiles[b].stag);
        const point = physical.guides[b][side][i];
        assert.equal(point.x, value.point.x); assert.equal(point.y, value.point.y);
      }
      for (const hit of physical.diagnostics.panelBlocks.constraints[side][b])
        assert.deepEqual(physical.guides[b][side][index(physical, hit.rank)], legacy.guides[b][side][index(legacy, hit.rank)]);
    }
    for (let i = 0; i < physical.diagnostics.x.length; i++) {
      if (i > body.leadingIndex && i < body.trailingIndex) continue;
      const lower = physical.guides[b].lower[i], upper = physical.guides[b].upper[i];
      assert.equal(lower.x, upper.x); assert.equal(lower.y, upper.y);
      assert.ok(Math.abs(upper.potential - lower.potential - (i < body.trailingIndex ? 0 : values[b].upper - values[b].lower)) < 1e-13);
    }
  }
  for (let side = 0; side < 2; side++) for (const hit of physical.diagnostics.panelBlocks.constraints.outer[side])
    assert.deepEqual(physical.outer[side][index(physical, hit.rank)], legacy.outer[side][index(legacy, hit.rank)]);
  assert.equal(physical.diagnostics.cutSpacing.length, 4);
  for (const cut of physical.diagnostics.cutSpacing) {
    assert.ok(cut.achievedEdgeChord > 0 && Number.isFinite(cut.achievedEdgeChord));
    assert.equal(cut.exactDiscreteJoinSpacing, false);
    assert.ok(cut.x.every((xx, i) => !i || xx > cut.x[i - 1]));
  }
  // Independently evaluate the declared continuous maps near stagnation.
  // Distance is linear to first order in rank and potential quadratic, not
  // the former square-root-distance map with finite potential slope.
  const fit = physical.diagnostics.facingOutlineSpacing;
  assert.equal(fit.exactModernMsetLaw, false);
  const maps = fit.outlines.map((outline, m) => createMonotoneCubicMap(outline.knots, outline.values,
    { derivatives: 'prescribed', slopes: fit.slopes[m] }));
  for (const [b, profile] of profiles.entries()) {
    const le = physical.diagnostics.panelBlocks.rank[`${b}:LE`], h = 1e-7, map = maps[2 * b];
    const d1 = -map.value(le - h), d2 = -map.value(le - 2 * h);
    assert.ok(Math.abs(d2 / d1 - 2) < 1e-5);
    const edge = profile.curve.branch('upper', 0, profile.stag).point;
    // Use larger steps when subtracting potentials so the quadratic signal
    // is resolved above double-precision cancellation at the O(1) origin.
    const phiErrors = [1e-3, 1e-4].map(step => {
      const phi1 = profile.phiStag - laws.get(profile.upstream).phi(edge.x + map.value(le - step));
      const phi2 = profile.phiStag - laws.get(profile.upstream).phi(edge.x + map.value(le - 2 * step));
      return Math.abs(phi2 / phi1 - 4);
    });
    assert.ok(phiErrors[1] < phiErrors[0] && phiErrors[1] < .01);
    for (const end of ['upstream', 'wake']) for (const cut of physical.diagnostics.cutSpacing.filter(c => c.body === b && c.end === end))
      for (const block of cut.blocks) {
        const a = physical.diagnostics.x[block.fromIndex], c = physical.diagnostics.x[block.toIndex];
        const rank = a + .37 * (c - a), upper = maps[2 * b], lower = maps[2 * b + 1];
        const offset = upper.value(a) - lower.value(a);
        assert.ok(Math.abs(upper.value(rank) - lower.value(rank) - offset) < 2e-14);
      }
  }
});
