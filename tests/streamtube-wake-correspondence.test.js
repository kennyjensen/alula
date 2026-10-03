import test from 'node:test';
import assert from 'node:assert/strict';
import { initializeStreamtubeWakeCorrespondence } from '../src/euler/streamtube-wake-correspondence.js';
import { prepareConvexWakeGridUpdate, requireConvexGridUpdate } from '../src/euler/streamtube-grid-update.js';

const close = (a, b, tolerance = 8e-14) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const pointClose = (a, b, tolerance) => { close(a.x, b.x, tolerance); close(a.y, b.y, tolerance); };
const midpoint = (grid, b, i) => {
  const a = grid[b][i].at(-1), c = grid[b + 1][i][0];
  return { x: (a.x + c.x) / 2, y: (a.y + c.y) / 2 };
};
// Independent vector projections; no production gap helper in the oracle.
function measure(grid, b, i, forward = false) {
  const a = midpoint(grid, b, forward ? i : i - 1), c = midpoint(grid, b, Math.min(i + 1, grid[0].length - 1));
  const length = Math.hypot(c.x - a.x, c.y - a.y), tx = (c.x - a.x) / length, ty = (c.y - a.y) / length;
  const lower = grid[b][i].at(-1), upper = grid[b + 1][i][0], dx = upper.x - lower.x, dy = upper.y - lower.y;
  return { tangent: { x: tx, y: ty }, tau: dx * tx + dy * ty, gap: dy * tx - dx * ty };
}
function fixture({ elements = 2, curved = true, tau = .3, thickness = .4 } = {}) {
  const nx = 8, tubes = elements === 1 ? [3, 2] : [3, 4, 2];
  const bodies = Array.from({ length: elements }, (_, b) => ({ leadingIndex: 1 + b, trailingIndex: 3 + 2 * b }));
  const massFractions = elements === 1 ? [[0, .125, .625, 1], [0, .375, 1]]
    : [[0, .125, .625, 1], [0, .0625, .3125, .75, 1], [0, .375, 1]];
  const cut = (b, i, sign) => {
    const te = bodies[b].trailingIndex, dx = i === te ? tau * (1 + .25 * b) : 0;
    return { x: i + sign * dx / 2, y: 3 * b + (curved ? .025 * i * i : 0) + sign * thickness / 2 };
  };
  const nodes = tubes.map((nt, g) => Array.from({ length: nx + 1 }, (_, i) => {
    const lower = g === 0 ? { x: i, y: -3 } : cut(g - 1, i, 1);
    const upper = g === elements ? { x: i, y: 3 * elements } : cut(g, i, -1);
    return Array.from({ length: nt + 1 }, (_, j) => ({ x: lower.x + j / nt * (upper.x - lower.x), y: lower.y + j / nt * (upper.y - lower.y) }));
  }));
  return { nodes, layout: { nx, elements, tubes, bodies, independentWakeBanks: true }, massFractions };
}

test('a first-wake trial can correct bank correspondence without changing its mean or gap', () => {
  const { nodes: previous, layout, massFractions } = fixture({ elements: 1, curved: false, tau: 0 });
  const proposed = structuredClone(previous), i = layout.bodies[0].trailingIndex + 1;
  proposed[0][i].at(-1).x = i - 1.1;
  proposed[1][i][0].x = i + 1.1;
  const before = structuredClone({ previous, proposed, layout, massFractions });
  const allocation = { groups: massFractions.map(row => row.slice(1).map((v, j) => ({ massFlow: v - row[j] }))) };
  assert.throws(() => requireConvexGridUpdate(previous, proposed), { code: 'streamtube-grid-step' });
  // Force the recovery path on this deliberately simple folded geometry.
  const trial = prepareConvexWakeGridUpdate(previous, proposed, layout, allocation, 1);
  assert.ok(trial.correction.fraction > 0 && trial.correction.fraction <= 1);
  assert.equal(trial.correction.equationsChanged, false);
  assert.equal(trial.correction.physicalAcceptanceRequired, true);
  assert.equal(requireConvexGridUpdate(previous, trial.nodes).limited, false);
  for (let station = i; station <= layout.nx; station++) {
    pointClose(midpoint(trial.nodes, 0, station), midpoint(proposed, 0, station));
    close(measure(trial.nodes, 0, station).gap, measure(proposed, 0, station).gap);
  }
  for (let g = 0; g < previous.length; g++) for (let station = 0; station < i; station++)
    assert.deepEqual(trial.nodes[g][station], proposed[g][station]);
  assert.deepEqual({ previous, proposed, layout, massFractions }, before);
});

test('safe and non-wake grid trials retain their original convexity policy', () => {
  const { nodes, layout } = fixture({ elements: 1, curved: false, tau: 0 });
  assert.equal(prepareConvexWakeGridUpdate(nodes, nodes, layout).nodes, nodes);
  const proposed = structuredClone(nodes);
  proposed[0][1][1].x = 4;
  let original;
  try { requireConvexGridUpdate(nodes, proposed); } catch (error) { original = error; }
  assert.ok(original);
  assert.throws(() => prepareConvexWakeGridUpdate(nodes, proposed, layout), error =>
    error.code === original.code && error.stepFraction === original.stepFraction);
});

test('a first-wake contact with a usable raw step preserves ordinary backtracking', () => {
  const { nodes, layout, massFractions } = fixture({ elements: 1, curved: false, tau: 0 });
  const proposed = structuredClone(nodes), i = layout.bodies[0].trailingIndex + 1;
  proposed[0][i].at(-1).x = i - 1.1;
  proposed[1][i][0].x = i + 1.1;
  const before = structuredClone(proposed);
  const allocation = { groups: massFractions.map(row => row.slice(1).map((v, j) => ({ massFlow: v - row[j] }))) };
  let original;
  try { requireConvexGridUpdate(nodes, proposed); } catch (error) { original = error; }
  assert.ok(original.stepFraction > 2 ** -12);
  assert.throws(() => prepareConvexWakeGridUpdate(nodes, proposed, layout, allocation), error =>
    error.code === original.code && error.stepFraction === original.stepFraction
    && !error.diagnostics.coordinateRepairBacktrack);
  assert.deepEqual(proposed, before);
});

test('an invalid paired wake requests Newton backtracking without accepting or clipping its gap', () => {
  const { nodes, layout, massFractions } = fixture({ elements: 1, curved: false, tau: 0 });
  const proposed = structuredClone(nodes), i = layout.bodies[0].trailingIndex + 1;
  proposed[0][i].at(-1).x = i - 1.1;
  proposed[1][i][0].x = i + 1.1;
  proposed[1][i][0].y = proposed[0][i].at(-1).y - 1;
  const before = structuredClone(proposed);
  const allocation = { groups: massFractions.map(row => row.slice(1).map((v, j) => ({ massFlow: v - row[j] }))) };
  assert.throws(() => prepareConvexWakeGridUpdate(nodes, proposed, layout, allocation, 1), error =>
    error.code === 'streamtube-grid-step' && error.stepFraction === .5
    && error.diagnostics.pairingFailure?.code === 'WAKE_CORRESPONDENCE_GEOMETRY'
    && Number.isFinite(error.diagnostics.rawGridStepFraction));
  assert.deepEqual(proposed, before);
});

test('curved multielement pairing preserves gaps and gives both banks the center-edge advance', () => {
  const input = fixture(), before = structuredClone(input), result = initializeStreamtubeWakeCorrespondence(input);
  assert.deepEqual(input, before);
  const { nodes, layout } = input, moved = result.nodes;
  for (let b = 0; b < layout.elements; b++) {
    const te = layout.bodies[b].trailingIndex;
    for (const [g, j] of [[b, layout.tubes[b]], [b + 1, 0]]) for (let i = 0; i <= te; i++)
      assert.deepEqual(moved[g][i][j], nodes[g][i][j]);
    for (let i = te + 1; i <= layout.nx; i++) {
      const old = measure(nodes, b, i), next = measure(moved, b, i);
      pointClose(midpoint(moved, b, i), midpoint(nodes, b, i));
      close(next.gap, old.gap);
      const left = midpoint(nodes, b, i - 1), right = midpoint(nodes, b, i);
      const dx = right.x - left.x, dy = right.y - left.y, length = Math.hypot(dx, dy);
      for (const [g, j] of [[b, layout.tubes[b]], [b + 1, 0]]) {
        const p = moved[g][i][j], q = moved[g][i - 1][j];
        close(((p.x - q.x) * dx + (p.y - q.y) * dy) / length, length);
        const shift = { x: p.x - nodes[g][i][j].x, y: p.y - nodes[g][i][j].y };
        close(shift.y * old.tangent.x - shift.x * old.tangent.y, 0);
      }
    }
  }
  for (let i = 0; i <= layout.nx; i++) {
    assert.deepEqual(moved[0][i][0], nodes[0][i][0]);
    assert.deepEqual(moved.at(-1)[i].at(-1), nodes.at(-1)[i].at(-1));
  }
  for (let g = 0; g < nodes.length; g++) assert.deepEqual(moved[g][0], nodes[g][0]);
  assert.ok(result.diagnostics.maximumNodeDisplacement > 0);
  assert.equal(result.diagnostics.initialGuessOnly, true);
  assert.equal(result.diagnostics.equationsChanged, false);
  assert.equal(result.diagnostics.geometryAcceptanceRequired, true);
  close(result.diagnostics.maximumGapChange, 0);
  close(result.diagnostics.maximumCenterChange, 0);
  close(result.diagnostics.maximumTangentialMismatch, 0);
  close(result.diagnostics.maximumForwardAdvanceMismatch, 0);
  assert.ok(result.diagnostics.minimumFrameProjection > 0);
  assert.ok(result.diagnostics.minimumBankForwardAdvance > 0);
  assert.equal(result.diagnostics.coordinateRule, 'equal-center-edge-projected-bank-advance');
  moved[0][0][0].x = 99; result.diagnostics.bodies[0].firstAfter.tangent.x = 99;
  assert.deepEqual(input, before);
});

test('rotating first-wake frame retains forward clearance with unequal displaced TE endpoints', () => {
  const input = fixture({ elements: 1, curved: false, tau: 0, thickness: 0 });
  const te = input.layout.bodies[0].trailingIndex;
  const centers = Array.from({ length: input.layout.nx + 1 }, (_, i) => i <= te
    ? { x: (i - te) * .1, y: 0 } : i === te + 1 ? { x: .00004, y: 0 }
      : { x: .001 * (i - te - 1), y: -.0001 * (i - te - 1) });
  for (let i = 0; i <= input.layout.nx; i++) {
    const c = centers[i];
    let d = { x: .006, y: .012 };
    if (i > te) {
      const a = centers[i - 1], z = centers[Math.min(i + 1, input.layout.nx)];
      const h = Math.hypot(z.x - a.x, z.y - a.y), tx = (z.x - a.x) / h, ty = (z.y - a.y) / h;
      d = { x: .0048 * tx - .0125 * ty, y: .0048 * ty + .0125 * tx };
    }
    input.nodes[0][i].at(-1).x = c.x - d.x / 2;
    input.nodes[0][i].at(-1).y = c.y - d.y / 2;
    input.nodes[1][i][0] = { x: c.x + d.x / 2, y: c.y + d.y / 2 };
  }
  const first = te + 1, original = measure(input.nodes, 0, first), tauTE = .006;
  // Copying the TE component into this rotating tangent would move the lower
  // first point behind its TE. This checks the replaced rule independently.
  const oldRuleLowerX = input.nodes[0][first].at(-1).x
    - .5 * (tauTE - original.tau) * original.tangent.x;
  assert.ok(oldRuleLowerX < input.nodes[0][te].at(-1).x);
  const { nodes } = initializeStreamtubeWakeCorrespondence(input);
  for (const [g, j] of [[0, 3], [1, 0]]) close(nodes[g][first][j].x - nodes[g][te][j].x, .00004);
  close(measure(nodes, 0, first).gap, .0125);
  // The appropriate local component differs from the TE component.
  assert.ok(Math.abs(measure(nodes, 0, first).tau - tauTE) > .0005);
});

test('curved pairing is idempotent to coordinate roundoff with staggered shared passages', () => {
  const input = fixture(), first = initializeStreamtubeWakeCorrespondence(input);
  const second = initializeStreamtubeWakeCorrespondence({ ...input, nodes: first.nodes });
  first.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) => pointClose(p, second.nodes[k][i][j]))));
  assert.ok(second.diagnostics.maximumNodeDisplacement < 8e-14);
});

test('interior displacement uses unequal cumulative mass fractions from both neighboring wakes', () => {
  const input = fixture(), { nodes: moved } = initializeStreamtubeWakeCorrespondence(input);
  const { nodes, layout, massFractions } = input;
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i <= layout.nx; i++) for (let j = 1; j < layout.tubes[g]; j++) {
    const eta = massFractions[g][j];
    for (const key of ['x', 'y']) close(moved[g][i][j][key] - nodes[g][i][j][key],
      (1 - eta) * (moved[g][i][0][key] - nodes[g][i][0][key])
      + eta * (moved[g][i].at(-1)[key] - nodes[g][i].at(-1)[key]));
  }
  assert.notEqual(moved[1][6][1].x - nodes[1][6][1].x,
    .75 * (moved[1][6][0].x - nodes[1][6][0].x) + .25 * (moved[1][6].at(-1).x - nodes[1][6].at(-1).x));
});

test('coordinate initialization commutes with rotation, translation and reference-length scaling', () => {
  const input = fixture(), expected = initializeStreamtubeWakeCorrespondence(input);
  for (const angle of [-1.2, .4, 2.7]) for (const scale of [.001, 1, 130]) {
    const move = p => ({ x: scale * (Math.cos(angle) * p.x - Math.sin(angle) * p.y) + 1.5,
      y: scale * (Math.sin(angle) * p.x + Math.cos(angle) * p.y) - .7 });
    const transformed = { ...input, nodes: input.nodes.map(g => g.map(row => row.map(move))) };
    const actual = initializeStreamtubeWakeCorrespondence(transformed);
    actual.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) =>
      pointClose(p, move(expected.nodes[k][i][j]), 3e-12 * Math.max(1, scale)))));
    actual.diagnostics.bodies.forEach((body, b) => close(body.tangentialOffset,
      scale * expected.diagnostics.bodies[b].tangentialOffset, 3e-12 * Math.max(1, scale)));
  }
});

test('zero tangential offset and zero-displacement gaps are exact detached no-ops', () => {
  for (const thickness of [0, .5]) {
    const input = fixture({ elements: 1, curved: false, tau: 0, thickness });
    input.nodes[0][0][0].x = -0;
    const { nodes, diagnostics } = initializeStreamtubeWakeCorrespondence(input);
    assert.deepEqual(nodes, input.nodes); assert.notEqual(nodes, input.nodes);
    assert.equal(diagnostics.maximumNodeDisplacement, 0);
    assert.equal(diagnostics.maximumGapChange, 0); assert.equal(diagnostics.maximumCenterChange, 0);
    assert.equal(diagnostics.maximumTangentialMismatch, 0);
    assert.equal(diagnostics.minimumOutputGap, thickness);
  }
});

test('aligned straight nonzero tau remains an exact no-op when binary representable', () => {
  const input = fixture({ elements: 1, curved: false, tau: .5, thickness: .5 });
  for (let i = input.layout.bodies[0].trailingIndex + 1; i <= input.layout.nx; i++) {
    input.nodes[0][i].at(-1).x -= .25; input.nodes[1][i][0].x += .25;
  }
  assert.deepEqual(initializeStreamtubeWakeCorrespondence(input).nodes, input.nodes);
});

test('coincident zero-gap bank aliases receive separate opposite shifts', () => {
  const input = fixture({ elements: 1, curved: false, tau: .5, thickness: .5 });
  for (let i = input.layout.bodies[0].trailingIndex + 1; i <= input.layout.nx; i++) {
    const center = { x: i, y: 0 };
    input.nodes[0][i][3] = center; input.nodes[1][i][0] = center;
  }
  const before = structuredClone(input), { nodes } = initializeStreamtubeWakeCorrespondence(input);
  for (let i = input.layout.bodies[0].trailingIndex + 1; i <= input.layout.nx; i++) {
    assert.notEqual(nodes[0][i][3], nodes[1][i][0]);
    assert.deepEqual(nodes[0][i][3], { x: i - .25, y: 0 });
    assert.deepEqual(nodes[1][i][0], { x: i + .25, y: 0 });
    assert.equal(measure(nodes, 0, i).gap, 0);
  }
  assert.deepEqual(input, before);
});

test('degenerate tangents, reversed gaps, invalid dimensions and malformed mass fractions reject without mutation', () => {
  const reject = (change, pattern) => {
    const input = fixture(); change(input); const before = structuredClone(input);
    assert.throws(() => initializeStreamtubeWakeCorrespondence(input), pattern); assert.deepEqual(input, before);
  };
  reject(x => { x.nodes[0][0][0].x = NaN; }, /coordinates/);
  reject(x => { x.nodes[0].pop(); }, /coordinates/);
  reject(x => { x.layout.independentWakeBanks = false; }, /layout/);
  reject(x => { x.layout.bodies[0].trailingIndex = x.layout.nx; }, /layout/);
  reject(x => { x.massFractions[0][1] = 0; }, /cumulative/);
  reject(x => { x.massFractions[0][1] = .8; }, /cumulative/);
  reject(x => { x.massFractions[0].pop(); }, /cumulative/);
  reject(x => { x.massFractions[1][0] = .01; }, /cumulative/);
  reject(x => {
    const te = x.layout.bodies[0].trailingIndex;
    x.nodes[0][te + 1][3] = { ...x.nodes[0][te][3] };
    x.nodes[1][te + 1][0] = { ...x.nodes[1][te][0] };
  }, /Degenerate first wake-center segment/);
  reject(x => {
    const te = x.layout.bodies[0].trailingIndex;
    for (const g of [0, 1]) x.nodes[g][te + 2][g ? 0 : 3] = { ...x.nodes[g][te][g ? 0 : 3] };
  }, /Degenerate wake centerline secant/);
  reject(x => { x.nodes[1][4][0].y = x.nodes[0][4][3].y - 1; }, /Negative normal wake gap/);
  reject(x => { x.nodes[1][3][0].y = x.nodes[0][3][3].y - 1; }, /Negative normal TE gap/);
});

test('singular and backward centered frames reject explicitly without clipping or mutation', () => {
  for (const endpoint of [{ x: 3, y: 1 }, { x: 2, y: 0 }]) {
    const input = fixture({ elements: 1, curved: false, tau: 0, thickness: 0 });
    // TE=(3,0), first=(4,0); the following center makes the centered tangent
    // perpendicular to, or backwards along, the positive first center edge.
    input.nodes[0][5][3] = { ...endpoint }; input.nodes[1][5][0] = { ...endpoint };
    const before = structuredClone(input);
    assert.throws(() => initializeStreamtubeWakeCorrespondence(input), error =>
      error.code === 'WAKE_CORRESPONDENCE_GEOMETRY' && error.body === 0 && error.station === 4
      && /Singular or backward/.test(error.message));
    assert.deepEqual(input, before);
  }
});
