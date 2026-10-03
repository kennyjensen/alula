import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { streamtubeEdgePressure } from '../src/euler/streamtube-edge-velocity.js';

const maximum = a => Math.max(...Array.from(a, Math.abs));
for (const edgeMatching of ['pressure', 'section-velocity']) {
const create = () => createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), { edgeMatching });
test(`${edgeMatching}: complete two-element Euler/BL Jacobian matches all independent state columns across geometry rebasing`, t => {
  const system = create(), { n, ne } = system;
  let x = system.initial.map((v, i) => i < ne ? v + 1e-7 * Math.sin(i) : v * (1 + .0005 * Math.sin(i)));
  assert.equal(system.bl.surfaces.length, 4); assert.equal(system.bl.wakes.length, 2);
  // Edge matching remains pointwise in BL velocity/pressure. The section-velocity
  // mode averages Euler section speeds, never adjacent BL unknowns.
  const base = system.residual(x), alternating = x.slice(), epsilon = 1e-5;
  [...system.bl.surfaces, ...system.bl.wakes].forEach(branch => branch.ids.forEach((id, k) => {
    alternating[ne + 4 * id + 3] += (k % 2 ? -1 : 1) * epsilon;
  }));
  const changed = system.residual(alternating);
  system.bl.stations.forEach(({ id }) => {
    const row = ne + 4 * id + 3;
    const expected = edgeMatching === 'pressure'
      ? streamtubeEdgePressure(alternating[row], system.euler.conditions).pressure - streamtubeEdgePressure(x[row], system.euler.conditions).pressure
      : alternating[row] - x[row];
    assert.ok(Math.abs(changed[row] - base[row] - expected) < 1e-13);
  });
  for (let chart = 0; chart < 2; chart++) {
    const jacobian = system.jacobian(x, { sparse: false });
    const errors = Array.from({ length: 4 }, () => ({ error: 0 }));
    const strength = [0, 0];
    for (let col = 0; col < n; col++) {
      const h = 2e-7 * Math.max(1, Math.abs(x[col])), plus = x.slice(), minus = x.slice(); plus[col] += h; minus[col] -= h;
      const p = system.residual(plus), m = system.residual(minus);
      for (let row = 0; row < n; row++) {
        const exact = jacobian[row * n + col], fd = (p[row] - m[row]) / (2 * h);
        const error = Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd));
        const block = 2 * Number(row >= ne) + Number(col >= ne);
        if (error > errors[block].error) errors[block] = { error, row, col, exact, fd };
        if (row < ne && col >= ne) strength[0] += Math.abs(exact);
        if (row >= ne && (row - ne) % 4 === 3 && col < ne) strength[1] += Math.abs(exact);
      }
    }
    for (const [block, worst] of errors.entries()) assert.ok(worst.error < 5e-6, JSON.stringify({ chart, block, ...worst }));
    // Thickness columns carry the 1/sqrt(Re) scale; their sum need not be O(1).
    assert.ok(strength.every(s => s > 1e-8), 'both displacement-to-Euler and Euler-to-edge coupling must be active');
    t.diagnostic(JSON.stringify({ chart, errors, couplingStrength: strength }));
    const before = system.evaluate(x), rebased = system.rebase(x), after = system.evaluate(rebased);
    assert.ok(maximum(before.residual.map((r, i) => r - after.residual[i])) < 2e-11);
    before.outer.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
      const q = after.outer.nodes[g][i][j]; assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 5e-13);
    })));
    x = rebased;
  }
});

test(`${edgeMatching}: one simultaneous two-element solve closes every BL and wake with independent per-region conservation`, t => {
  const system = create(), r = solveCoupledStreamtubeBody(system, { maxIterations: 8, tolerance: 1e-10 });
  assert.equal(r.converged, true, r.reason); assert.ok(maximum(r.residual) <= 1e-10);
  assert.equal(r.mesh.initialization.flowSolved, true);
  assert.ok(r.history.length > 1); assert.ok(r.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert.equal(r.boundaryLayer.stations.length, 31); assert.equal(r.boundaryLayer.surfaces.length, 4); assert.equal(r.boundaryLayer.wakes.length, 2);
  assert.equal(r.cl, null); assert.equal(r.cd, null); assert.equal(r.cm, null);
  t.diagnostic(JSON.stringify({ unknowns: system.n, iterations: r.history.length - 1, families: r.families, linear: r.linearDiagnostics }));
  const { flow } = r;
  for (let g = 0; g < flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: flow.nodes[g], sections: flow.sections.map(row => row[g]), cells: flow.cells.map(row => row[g]) }, system.euler.conditions.gamma);
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(maximum(c[key]) < 2e-9, `region ${g}, ${key}: ${c[key]}`);
  }
  for (const w of r.boundaryLayer.wakes) {
    const sides = r.boundaryLayer.surfaces.filter(s => s.body === w.body).map(s => r.boundaryLayer.stations[s.ids.at(-1)]);
    const first = r.boundaryLayer.stations[w.ids[0]], [a, b] = sides;
    assert.ok(Math.abs(first.theta - a.theta - b.theta) < 1e-12);
    assert.ok(Math.abs(first.deltaStar - a.deltaStar - b.deltaStar) < 1e-12);
    assert.ok(Math.abs(first.aux * first.theta - a.aux * a.theta - b.aux * b.theta) < 1e-12);
    for (const id of w.ids.slice(1)) {
      const s = r.boundaryLayer.stations[id], lower = flow.nodes[w.body][s.i].at(-1), upper = flow.nodes[w.body + 1][s.i][0];
      assert.ok(Math.abs(Math.hypot(upper.x - lower.x, upper.y - lower.y) - system.euler.conditions.lengthScale * s.deltaStar) < 2e-12);
    }
  }
  const restart = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    { edgeMatching, initialEuler: r.flow, initialBL: r.x.slice(system.ne) });
  const again = restart.evaluate(restart.initial);
  assert.ok(maximum(again.residual) < 1e-10);
  flow.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    const q = again.outer.nodes[g][i][j]; assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 5e-13);
  })));
});
}
