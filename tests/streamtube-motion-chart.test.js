import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeMotionDirections } from '../src/euler/streamtube-geometry.js';
import { createStreamtubeBodyLayout } from '../src/euler/streamtube-body-layout.js';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directBodyConservation } from './oracles/streamtube-body.js';

const near = (a, b, tolerance = 3e-13) => assert.ok(Math.abs(a - b) < tolerance, `${a} != ${b}`);
const eachPoint = (nodes, fn) => nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => fn(p, g, i, j))));

test('NCALC station stencil follows independent quadratic streamline secants and preserves rotation', () => {
  const layout = createStreamtubeBodyLayout({ segments: 6, tubes: [3, 3], bodies: [{ leadingIndex: 2, trailingIndex: 4 }] });
  const xs = [-2, -1.2, -.5, 0, .4, 1.3, 2.5];
  const grid = layout.nodes.map((group, g) => group.map((row, i) => row.map((_, j) => ({ x: xs[i], y: .1 * xs[i] ** 2 + g + j / 3 }))));
  for (const mode of ['centered', 'body-stations']) {
    const normals = streamtubeMotionDirections(layout, grid, mode);
    for (const node of layout.positions) {
      let a = xs[Math.max(0, node.i - 1)], b = xs[Math.min(layout.nx, node.i + 1)];
      if (mode === 'body-stations' && node.i === 2) b = xs[node.i];
      if (mode === 'body-stations' && node.i === 4) a = xs[node.i];
      // The secant slope of y=0.1*x^2 is exactly 0.1*(a+b).
      const slope = .1 * (a + b), n = normals.get(node.column);
      near(n.x, -slope / Math.hypot(1, slope)); near(n.y, 1 / Math.hypot(1, slope));
    }
    const angle = .7, c = Math.cos(angle), s = Math.sin(angle);
    const rotated = grid.map(group => group.map(row => row.map(p => ({ x: 3 * (c * p.x - s * p.y) + 4, y: 3 * (s * p.x + c * p.y) - 2 }))));
    const other = streamtubeMotionDirections(layout, rotated, mode);
    for (const [col, n] of normals) { near(other.get(col).x, c * n.x - s * n.y); near(other.get(col).y, s * n.x + c * n.y); }
  }
});

test('multielement station extension affects adjacent passages while unrelated cuts keep one centered direction', () => {
  const layout = createStreamtubeBodyLayout({ segments: 8, tubes: [2, 2, 2], bodies: [{ leadingIndex: 2, trailingIndex: 4 }, { leadingIndex: 5, trailingIndex: 7 }] });
  const grid = layout.nodes.map((group, g) => group.map((row, i) => row.map((_, j) => ({ x: i, y: .1 * i * i + g + j / 2 }))));
  const a = streamtubeMotionDirections(layout, grid), b = streamtubeMotionDirections(layout, grid, 'body-stations');
  assert.notDeepEqual(a.get(layout.nodes[1][2][1].column), b.get(layout.nodes[1][2][1].column));
  assert.deepEqual(a.get(layout.nodes[2][2][1].column), b.get(layout.nodes[2][2][1].column));
  const cut = layout.nodes[1][2][2].column;
  assert.equal(cut, layout.nodes[2][2][0].column); assert.deepEqual(a.get(cut), b.get(cut));
  const conflict = createStreamtubeBodyLayout({ segments: 8, tubes: [2, 2, 2], bodies: [{ leadingIndex: 2, trailingIndex: 4 }, { leadingIndex: 4, trailingIndex: 7 }] });
  assert.throws(() => streamtubeMotionDirections(conflict, grid, 'body-stations'), /Conflicting/);
});

test('walls-only stagnation motion leaves free physical nodes fixed and rebase/restart preserves their positions', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'walls-only', normalStencil: 'body-stations' };
  const system = createStreamtubeBodySystem(input), { layout } = system;
  const state = system.initial.slice(); state[layout.globals.stagnation[0]] = 1e-4; state[layout.globals.stagnation[1]] = -1e-4;
  const initial = system.decode(system.initial).nodes, moved = system.decode(state).nodes;
  let changedWall = 0;
  eachPoint(moved, (p, g, i, j) => {
    const q = initial[g][i][j], node = layout.nodes[g][i][j];
    if (node.column !== null) assert.deepEqual(p, q);
    else changedWall = Math.max(changedWall, Math.hypot(p.x - q.x, p.y - q.y));
  });
  assert.ok(changedWall > 1e-5);
  const rebased = system.rebase(state), after = system.decode(rebased).nodes;
  eachPoint(after, (p, g, i, j) => { near(p.x, moved[g][i][j].x); near(p.y, moved[g][i][j].y); });
  const restarted = createStreamtubeBodySystem(input), restored = restarted.adoptGeometry(rebased, after);
  eachPoint(restarted.decode(restored).nodes, (p, g, i, j) => { near(p.x, moved[g][i][j].x); near(p.y, moved[g][i][j].y); });
  const chartBefore = system.geometryChart(), invalid = structuredClone(after); invalid[0][layout.bodies[0].leadingIndex][2].x += .01;
  assert.throws(() => system.adoptGeometry(rebased, invalid), /wall/);
  assert.deepEqual(system.geometryChart(), chartBefore);
});

test('walls-only charts retain every residual Jacobian column for both equations, both stencils and rebasing', () => {
  for (const streamwiseMode of ['momentum', 'isentropic']) for (const normalStencil of ['centered', 'body-stations']) {
    const system = createStreamtubeBodySystem({ ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, alpha: .25 }), streamwiseMode, normalStencil, stagnationMotion: 'walls-only' });
    const { layout } = system, n = layout.n;
    let state = system.initial.map((_, i) => (i < layout.densityCount ? 1e-3 : 1e-6) * Math.sin(i + .3));
    for (let chart = 0; chart < 2; chart++) {
      const matrix = system.jacobian(state), h = 1e-7;
      let worst = 0;
      for (let col = 0; col < n; col++) {
        const plus = state.slice(), minus = state.slice(); plus[col] += h; minus[col] -= h;
        const a = system.residual(plus), b = system.residual(minus);
        for (let row = 0; row < n; row++) {
          const exact = matrix[row * n + col], fd = (a[row] - b[row]) / (2 * h);
          worst = Math.max(worst, Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd)));
        }
      }
      assert.ok(worst < 4e-7, JSON.stringify({ streamwiseMode, normalStencil, chart, worst }));
      state = system.rebase(state);
    }
  }
});

test('source-comparison motion with additive density closes a conservative two-element root', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'walls-only', normalStencil: 'body-stations' };
  const system = createStreamtubeBodySystem(input), r = solveStreamtubeBody(system, { stepMethod: 'density-newton', maxIterations: 16, tolerance: 1e-11 });
  assert.equal(r.converged, true, `${r.reason}, R=${r.diagnostics.residual}`);
  const balance = directBodyConservation(r, input.bodies, system.conditions);
  for (const v of [...balance.balance, ...balance.cutTraction]) assert.ok(Math.abs(v) < 2e-9);
});

test('independent motion retains displaced wall/wake geometry derivatives and restores separated banks', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), stagnationMotion: 'walls-only', normalStencil: 'body-stations' };
  input.displacement = {
    surfaces: input.bodies.map(b => Object.fromEntries(['lower', 'upper'].map(side => [side, Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, k) => 1e-4 * (1 + k * (side === 'upper' ? .3 : .2)))]))),
    wakes: input.bodies.map(b => Array.from({ length: input.outerLower.length - 1 - b.trailingIndex }, () => 5e-4)),
  };
  const system = createStreamtubeBodySystem(input), { layout } = system;
  let state = system.initial.map((_, i) => i < layout.densityCount ? 0 : 1e-6 * Math.sin(i + .3));
  for (let chart = 0; chart < 2; chart++) {
    const geometry = system.geometryDerivatives(state), h = 1e-7;
    for (let col = layout.densityCount; col < layout.n; col++) {
      const plus = state.slice(), minus = state.slice(); plus[col] += h; minus[col] -= h;
      const a = system.decode(plus).nodes, b = system.decode(minus).nodes;
      eachPoint(a, (p, g, i, j) => {
        const d = geometry[g][i][j].get(col) ?? { x: 0, y: 0 };
        for (const key of ['x', 'y']) near(d[key], (p[key] - b[g][i][j][key]) / (2 * h), 1e-7);
      });
    }
    const before = system.decode(state).nodes;
    state = system.rebase(state);
    eachPoint(system.decode(state).nodes, (p, g, i, j) => { near(p.x, before[g][i][j].x); near(p.y, before[g][i][j].y); });
  }
});
