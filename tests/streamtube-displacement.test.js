import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { initializeStreamtubeBodyFromGrid } from '../src/euler/tests/streamtube-body-restart.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { streamtubeMeshConnectivity } from '../src/geometry/streamtube-mesh-connectivity.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { extendStreamtubeDisplacement } from '../src/euler/streamtube-displacement.js';

const close = (a, b, tol = 2e-11) => assert.ok(Math.abs(a - b) <= tol, `${a} != ${b}`);
const max = a => Math.max(...Array.from(a, Math.abs));
const pointClose = (a, b, tol) => { close(a.x, b.x, tol); close(a.y, b.y, tol); };

test('a displaced stagnation point does not overtake a clustered inlet cut', () => {
  const input = intrinsicBodyFixture({ elements: 1, bodySegments: 8, tubes: 3 });
  const body = input.bodies[0], le = body.leadingIndex;
  for (const row of [input.cutPaths[0], input.outerLower, input.outerUpper]) row[le - 1].x = -.00015;
  input.displacement = { surfaces: [{ upper: Array(body.trailingIndex - le + 1).fill(.0003),
    lower: Array(body.trailingIndex - le + 1).fill(.0003) }],
    wakes: [Array(input.outerLower.length - 1 - body.trailingIndex).fill(.0006)] };
  const system = createStreamtubeBodySystem(input), before = system.decode(system.initial);
  const original = extendStreamtubeDisplacement(system, system.initial);
  assert.equal(streamtubeMeshSnapshot({ system, nodes: original }).quality.valid, false);
  const repaired = extendStreamtubeDisplacement(system, system.initial, { inletDisplacement: true });
  assert.equal(streamtubeMeshSnapshot({ system, nodes: repaired }).quality.valid, true);
  assert.deepEqual(system.decode(system.initial), before);
  for (let g = 0; g < repaired.length; g++) {
    assert.deepEqual(repaired[g][0], original[g][0]);
    assert.deepEqual(repaired[g].slice(le), original[g].slice(le));
  }
  for (let i = 0; i <= le; i++) assert.deepEqual(repaired[0][i].at(-1), repaired[1][i][0]);
});
function fixture({ amplitude = .0003, ...options } = {}) {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, ...options });
  input.displacement = { surfaces: input.bodies.map((b, k) => {
    const count = b.trailingIndex - b.leadingIndex + 1;
    return Object.fromEntries(['upper', 'lower'].map((side, s) => [side, Array.from({ length: count }, (_, i) =>
      amplitude * (1 + .2 * k + .1 * (s + 1) * i / (count - 1)))]));
  }), wakes: input.bodies.map((b, k) => Array.from({ length: input.outerLower.length - 1 - b.trailingIndex }, (_, i) =>
    amplitude * (2 + .4 * k) * (1 + .03 * i))) };
  return input;
}

test('all element displacements point outward and wake banks have their prescribed total normal gap', () => {
  const input = fixture(), before = structuredClone(input), system = createStreamtubeBodySystem(input), { layout } = system;
  const initialNodes = system.decode(system.initial).nodes, chart = new Map(system.geometryChart().map(p => [p.column, p.normal]));
  for (const node of layout.positions.filter(p => p.kind === 'cut')) {
    const left = Math.max(0, node.i - 1), right = Math.min(layout.nx, node.i + 1), b = node.body;
    const a = initialNodes[b][left].at(-1), c = initialNodes[b][right].at(-1), u = initialNodes[b + 1][left][0], v = initialNodes[b + 1][right][0];
    const tx = c.x - a.x + v.x - u.x, ty = c.y - a.y + v.y - u.y, n = chart.get(node.column);
    close((tx * n.x + ty * n.y) / Math.hypot(tx, ty), 0, 1e-12);
  }
  const state = system.initial.map((v, i) => i < layout.densityCount ? v : 1e-5 * Math.sin(i));
  const { nodes, undisplacedNodes, stagnation } = system.decode(state), mesh = streamtubeMeshSnapshot({ system, nodes });
  assert.deepEqual(input, before); assert.equal(mesh.quality.valid, true);
  const { indices, fixed } = streamtubeMeshConnectivity(layout, nodes);
  for (let b = 0; b < layout.elements; b++) {
    const body = layout.bodies[b], te = body.trailingIndex;
    for (const side of ['upper', 'lower']) {
      const g = side === 'lower' ? b : b + 1, j = side === 'lower' ? layout.tubes[b] : 0;
      for (let i = body.leadingIndex; i <= te; i++) {
        const k = i - body.leadingIndex, raw = undisplacedNodes[g][i][j], p = nodes[g][i][j];
        const v = system.curves[b].branch(side, system.fractions[b][side][k], stagnation[b]).derivative;
        const dx = p.x - raw.x, dy = p.y - raw.y, length = Math.hypot(v.x, v.y);
        close(Math.hypot(dx, dy), input.displacement.surfaces[b][side][k]);
        close((dx * v.x + dy * v.y) / length, 0);
        assert.ok((dx * v.y - dy * v.x) / length > 0, 'A displacement must move into fluid.');
      }
    }
    assert.equal(indices[b][body.leadingIndex].at(-1), indices[b + 1][body.leadingIndex][0]);
    assert.notEqual(indices[b][te].at(-1), indices[b + 1][te][0], 'Separate TE displacement endpoints.');
    for (let i = te + 1; i <= layout.nx; i++) {
      const lo = nodes[b][i].at(-1), hi = nodes[b + 1][i][0], raw = undisplacedNodes[b][i].at(-1);
      close(Math.hypot(hi.x - lo.x, hi.y - lo.y), input.displacement.wakes[b][i - te - 1]);
      pointClose({ x: .5 * (lo.x + hi.x), y: .5 * (lo.y + hi.y) }, raw);
      assert.notEqual(indices[b][i].at(-1), indices[b + 1][i][0]);
      assert.ok(fixed.has(indices[b][i].at(-1)) && fixed.has(indices[b + 1][i][0]), 'Geometry-only repair must preserve the wake gap.');
    }
  }
  assert.throws(() => createStreamtubeBodySystem({ ...input, displacement: { ...input.displacement, wakes: [] } }), /every body/);
  const bad = structuredClone(input); bad.displacement.surfaces[0].upper[0] *= 2;
  assert.throws(() => createStreamtubeBodySystem(bad), /leading-edge/);
  bad.displacement.surfaces[0].upper[0] = -1;
  assert.throws(() => createStreamtubeBodySystem(bad), /negative/);
});

test('zero displacement recovers the inviscid residual and complete Jacobian', () => {
  const input = fixture({ amplitude: 0 }), a = createStreamtubeBodySystem(input);
  const { displacement, ...inviscid } = input, b = createStreamtubeBodySystem(inviscid);
  const state = a.initial.map((_, i) => 1e-6 * Math.sin(i));
  assert.deepEqual(a.decode(state).nodes, b.decode(state).nodes);
  assert.deepEqual(a.residual(state), b.residual(state));
  const ja = a.jacobian(state), jb = b.jacobian(state);
  ja.forEach((v, i) => close(v, jb[i], 2e-13));
});

test('the displaced two-body full Jacobian includes curved wall normals, both wake banks and centerline end tangency', () => {
  const system = createStreamtubeBodySystem(fixture()), { layout } = system;
  let state = system.initial.map((_, i) => 1e-5 * Math.sin(i));
  for (let chart = 0; chart < 2; chart++) {
    const jacobian = system.jacobian(state), geometry = system.geometryDerivatives(state);
    for (let col = 0; col < layout.n; col++) {
      const h = 2e-7, plus = state.slice(), minus = state.slice(); plus[col] += h; minus[col] -= h;
      const a = system.evaluate(plus), b = system.evaluate(minus);
      for (let row = 0; row < layout.n; row++) {
        const fd = (a.residual[row] - b.residual[row]) / (2 * h), exact = jacobian[row * layout.n + col];
        assert.ok(Math.abs(fd - exact) < 3e-6 * Math.max(1, Math.abs(exact)), `row ${row}, col ${col}: ${fd} != ${exact}`);
      }
      for (let g = 0; g <= layout.elements; g++) for (let i = 0; i <= layout.nx; i++) for (let j = 0; j <= layout.tubes[g]; j++) {
        const d = geometry[g][i][j].get(col) ?? { x: 0, y: 0 };
        for (const key of ['x', 'y']) close((a.nodes[g][i][j][key] - b.nodes[g][i][j][key]) / (2 * h), d[key], 1e-7);
      }
    }
    state = system.rebase(state);
  }
});

test('displacement geometry is invariant under rotation and scaling; physical restarts preserve the grid and reject changed wake gaps', () => {
  const input = fixture(), system = createStreamtubeBodySystem(input);
  const state = system.initial.map((_, i) => 1e-5 * Math.sin(i)), original = system.evaluate(state);
  const rebased = system.rebase(state), again = system.evaluate(rebased);
  original.residual.forEach((r, i) => close(r, again.residual[i], 2e-12));
  original.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) => pointClose(p, again.nodes[k][i][j]))));
  const restart = initializeStreamtubeBodyFromGrid(input, { ...again, strengths: again.strengths });
  const restored = restart.system.evaluate(restart.initial);
  restored.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) => pointClose(p, again.nodes[k][i][j]))));
  const chart = system.geometryChart(), bad = structuredClone(again.nodes), i = input.bodies[0].trailingIndex + 1;
  bad[1][i][0].y += .0001;
  assert.throws(() => system.adoptGeometry(rebased, bad), /wake banks/);
  assert.deepEqual(system.geometryChart(), chart);
  const angle = 23 * Math.PI / 180, factor = 3.7;
  const move = p => ({ x: 2 + factor * (Math.cos(angle) * p.x - Math.sin(angle) * p.y), y: -1 + factor * (Math.sin(angle) * p.x + Math.cos(angle) * p.y) });
  const moved = { ...input, alpha: input.alpha + 23, bodies: input.bodies.map(b => ({ ...b, points: b.points.map(move), stagnationParameter: b.stagnationParameter * factor })),
    outerLower: input.outerLower.map(move), outerUpper: input.outerUpper.map(move), cutPaths: input.cutPaths.map(row => row.map(move)),
    displacement: { surfaces: input.displacement.surfaces.map(s => Object.fromEntries(Object.entries(s).map(([side, row]) => [side, row.map(d => factor * d)]))),
      wakes: input.displacement.wakes.map(row => row.map(d => factor * d)) } };
  const target = createStreamtubeBodySystem(moved), value = target.evaluate(state);
  original.nodes.forEach((g, k) => g.forEach((row, i) => row.forEach((p, j) => pointClose(move(p), value.nodes[k][i][j], 1e-10))));
  // The least-squares farfield rows carry sqrt(length) from normalization.
  original.residual.forEach((r, i) => close(r * (system.layout.rows[i].kind === 'farfieldMatch' ? Math.sqrt(factor) : 1), value.residual[i], 2e-10));
});

test('analytic Euler-to-displacement columns include every element side and wake, including TE-driven wake-normal changes', () => {
  for (const flowModel of ['compressible', 'incompressible']) {
  const input = { ...fixture(), flowModel, mach: flowModel === 'incompressible' ? 0 : .2 }, system = createStreamtubeBodySystem(input), state = system.initial.map((_, i) => 1e-5 * Math.sin(i));
  const jacobian = system.jacobian(state, { includeDisplacement: true }), original = system.jacobian(state);
  assert.deepEqual(jacobian.state, original);
  assert.ok(jacobian.parameters.some(p => p.body === 0 && p.kind === 'wake'));
  assert.ok(jacobian.parameters.some(p => p.body === 1 && p.kind === 'wake'));
  for (const [col, p] of jacobian.parameters.entries()) {
    const h = 1e-7, plus = structuredClone(input), minus = structuredClone(input);
    const change = (target, delta) => {
      if (p.kind === 'wake') target.displacement.wakes[p.body][p.index] += delta;
      else for (const side of p.side === 'both' ? ['upper', 'lower'] : [p.side]) target.displacement.surfaces[p.body][side][p.index] += delta;
    };
    change(plus, h); change(minus, -h);
    // A displacement derivative belongs to one frozen Newton chart. A
    // freshly initialized geometry at each thickness would also rotate its
    // initial directions, differentiating a different coordinate problem.
    const plusSystem = createStreamtubeBodySystem(input), minusSystem = createStreamtubeBodySystem(input);
    plusSystem.setDisplacement(plus.displacement); minusSystem.setDisplacement(minus.displacement);
    const a = plusSystem.residual(state), b = minusSystem.residual(state);
    let effect = 0;
    for (let row = 0; row < system.layout.n; row++) {
      const exact = jacobian.displacement[row].get(col) ?? 0, fd = (a[row] - b[row]) / (2 * h);
      effect = Math.max(effect, Math.abs(exact));
      assert.ok(Math.abs(fd - exact) < 3e-6 * Math.max(1, Math.abs(exact)), `${JSON.stringify(p)}, row ${row}: ${fd} != ${exact}`);
    }
    assert.ok(effect > 1e-5, `Missing displacement feedback for ${JSON.stringify(p)}`);
  }
  }
});

test('two-element conservative Euler responds to every prescribed wall/wake displacement and conserves each fluid region', () => {
  const input = fixture({ bodySegments: 8, tubes: 3 }), system = createStreamtubeBodySystem(input);
  const result = solveStreamtubeBody(system, { tolerance: 1e-11, maxIterations: 20 });
  assert.equal(result.converged, true, result.reason); assert.ok(result.diagnostics.residual < 1e-11);
  assert.ok(result.linearDiagnostics.maxRelativeResidual <= 1e-10);
  assert.equal(streamtubeMeshSnapshot({ system, nodes: result.nodes }).quality.valid, true);
  assert.equal(result.surfaces.length, 4); assert.equal(result.displacement.wakes.length, 2);
  assert.match(result.forceStatus, /not solid-wall aerodynamic forces/); assert.match(result.formulation, /no BL coupling/);
  for (let g = 0; g < result.nodes.length; g++) {
    const c = directChannelConservation({ nodes: result.nodes[g], sections: result.sections.map(row => row[g]), cells: result.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'external', 'internalCancellation']) assert.ok(max(c[key]) < 2e-9, `${g}: ${key} ${c[key]}`);
  }
  // The same physical grid cannot solve the zero-displacement problem:
  // geometry feeds back into the Euler equations, not just into its plot.
  const { displacement, ...base } = input, plain = createStreamtubeBodySystem(base);
  assert.throws(() => plain.adoptGeometry(result.x, result.nodes), /wall does not match|disconnected/);
  const withoutDisplacement = plain.adoptGeometry(result.x, result.undisplacedNodes);
  assert.ok(max(plain.residual(withoutDisplacement)) > 1e-4);
});
