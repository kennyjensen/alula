import test from 'node:test';
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { createStreamtubeBoundaryLayers } from '../src/euler/streamtube-boundary-layers.js';
import { seedStreamtubeWakeBanks, extendStreamtubeDisplacement } from '../src/euler/streamtube-displacement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { twoActiveFiniteBaseWakes } from './fixtures/two-active-finite-base-wakes.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { initializeStreamtubeBodyFromFlow } from '../src/euler/tests/streamtube-body-flow-restart.js';

const close = (a, b, tolerance = 2e-13) => assert.ok(Math.abs(a - b) <= tolerance, `${a} != ${b}`);
const max = a => a.reduce((v, x) => Math.max(v, Math.abs(x)), 0);
const center = (nodes, b, i, key) => .5 * (nodes[b][i].at(-1)[key] + nodes[b + 1][i][0][key]);
function fixture(amplitude = .0003) {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, contourPanels: 40 });
  return { ...input, wakeGeometry: 'independent-banks', wakeDisplacementMotion: 'te-center', displacement: {
    surfaces: input.bodies.map((b, k) => Object.fromEntries(['upper', 'lower'].map((side, s) => [side,
      Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => amplitude * (1 + .2 * k + .1 * (s + 1) * i))]))),
    wakes: input.bodies.map(b => Array(input.outerLower.length - 1 - b.trailingIndex).fill(amplitude * 3)) } };
}
function nodesClose(a, b, tolerance = 2e-14) {
  a.forEach((g, k) => g.forEach((r, i) => r.forEach((p, j) => {
    close(p.x, b[k][i][j].x, tolerance); close(p.y, b[k][i][j].y, tolerance);
  })));
}
function coupledFixture() {
  const { displacement, ...input } = fixture(0), euler = createStreamtubeBodySystem({ ...input, displacement });
  const eulerState = euler.initial.slice(), options = { transitionMode: 'automatic', tripFractions: [[1, 1], [1, 1]], edgeMatching: 'section-velocity' };
  const bl = createStreamtubeBoundaryLayers(euler, eulerState, options);
  const initialBL = new Float64Array(4 * bl.stations.length);
  for (const s of bl.stations) initialBL.set([['similarity', 'laminar'].includes(s.regime) ? .2 : .03,
    (s.kind === 'wake' ? .0004 : .0002) / bl.scale,
    (s.kind === 'wake' ? .001 : .0005) / bl.scale, 1 + .005 * Math.sin(s.id)], 4 * s.id);
  options.transitionState = bl.snapshotActive();
  const system = createCoupledStreamtubeBody(input, { ...options, initialBL,
    initialEuler: { x: eulerState, nodes: euler.decode(eulerState).nodes } });
  const x = system.initial.subarray(0, system.ne);
  system.initial.set(system.euler.adoptGeometry(x, seedStreamtubeWakeBanks(system.euler, x)));
  return { input, options, system };
}

test('the explicit TE chart translates both free wake banks once and preserves gap freedom', () => {
  const input = fixture(), system = createStreamtubeBodySystem(input), x = system.initial;
  assert.equal(system.layout.wakeDisplacementMotion, 'te-center');
  const before = system.decode(x), changed = structuredClone(input.displacement);
  changed.surfaces[0].upper[changed.surfaces[0].upper.length - 1] *= 2;
  system.setDisplacement(changed);
  const after = system.decode(x);
  for (let b = 0; b < system.layout.elements; b++) {
    const te = system.layout.bodies[b].trailingIndex;
    for (const key of ['x', 'y']) {
      const delta = center(after.nodes, b, te, key) - center(before.nodes, b, te, key);
      for (let i = te + 1; i <= system.layout.nx; i++) {
        close(after.nodes[b][i].at(-1)[key] - before.nodes[b][i].at(-1)[key], delta);
        close(after.nodes[b + 1][i][0][key] - before.nodes[b + 1][i][0][key], delta);
        close(after.nodes[b + 1][i][0][key] - after.nodes[b][i].at(-1)[key],
          before.nodes[b + 1][i][0][key] - before.nodes[b][i].at(-1)[key]);
      }
    }
  }
  changed.wakes[0][0] *= 1.5; system.setDisplacement(changed);
  assert.deepEqual(system.decode(x).nodes, after.nodes, 'Wake thickness controls its gap equation, not bank coordinates.');
  const node = system.layout.positions.find(p => p.kind === 'cut' && p.side === 'upper' && p.body === 0);
  const y = x.slice(); y[node.column] += .0001;
  assert.notDeepEqual(system.decode(y).nodes[node.body + 1][node.i][0], after.nodes[node.body + 1][node.i][0]);
});

test('cold seed, displacement extension and geometry adoption carry the mean TE shift exactly once', () => {
  for (const operation of [seedStreamtubeWakeBanks, extendStreamtubeDisplacement]) {
    const input = fixture(), system = createStreamtubeBodySystem(input), x = system.initial;
    const before = system.decode(x), moved = operation(system, x);
    for (let b = 0; b < system.layout.elements; b++) for (const key of ['x', 'y']) {
      const te = system.layout.bodies[b].trailingIndex;
      const delta = center(before.nodes, b, te, key) - center(before.undisplacedNodes, b, te, key);
      for (let i = te + 1; i <= system.layout.nx; i++) close(center(moved, b, i, key),
        center(before.undisplacedNodes, b, i, key) + delta, 2e-15);
    }
    const adopted = system.adoptGeometry(x, moved); nodesClose(system.decode(adopted).nodes, moved);
    const again = system.rebase(adopted); nodesClose(system.decode(again).nodes, moved);
    const restored = createStreamtubeBodySystem(input), y = restored.adoptGeometry(again, moved);
    nodesClose(restored.decode(y).nodes, moved);
  }
});

test('zero translation and omitted legacy mode retain their original residuals and wake bases', () => {
  const input = fixture(0), te = createStreamtubeBodySystem(input);
  const { wakeDisplacementMotion, ...legacy } = input, old = createStreamtubeBodySystem(legacy);
  assert.deepEqual(te.decode(te.initial).nodes, old.decode(old.initial).nodes);
  assert.deepEqual(te.residual(te.initial), old.residual(old.initial));
  assert.throws(() => createStreamtubeBodySystem({ ...input, wakeGeometry: 'centerline' }), /independent wake/);
  assert.throws(() => createStreamtubeBodySystem({ ...input, wakeDisplacementMotion: 'unknown' }), /TE-center/);
});

test('remote wake Euler rows retain both TE thickness derivatives through sparse assembly and rebasing', () => {
  const input = { ...fixture(), streamwiseMode: 'hybrid', hybrid: { ismom: 4, epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const system = createStreamtubeBodySystem(input), { layout } = system;
  let x = system.initial.map((_, i) => 1e-6 * Math.sin(i));
  for (let chart = 0; chart < 2; chart++) {
    const matrix = system.jacobian(x, { sparse: true, includeDisplacement: true });
    const dense = system.jacobian(x, { includeDisplacement: true });
    for (let r = 0; r < layout.n; r++) for (let p = matrix.state.rowPtr[r]; p < matrix.state.rowPtr[r + 1]; p++)
      close(matrix.state.values[p], dense.state[r * layout.n + matrix.state.colIndex[p]], 0);
    for (const [col, p] of system.displacementParameters.entries()) {
      const body = layout.bodies[p.body];
      if (p.kind !== 'surface' || p.index !== body.trailingIndex - body.leadingIndex) continue;
      const h = 2e-7, samples = [1, -1].map(sign => {
        const d = structuredClone(input.displacement); d.surfaces[p.body][p.side][p.index] += sign * h;
        system.setDisplacement(d); return system.residual(x);
      });
      let remote = 0;
      for (let r = 0; r < layout.n; r++) {
        const exact = matrix.displacement[r].get(col) ?? 0, fd = (samples[0][r] - samples[1][r]) / (2 * h);
        assert.ok(Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd)) < 3e-6, `TE ${p.body}/${p.side}, row ${r}`);
        if (layout.rows[r].i > body.trailingIndex + 2) remote = Math.max(remote, Math.abs(exact));
      }
      assert.ok(remote > 1e-6, 'The TE displacement column must reach remote wake rows.');
    }
    system.setDisplacement(input.displacement); x = system.rebase(x);
  }
});

test('full automatic coupled derivatives and checkpoint reconstruction include TE-driven wake arc motion', () => {
  const { input, options, system } = coupledFixture(), x = system.initial.slice();
  const initial = system.evaluate(x), matrix = system.jacobian(x);
  const d = x.map((v, k) => Math.sin(.53 * k + .2) * (k < system.ne ? .001 : .02 * Math.max(.1, Math.abs(v))));
  const exact = sparseProduct(matrix, d), h = 2e-6;
  const a = system.residual(x.map((v, k) => v + h * d[k])), b = system.residual(x.map((v, k) => v - h * d[k]));
  let error = 0;
  for (let k = 0; k < x.length; k++) {
    const fd = (a[k] - b[k]) / (2 * h);
    error = Math.max(error, Math.abs(exact[k] - fd) / Math.max(1, Math.abs(exact[k]), Math.abs(fd)));
  }
  assert.ok(error < 5e-6, `Full coupled directional error ${error}`);
  system.evaluate(x);
  const rebased = system.rebase(x), value = system.evaluate(rebased);
  assert.ok(max(initial.residual.map((v, k) => v - value.residual[k])) < 2e-11);
  const restart = JSON.parse(JSON.stringify({ input, options, initialEuler: { x: Array.from(rebased.subarray(0, system.ne)),
    nodes: value.outer.nodes, undisplacedNodes: value.outer.undisplacedNodes }, initialBL: Array.from(rebased.subarray(system.ne)) }));
  const replay = createCoupledStreamtubeBody(restart.input, { ...restart.options,
    initialEuler: restart.initialEuler, initialBL: restart.initialBL });
  assert.deepEqual(replay.evaluate(replay.initial).families, value.families);
  assert.deepEqual(replay.evaluate(replay.initial).outer.nodes, value.outer.nodes);
});

test('TE thickness sensitivities reach both finite-base BL wake rows and exact fluid-domain constraints', () => {
  const f = twoActiveFiniteBaseWakes(), before = f.system.evaluate(f.x);
  const system = createCoupledStreamtubeBody({ ...f.input, wakeDisplacementMotion: 'te-center' }, {
    initialEuler: { x: f.x.subarray(0, f.system.ne), undisplacedNodes: before.outer.undisplacedNodes },
    initialBL: f.x.subarray(f.system.ne), reynolds: 1e6, ncrit: 9, edgeMatching: 'section-velocity' });
  const x = system.initial.slice(), j = system.jacobian(x), constraints = system.stepConstraints(x);
  for (const surface of system.bl.surfaces) {
    const col = system.ne + 4 * surface.ids.at(-1) + 2, d = new Float64Array(system.n); d[col] = .1;
    const exact = sparseProduct(j, d), h = 2e-6;
    const plus = x.map((v, k) => v + h * d[k]), minus = x.map((v, k) => v - h * d[k]);
    const a = system.residual(plus), b = system.residual(minus);
    let wakeEffect = 0;
    for (let r = 0; r < system.n; r++) {
      const fd = (a[r] - b[r]) / (2 * h);
      assert.ok(Math.abs(exact[r] - fd) / Math.max(1, Math.abs(exact[r]), Math.abs(fd)) < 5e-6,
        `${surface.body}/${surface.side} row ${r}: ${exact[r]} != ${fd}`);
      if (r >= system.ne && system.bl.stations[Math.floor((r - system.ne) / 4)].kind === 'wake')
        wakeEffect = Math.max(wakeEffect, Math.abs(fd));
    }
    assert.ok(wakeEffect > 1e-6);
    const ca = system.constraintValues(plus), cb = system.constraintValues(minus);
    constraints.forEach((c, r) => {
      if (c.kind !== 'kinematic-shape') return;
      const exact = .1 * (c.gradient.get(col) ?? 0), fd = (ca[r] - cb[r]) / (2 * h);
      assert.ok(Math.abs(exact - fd) / Math.max(1, Math.abs(exact), Math.abs(fd)) < 5e-7,
        `${surface.body}/${surface.side} shape ${r}: ${exact} != ${fd}`);
    });
  }
});

test('normal and streamwise refinement retain physical parent banks without reapplying TE translation', () => {
  for (const streamwiseFactor of [1, 2]) {
    const { input, system } = coupledFixture(), before = system.evaluate(system.initial);
    const refined = refineCoupledStreamtubeBody(input, system, { streamwiseFactor, normalFactor: 2 });
    assert.equal(refined.input.wakeDisplacementMotion, 'te-center');
    const after = refined.system.evaluate(refined.system.initial);
    for (let g = 0; g < before.outer.nodes.length; g++) for (let i = 0; i < before.outer.nodes[g].length; i++)
      for (let j = 0; j < before.outer.nodes[g][i].length; j++) {
        const p = before.outer.nodes[g][i][j], q = after.outer.nodes[g][streamwiseFactor * i][2 * j];
        close(p.x, q.x); close(p.y, q.y);
      }
  }
});

test('TE-center displacement coordinates transform with rigid rotation and dimensional scaling', () => {
  const input = fixture(), original = createStreamtubeBodySystem(input), value = original.decode(original.initial);
  const angle = .63, factor = 2.7, c = Math.cos(angle), s = Math.sin(angle);
  const point = p => ({ x: .21 + factor * (c * p.x - s * p.y), y: -.37 + factor * (s * p.x + c * p.y) });
  const mapped = { ...input, alpha: input.alpha + angle * 180 / Math.PI,
    bodies: input.bodies.map(b => ({ ...b, points: b.points.map(point), stagnationParameter: factor * b.stagnationParameter })),
    outerLower: input.outerLower.map(point), outerUpper: input.outerUpper.map(point), cutPaths: input.cutPaths.map(row => row.map(point)),
    displacement: { surfaces: input.displacement.surfaces.map(b => Object.fromEntries(Object.entries(b).map(([side, row]) => [side, row.map(d => factor * d)]))),
      wakes: input.displacement.wakes.map(row => row.map(d => factor * d)) } };
  const target = createStreamtubeBodySystem(mapped), actual = target.decode(target.initial).nodes;
  const expected = value.nodes.map(group => group.map(row => row.map(point)));
  nodesClose(actual, expected, 2e-12);
});

test('same-flow warm restart preserves translated physical nodes and rejects a changed chart policy', () => {
  const input = fixture(), source = createStreamtubeBodySystem(input), x = source.initial;
  const before = source.evaluate(x), target = initializeStreamtubeBodyFromFlow(input, source, { initial: x });
  nodesClose(target.flow.nodes, before.nodes);
  assert.ok(max(target.flow.residual.map((v, k) => v - before.residual[k])) < 2e-12);
  assert.throws(() => initializeStreamtubeBodyFromFlow({ ...input, wakeDisplacementMotion: 'fixed' }, source, { initial: x }),
    /wakeDisplacementMotion/);
});
