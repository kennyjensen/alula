import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody, coupledStreamtubeResult } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { independentWakeBankRestart } from '../src/euler/tests/streamtube-independent-wake-restart.js';
import { streamtubeWakeGap } from '../src/euler/streamtube-wake-geometry.js';
import { extendStreamtubeDisplacement } from '../src/euler/streamtube-displacement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { directChannelConservation } from './oracles/streamtube.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

const maximum = a => Math.max(0, ...Array.from(a, Math.abs));
function prescribed() {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  return { ...input, wakeGeometry: 'independent-banks', displacement: {
    surfaces: input.bodies.map(b => Object.fromEntries(['upper', 'lower'].map((s, k) => [s,
      Array.from({ length: b.trailingIndex - b.leadingIndex + 1 }, (_, i) => .0003 * (1 + i * (k + 1) / 10))]))),
    wakes: input.bodies.map(b => Array(input.outerLower.length - 1 - b.trailingIndex).fill(.001)) } };
}
function checkDirections(euler, nodes) {
  const chart = new Map(euler.geometryChart().map(p => [p.column, p.normal]));
  for (const p of euler.layout.positions.filter(p => p.kind === 'cut' && p.side)) {
    const point = i => p.side === 'lower' ? nodes[p.body][i].at(-1) : nodes[p.body + 1][i][0];
    const a = point(p.i - 1), b = point(Math.min(euler.layout.nx, p.i + 1)), n = chart.get(p.column);
    assert.ok(Math.abs((b.x - a.x) * n.x + (b.y - a.y) * n.y) / Math.hypot(b.x - a.x, b.y - a.y) < 1e-12);
  }
}

test('initial displacement extension translates every wake center by the TE shift exactly once', () => {
  for (const wakeGeometry of ['centerline', 'independent-banks']) {
    const system = createStreamtubeBodySystem({ ...prescribed(), wakeGeometry });
    const decoded = system.decode(system.initial), moved = extendStreamtubeDisplacement(system, system.initial);
    const center = (nodes, b, i, key) => .5 * (nodes[b][i].at(-1)[key] + nodes[b + 1][i][0][key]);
    for (let b = 0; b < system.layout.elements; b++) {
      const te = system.layout.bodies[b].trailingIndex;
      for (const key of ['x', 'y']) {
        const delta = center(decoded.nodes, b, te, key) - center(decoded.undisplacedNodes, b, te, key);
        for (let i = te + 1; i <= system.layout.nx; i++) assert.ok(Math.abs(center(moved, b, i, key)
          - center(decoded.undisplacedNodes, b, i, key) - delta) < 2e-15, `${wakeGeometry}, body ${b}, station ${i}, ${key}`);
      }
    }
  }
});

test('thickness changes only the independent gap row; full Euler and thickness derivatives pass', () => {
  const input = prescribed(), system = createStreamtubeBodySystem(input), x = system.initial.map((_, i) => 1e-6 * Math.sin(i));
  checkDirections(system, system.decode(system.initial).nodes);
  for (let chart = 0; chart < 2; chart++) {
    const a = system.evaluate(x), jac = system.jacobian(x, { includeDisplacement: true });
    for (let col = 0; col < system.layout.n; col++) {
      const h = 2e-7, p = x.slice(), m = x.slice(); p[col] += h; m[col] -= h;
      const plus = system.residual(p), minus = system.residual(m);
      for (let row = 0; row < system.layout.n; row++) {
        const d = jac.state[row * system.layout.n + col], fd = (plus[row] - minus[row]) / (2 * h);
        assert.ok(Math.abs(d - fd) / Math.max(1, Math.abs(d), Math.abs(fd)) < 3e-6, `Euler ${row}, ${col}`);
      }
    }
    for (const [col, p] of system.displacementParameters.entries()) {
      const h = 2e-7, values = [1, -1].map(sign => {
        const d = structuredClone(input.displacement);
        if (p.kind === 'wake') d.wakes[p.body][p.index] += sign * h;
        else for (const side of p.side === 'both' ? ['lower', 'upper'] : [p.side]) d.surfaces[p.body][side][p.index] += sign * h;
        system.setDisplacement(d); return system.evaluate(x);
      });
      for (let row = 0; row < system.layout.n; row++) {
        const d = jac.displacement[row].get(col) ?? 0, fd = (values[0].residual[row] - values[1].residual[row]) / (2 * h);
        assert.ok(Math.abs(d - fd) / Math.max(1, Math.abs(d), Math.abs(fd)) < 3e-6, `displacement ${row}, ${col}`);
      }
      if (p.kind === 'wake') assert.deepEqual(values[0].nodes, values[1].nodes);
    }
    system.setDisplacement(input.displacement);
    const next = system.rebase(x); x.set(next);
    assert.ok(maximum(a.residual.map((v, i) => v - system.residual(x)[i])) < 1e-10);
    checkDirections(system, system.decode(x).nodes);
  }
  const invalid = system.decode(x).nodes, body = 0, i = system.layout.bodies[0].trailingIndex + 2;
  const lower = invalid[body][i].at(-1); invalid[body + 1][i][0] = { x: lower.x, y: lower.y - .01 };
  const bad = system.adoptGeometry(x, invalid); assert.equal(system.admissible(bad), false);
});

test('the complete saved default converts without changing its geometry, BL state or existing residuals', t => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/default-coupled-ises-wake-fold.json', import.meta.url)));
  const c = independentWakeBankRestart(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
  assert.equal(c.system.n, 8278); assert.equal(c.diagnostics.addedWakeEquations, 83);
  assert.ok(c.diagnostics.retainedEulerResidualChange < 1e-11); assert.ok(c.diagnostics.boundaryLayerResidualChange < 1e-11);
  assert.ok(c.diagnostics.maxGeometryChange < 1e-13); assert.ok(c.diagnostics.gapResidual < 1e-13);
  assert.equal(c.system.admissible(c.system.initial), false);
  assert.equal(c.system.admissible(c.system.initial, { requireConvex: false }), true);
  checkDirections(c.system.euler, c.system.evaluate(c.system.initial).outer.nodes);
  assert.deepEqual(c.system.bl.activeTargets(c.system.initial.subarray(0, c.system.ne)).filter(p => p.from !== p.to), []);
  t.diagnostic(JSON.stringify(c.diagnostics));
});

test('the complete default root closes all coupled equations but cannot bypass its concave final quad', () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/default-independent-wake-root.json', import.meta.url)));
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: Float64Array.from(f.initialBL) });
  const r = coupledStreamtubeResult(system, system.initial, { tolerance: 1e-8 });
  assert.equal(system.n, 8278); assert.ok(maximum(r.residual) < 1e-8);
  assert.equal(r.converged, false); assert.equal(r.mesh.quality.valid, false); assert.equal(r.mesh.initialization.flowSolved, false);
  assert.deepEqual(r.mesh.quality.invalidCells, [463]); assert.equal(directStreamtubeVolumeGeometry(r.flow.nodes).valid, true);
  assert.equal(r.cl, null); assert.equal(r.cd, null);
});

test('independent banks solve with all four BLs and both wakes, conserve, and retain their complete restart Jacobian', t => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), wakeGeometry: 'independent-banks' };
  const snapshots = [], r = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', maxIterations: 12, tolerance: 1e-10,
    stepAcceptance: 'admissible', onMesh: m => snapshots.push(m.nodes) });
  assert.equal(r.converged, true, r.reason); assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.boundaryLayer.surfaces.length, 4); assert.equal(r.boundaryLayer.wakes.length, 2);
  assert.deepEqual(snapshots.at(-1), r.flow.nodes); assert.equal(directStreamtubeVolumeGeometry(r.flow.nodes).valid, true);
  const ne = r.x.length - 4 * r.boundaryLayer.stations.length;
  const system = createCoupledStreamtubeBody(r.solverInput, { ...r.coupledOptions, initialEuler: r.flow, initialBL: r.x.slice(ne) });
  const x = system.initial, value = system.evaluate(x), jac = system.jacobian(x), errors = [];
  assert.ok(maximum(value.residual) < 1e-10); assert.deepEqual(value.outer.nodes, r.flow.nodes);
  checkDirections(system.euler, value.outer.nodes);
  for (const kind of ['euler', 'bl', 'both']) {
    const d = x.map((v, i) => i < ne ? kind === 'bl' ? 0 : .001 * Math.sin(.73 * i + .2)
      : kind === 'euler' ? 0 : Math.max(.01, Math.abs(v)) * Math.sin(.43 * i + .4));
    const exact = sparseProduct(jac, d), h = 2e-6;
    const a = system.residual(x.map((v, i) => v + h * d[i])), b = system.residual(x.map((v, i) => v - h * d[i]));
    const error = maximum(exact.map((v, i) => { const fd = (a[i] - b[i]) / (2 * h); return (fd - v) / Math.max(1, Math.abs(v), Math.abs(fd)); }));
    assert.ok(error < 5e-6, `${kind}: ${error}`); errors.push({ kind, error });
  }
  system.evaluate(x);
  for (let g = 0; g < r.flow.nodes.length; g++) {
    const c = directChannelConservation({ nodes: r.flow.nodes[g], sections: r.flow.sections.map(row => row[g]), cells: r.flow.cells.map(row => row[g]) });
    for (const key of ['maxLocal', 'total', 'internalCancellation']) assert.ok(maximum(c[key]) < 2e-9, `${g}, ${key}`);
  }
  for (const w of r.boundaryLayer.wakes) for (const id of w.ids.slice(1)) {
    const s = r.boundaryLayer.stations[id], indices = [s.i - 1, s.i, Math.min(system.euler.layout.nx, s.i + 1)];
    const gap = streamtubeWakeGap(indices.map(i => r.flow.nodes[w.body][i].at(-1)), indices.map(i => r.flow.nodes[w.body + 1][i][0])).gap;
    assert.ok(Math.abs(gap - system.euler.conditions.lengthScale * s.deltaStar) < 1e-10);
  }
  t.diagnostic(JSON.stringify({ unknowns: r.x.length, updates: r.history.length - 1, families: r.families, errors }));
});
