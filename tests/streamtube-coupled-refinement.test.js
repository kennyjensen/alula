import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createCoupledStreamtubeBody, solveCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { refineCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-refinement.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directChannelConservation } from './oracles/streamtube.js';

const retained = () => {
  const f = JSON.parse(fs.readFileSync(new URL('fixtures/streamtube-isentropic-coupled-root.json', import.meta.url)));
  return { input: f.input, system: createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL }) };
};
function checkTransfer(input, parent, initial = parent.initial, normalFactor = 2) {
  const before = structuredClone(input), x = initial.slice(), old = parent.evaluate(initial);
  const refined = refineCoupledStreamtubeBody(input, parent, { initial, normalFactor }), child = refined.system, v = child.evaluate(child.initial);
  assert.deepEqual(input, before); assert.deepEqual(initial, x); assert.equal(child.admissible(child.initial), true);
  assert.equal(refined.diagnostics.quality.valid, true); assert.equal(child.bl.surfaces.length, parent.bl.surfaces.length);
  assert.equal(child.bl.wakes.length, parent.bl.wakes.length);
  assert.equal(child.conditions.referenceChord, parent.conditions.referenceChord);
  assert.equal(child.conditions.reynolds, parent.conditions.reynolds);
  assert.equal(child.conditions.edgeMatching, parent.conditions.edgeMatching);
  for (let g = 0; g < old.outer.nodes.length; g++) {
    old.outer.nodes[g].forEach((row, i) => row.forEach((p, j) => {
      const q = v.outer.nodes[g][2 * i][normalFactor * j];
      assert.ok(Math.hypot(q.x - p.x, q.y - p.y) < 2e-12, `physical parent node ${g}/${i}/${j}`);
    }));
    old.outer.allocation.groups[g].forEach((q, j) => {
      const masses = v.outer.allocation.groups[g].slice(normalFactor * j, normalFactor * (j + 1)).map(t => t.massFlow);
      assert.ok(masses.every(m => Math.abs(m / q.massFlow - 1 / normalFactor) < 1e-14));
      assert.ok(Math.abs(masses.reduce((a, b) => a + b, 0) / q.massFlow - 1) < 1e-14);
    });
  }
  for (const branch of parent.bl.surfaces) {
    const next = child.bl.surfaces.find(b => b.body === branch.body && b.side === branch.side);
    assert.equal(next.tripParameter, branch.tripParameter);
    const oldFirst = old.layers.states[branch.ids[0]], newFirst = v.layers.states[next.ids[0]];
    assert.ok(newFirst.s / oldFirst.s < .501 && newFirst.s / oldFirst.s > .499);
  }
  for (const station of parent.bl.stations) {
    const next = child.bl.stations.find(s => s.kind === station.kind && s.body === station.body && s.side === station.side && s.i === 2 * station.i);
    assert.ok(next);
    for (const key of ['aux', 'theta', 'deltaStar', 'ue']) assert.equal(v.layers.states[next.id][key], old.layers.states[station.id][key]);
  }
  for (let b = 0; b < parent.bl.wakes.length; b++) {
    const wake = child.bl.wakes[b];
    for (const id of wake.ids.slice(1)) {
      const station = child.bl.stations[id], i = station.i, a = v.outer.nodes[b][i].at(-1), q = v.outer.nodes[b + 1][i][0];
      const gap = v.layers.states[id].deltaStar * child.conditions.referenceChord;
      assert.ok(Math.abs(Math.hypot(a.x - q.x, a.y - q.y) - gap) < 1e-12);
    }
  }
  return refined;
}

test('nested attached-root refinement preserves parent fluid nodes, BL values and material trips while halving intervals and masses', t => {
  const { input, system } = retained(), r = checkTransfer(input, system);
  assert.equal(r.system.n, 11929);
  // A prolonged guess has an independently evaluated equation defect and
  // cannot inherit the parent's converged status.
  const stopped = solveCoupledStreamtubeBody(r.system, { maxIterations: 0 });
  assert.equal(stopped.converged, false); assert.ok(stopped.families.euler > 1e-4);
  t.diagnostic(JSON.stringify({ unknowns: r.system.n, families: stopped.families, quality: r.diagnostics.quality }));
});

for (const edgeMatching of ['pressure', 'section-velocity']) test(`${edgeMatching}: nested refinement transfers all four surfaces, two wakes and capture globals on a coupled two-element root`, () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }), parent = createCoupledStreamtubeBody(input, { edgeMatching });
  const root = solveCoupledStreamtubeBody(parent, { tolerance: 1e-10 }); assert.equal(root.converged, true);
  const r = checkTransfer(input, parent, root.x);
  assert.equal(r.system.bl.surfaces.length, 4); assert.equal(r.system.bl.wakes.length, 2);
});

test('stagnation-scaled refinement halves surface intervals and quarters every streamtube mass while retaining parent data', t => {
  const { input, system } = retained(), r = checkTransfer(input, system, system.initial, 4);
  assert.equal(r.diagnostics.normalFactor, 4); assert.equal(r.diagnostics.streamwiseFactor, 2);
  t.diagnostic(JSON.stringify({ unknowns: r.system.n, quality: r.diagnostics.quality, families: r.system.evaluate(r.system.initial).families }));
});

test('refinement rejects invalid factors and node budgets without changing the parent state', () => {
  const { input, system } = retained(), before = system.initial.slice();
  for (const controls of [{ streamwiseFactor: 1, normalFactor: 1 }, { streamwiseFactor: 0 }, { normalFactor: 1.5 }, { maxNodes: 1 },
    { normalSubdivisions: [] }, { normalSubdivisions: [[1], [1]] }, { normalSubdivisions: null },
    { normalFactor: 2, normalSubdivisions: input.weights.map(row => row.map(() => 2)) },
    ...[0, 1.5, 5].map(n => ({ normalSubdivisions: input.weights.map(row => row.map(() => n)) }))])
    assert.throws(() => refineCoupledStreamtubeBody(input, system, controls), /refinement/);
  assert.deepEqual(system.initial, before);
});

test('selective normal refinement conserves each parent mass and preserves all four surface BLs and both wakes', t => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const parent = createCoupledStreamtubeBody(input, { edgeMatching: 'section-velocity' });
  const root = solveCoupledStreamtubeBody(parent, { tolerance: 1e-10 }); assert.equal(root.converged, true);
  const before = structuredClone(input), x = root.x.slice(), old = parent.evaluate(root.x);
  const counts = [[1, 4], [4, 4], [4, 1]];
  const r = refineCoupledStreamtubeBody(input, parent, { initial: root.x, streamwiseFactor: 1, normalSubdivisions: counts });
  const child = r.system, value = child.evaluate(child.initial);
  assert.deepEqual(input, before); assert.deepEqual(root.x, x);
  assert.deepEqual(child.euler.layout.tubes, [5, 8, 5]);
  assert.equal(child.n < 4 * parent.n, true); assert.equal(r.diagnostics.normalFactor, null);
  assert.equal(r.diagnostics.quality.valid, true); assert.equal(child.admissible(child.initial), true);
  for (let g = 0; g < counts.length; g++) {
    const positions = [0]; for (const n of counts[g]) positions.push(positions.at(-1) + n);
    old.outer.nodes[g].forEach((row, i) => row.forEach((p, j) => {
      const q = value.outer.nodes[g][i][positions[j]];
      assert.ok(Math.hypot(q.x - p.x, q.y - p.y) < 2e-12);
    }));
    old.outer.allocation.groups[g].forEach((q, j) => {
      const masses = value.outer.allocation.groups[g].slice(positions[j], positions[j + 1]).map(tube => tube.massFlow);
      assert.ok(masses.every(m => Math.abs(m / q.massFlow - 1 / counts[g][j]) < 1e-14));
      assert.ok(Math.abs(masses.reduce((a, b) => a + b, 0) / q.massFlow - 1) < 1e-14);
    });
  }
  assert.deepEqual(value.layers.states, old.layers.states);
  assert.deepEqual(child.bl.surfaces.map(b => b.tripParameter), parent.bl.surfaces.map(b => b.tripParameter));
  assert.deepEqual(value.outer.captured, old.outer.captured);
  const solved = solveCoupledStreamtubeBody(child, { tolerance: 1e-10, maxIterations: 12 });
  assert.equal(solved.converged, true);
  assert.ok(Object.values(solved.families).every(v => v < 1e-10));
  solved.flow.nodes.forEach((nodes, g) => {
    const c = directChannelConservation({ nodes, sections: solved.flow.sections.map(row => row[g]),
      cells: solved.flow.cells.map(row => row[g]) }, child.euler.conditions.gamma);
    for (const key of ['maxLocal', 'total', 'internalCancellation'])
      assert.ok(c[key].every(v => Math.abs(v) < 2e-9), `independent ${key} balance in passage ${g}`);
  });
  t.diagnostic(JSON.stringify({ parent: parent.n, refined: child.n, iterations: solved.history.length - 1, families: solved.families }));
});
