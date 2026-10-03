import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { respondToCoupledProjectionGeometry as respond } from '../src/euler/streamtube-coupled-projection-geometry.js';

function fixture() {
  const stations = [
    { id: 0, kind: 'surface', body: 0, side: 'upper', i: 1 },
    { id: 1, kind: 'surface', body: 1, side: 'lower', i: 1 },
    { id: 2, kind: 'surface', body: 0, side: 'upper', i: 2 },
    { id: 3, kind: 'wake', body: 0, i: 3, regime: 'wake' },
  ];
  const x = Float64Array.from([7, 11, ...stations.flatMap((_, id) => [.03, .4, [1.3, 1.7, 1.2, .8][id], 1])]);
  const base = Array.from({ length: 3 }, (_, g) => Array.from({ length: 4 }, (_, i) =>
    Array.from({ length: 4 }, (_, j) => ({ x: i, y: 10 * g + 2 * j }))));
  const masses = [[2, 4, 2], [1, 3, 6], [3, 2, 5]];
  let displacement = Array.from(x.slice(2)); const calls = { decode: 0, thicknesses: 0, set: 0 };
  const system = { ne: 2, n: x.length, bl: { stations, phase: [2, 2],
    thicknesses: v => { calls.thicknesses++; return Array.from(v); },
    evaluate: () => { throw new Error('Raw BL must never be evaluated.'); } },
    euler: { layout: { bodies: [{ trailingIndex: 2 }, { trailingIndex: 2 }], independentWakeBanks: true },
      setDisplacement: value => { calls.set++; displacement = structuredClone(value); },
      decode: proposed => {
        calls.decode++; assert.deepEqual(Array.from(proposed), [7, 11]);
        const nodes = structuredClone(base);
        nodes[1][1][0].y += displacement[2]; nodes[1][1][3].y -= displacement[6];
        nodes[1][2][0].y += displacement[10]; // TE moves; independent wake does not.
        return { nodes, allocation: { groups: masses.map(g => g.map(massFlow => ({ massFlow }))) } };
      } } };
  const decoded = system.euler.decode(x.slice(0, 2)); calls.decode = 0;
  const changes = stations.map((station, id) => ({ ...station, deltaStar: x[2 + 4 * id + 2],
    beforeDeltaStar: [1, 1.2, 1, .7][id] }));
  return { system, x, decoded, changes, calls, displacement: () => displacement };
}

test('shared-passage opposite-bank and TE increments use proposed nonuniform physical masses', () => {
  const f = fixture(), before = { x: Array.from(f.x), decoded: structuredClone(f.decoded), phase: [...f.system.bl.phase] };
  const projection = { displacementChanges: f.changes }, unchangedProjection = structuredClone(projection);
  const result = respond(f.system, f.x, projection, f.decoded);
  const close = (a, b) => assert.ok(Math.abs(a - b) < 2e-15, `${a} != ${b}`);
  // Shared passage eta=[0,.1,.4,1], not index thirds. Both endpoints move.
  close(result.nodes[1][1][1].y - f.decoded.nodes[1][1][1].y, .9 * .3 - .1 * .5);
  close(result.nodes[1][1][2].y - f.decoded.nodes[1][1][2].y, .6 * .3 - .4 * .5);
  close(result.nodes[1][2][1].y - f.decoded.nodes[1][2][1].y, .9 * .2);
  for (let g = 0; g < 3; g++) for (let i = 0; i < 4; i++)
    for (const j of [0, 3]) assert.deepEqual(result.nodes[g][i][j], before.decoded.nodes[g][i][j]);
  assert.deepEqual(result.nodes[0], before.decoded.nodes[0]); assert.deepEqual(result.nodes[2], before.decoded.nodes[2]);
  assert.deepEqual(result.nodes[1][3], before.decoded.nodes[1][3]);
  assert.equal(result.diagnostics.surfaceCorrections, 3); assert.equal(result.diagnostics.surfaceTrailingEdgeCorrections, 1);
  assert.equal(result.diagnostics.wakeCorrections, 1); assert.equal(result.diagnostics.rawToProjectedBoundaryNodes, 3);
  assert.equal(result.diagnostics.changedInteriorNodes, 4); assert.equal(result.diagnostics.geometryDecodes, 1);
  assert.equal(result.diagnostics.boundaryExtensions, 1); assert.equal(result.diagnostics.rawBLAccepted, false);
  assert.deepEqual(f.calls, { decode: 1, thicknesses: 2, set: 2 });
  assert.deepEqual(f.displacement(), before.x.slice(2)); assert.deepEqual(Array.from(f.x), before.x);
  assert.deepEqual(f.decoded, before.decoded); assert.deepEqual(f.system.bl.phase, before.phase);
  assert.deepEqual(projection, unchangedProjection);
});

test('inactive and auxiliary-only projections are exact no-op paths without system access', () => {
  const nodes = [[[{ x: 0, y: 0 }]]];
  for (const projection of [undefined, { active: false }, { displacementChanges: [], auxiliaryChanges: [{ id: 3, after: .25 }] }]) {
    const result = respond(null, null, projection, { nodes });
    assert.equal(result.nodes, nodes); assert.equal(result.diagnostics.active, false);
    assert.equal(result.diagnostics.geometryDecodes, 0); assert.equal(result.diagnostics.boundaryExtensions, 0);
  }
});

test('independent wake delta corrections are counted without inventing bank movement', () => {
  const f = fixture(), result = respond(f.system, f.x, { displacementChanges: [f.changes[3]] }, f.decoded);
  assert.equal(result.nodes, f.decoded.nodes); assert.equal(result.diagnostics.wakeCorrections, 1);
  assert.equal(result.diagnostics.independentWakeBanks, true); assert.equal(result.diagnostics.geometryDecodes, 1);
  assert.equal(result.diagnostics.rawToProjectedBoundaryNodes, 0); assert.equal(result.diagnostics.boundaryExtensions, 0);
  assert.deepEqual(f.displacement(), Array.from(f.x.slice(2)));
});

test('invalid projection records reject before mutation and raw decode failure restores projected displacement', () => {
  for (const change of [{ id: 99 }, { deltaStar: 1.31 }, { beforeDeltaStar: -1 }, { side: 'lower' }]) {
    const f = fixture();
    assert.throws(() => respond(f.system, f.x, { displacementChanges: [{ ...f.changes[0], ...change }] }, f.decoded), /Projection geometry/);
    assert.deepEqual(f.calls, { decode: 0, thicknesses: 0, set: 0 });
  }
  const duplicate = fixture();
  assert.throws(() => respond(duplicate.system, duplicate.x,
    { displacementChanges: [duplicate.changes[0], duplicate.changes[0]] }, duplicate.decoded), /duplicate/);
  assert.equal(duplicate.calls.set, 0);
  const f = fixture(), error = Object.freeze(new Error('manufactured raw geometry failure'));
  f.system.euler.decode = () => { throw error; };
  assert.throws(() => respond(f.system, f.x, { displacementChanges: f.changes }, f.decoded), e => e === error);
  assert.deepEqual(f.displacement(), Array.from(f.x.slice(2))); assert.equal(f.calls.set, 2);
});

test('a changed raw interior is rejected after restoring projected displacement', () => {
  const f = fixture(), decode = f.system.euler.decode;
  f.system.euler.decode = x => { const result = decode(x); result.nodes[1][1][1].x += .01; return result; };
  assert.throws(() => respond(f.system, f.x, { displacementChanges: f.changes }, f.decoded), /free interior/);
  assert.deepEqual(f.displacement(), Array.from(f.x.slice(2)));
});

test('the retained exact i104 local projection response reproduces its independently qualified geometry', () => {
  const f = JSON.parse(fs.readFileSync(new URL('./fixtures/projection-geometry-i104.json', import.meta.url)));
  const stations = Array.from({ length: 40 }, (_, id) => ({ id, kind: 'surface', body: 0, side: 'upper', i: 65 + id }));
  const x = new Float64Array(1 + 4 * stations.length); x[0] = 9; x[1 + 4 * 39 + 2] = f.change.deltaStar;
  let thickness = Array.from(x.slice(1)), calls = 0;
  const system = { n: x.length, ne: 1, bl: { stations, thicknesses: x => Array.from(x) },
    euler: { layout: { bodies: [{ trailingIndex: 127 }], independentWakeBanks: true },
      setDisplacement: t => { thickness = t; },
      decode: euler => { calls++; assert.equal(euler[0], 9); assert.equal(thickness[4 * 39 + 2], f.change.beforeDeltaStar);
        return { nodes: [structuredClone(f.raw)] }; } } };
  const decoded = { nodes: [structuredClone(f.projected)], allocation: { groups: [f.masses.map(massFlow => ({ massFlow }))] } };
  const response = respond(system, x, { displacementChanges: [f.change] }, decoded);
  assert.deepEqual(response.nodes[0], f.extended); assert.equal(calls, 1);
  assert.equal(response.diagnostics.changedInteriorNodes, 13);
  assert.deepEqual(thickness, Array.from(x.slice(1)));
});
