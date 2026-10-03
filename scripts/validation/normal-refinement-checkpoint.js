// SPDX-License-Identifier: GPL-2.0-or-later
// A controlled normal-refinement continuation: transfer existing maintenance
// history without invoking a new initial SMOVE. This is a refinement-policy
// experiment, not a claim about an unavailable MSES refinement listing.
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';
const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const create = f => createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });

export function normalRefinementCheckpoint(parentCheckpoint, refinedRestart, normalSubdivisions) {
  assert.equal(parentCheckpoint.version, 1);
  const a = parentCheckpoint.restart, b = serialize(refinedRestart), parent = create(a), child = create(b);
  const before = parent.evaluate(parent.initial), after = child.evaluate(child.initial);
  assert.deepEqual(before.families, parentCheckpoint.families);
  assert.equal(child.euler.layout.nx, parent.euler.layout.nx, 'Only normal refinement can retain this inlet history.');
  // Refinement writes explicit defaults that a serialized cold checkpoint
  // may omit. Compare the effective model and material trips, not spelling.
  assert.deepEqual(child.conditions, parent.conditions); assert.deepEqual(child.bl.trips, parent.bl.trips);
  assert.equal(child.bl.transitionMode, parent.bl.transitionMode);
  assert.deepEqual(child.euler.conditions, parent.euler.conditions);
  assert.deepEqual(child.euler.fractions, parent.euler.fractions);
  assert.deepEqual(child.bl.snapshotActive(), parent.bl.snapshotActive());
  for (const key of new Set([...Object.keys(a.input), ...Object.keys(b.input)]))
    if (!['weights', 'gridSpacing', 'bodies'].includes(key)) assert.deepEqual(b.input[key], a.input[key], `Changed input: ${key}`);
  const bodyGeometry = bodies => bodies.map(({ surfaceFractions, ...body }) => body);
  assert.deepEqual(bodyGeometry(b.input.bodies), bodyGeometry(a.input.bodies));
  assert.equal(normalSubdivisions.length, parent.euler.layout.tubes.length);
  const retained = normalSubdivisions.map((row, g) => {
    assert.equal(row.length, parent.euler.layout.tubes[g]);
    const indices = [0];
    for (const count of row) { assert.ok(Number.isInteger(count) && count >= 1 && count <= 4); indices.push(indices.at(-1) + count); }
    assert.equal(indices.at(-1), child.euler.layout.tubes[g]); return indices;
  });
  let nodeError = 0, massRelativeError = 0, physicalBLError = 0;
  before.outer.nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    const q = after.outer.nodes[g][i][retained[g][j]]; nodeError = Math.max(nodeError, Math.hypot(q.x - p.x, q.y - p.y));
  })));
  before.outer.allocation.groups.forEach((group, g) => group.forEach((tube, j) => {
    const mass = after.outer.allocation.groups[g].slice(retained[g][j], retained[g][j + 1]).reduce((n, t) => n + t.massFlow, 0);
    massRelativeError = Math.max(massRelativeError, Math.abs(mass / tube.massFlow - 1));
  }));
  assert.equal(parent.bl.stations.length, child.bl.stations.length);
  for (const p of parent.bl.stations) {
    assert.deepEqual(child.bl.stations[p.id], p);
    for (const key of ['s', 'aux', 'theta', 'deltaStar', 'ue'])
      physicalBLError = Math.max(physicalBLError, Math.abs(before.layers.states[p.id][key] - after.layers.states[p.id][key]));
  }
  for (const [key, indices] of Object.entries(parent.euler.layout.globals)) {
    const aCols = Array.isArray(indices) ? indices : [indices], target = child.euler.layout.globals[key];
    const bCols = Array.isArray(target) ? target : [target];
    aCols.forEach((col, i) => { if (col !== null) assert.equal(parent.initial[col], child.initial[bCols[i]], `Changed global ${key}`); });
  }
  assert.ok(nodeError < 2e-12 * parent.euler.conditions.lengthScale && massRelativeError < 1e-13 && physicalBLError < 1e-14);
  const checkpoint = { version: 1, families: after.families, restart: b, continuation: serialize(parentCheckpoint.continuation) };
  return { checkpoint, diagnostics: { nodeError, massRelativeError, physicalBLError,
    retainedNormalStations: retained, parentFamilies: before.families, refinedFamilies: after.families,
    maintenanceHistoryPreserved: true, initialSMOVERepeated: false } };
}
