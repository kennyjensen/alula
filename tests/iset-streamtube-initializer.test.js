import test from 'node:test';
import assert from 'node:assert/strict';
import { createIsetStreamtubeGrid } from '../src/euler/tests/iset-initializer.js';
import { naca4 } from '../src/geometry/airfoil.js';
import { initializeStreamtubeDensities } from '../src/euler/streamtube-initial-state.js';

const fixture = () => ({ points: naca4('0012', 40), inletStations: [0, .3, .7, .9, .97, 1],
  outletStations: [0, .03, .1, .3, .7, 1],
  surfaceStations: Object.fromEntries(['upper', 'lower'].map(side => [side,
    Array.from({ length: 9 }, (_, i) => .5 * (1 - Math.cos(Math.PI * i / 8)))])),
  transverseStations: [0, .03, .15, .5, .5, .85, .97, 1] });

test('the ISET seed and RESPLI walls connect to the actual quadrilateral chart and mesh preview', () => {
  const input = fixture(), before = structuredClone(input), notifications = [];
  const r = createIsetStreamtubeGrid(input, { onMesh: (mesh, stage) => notifications.push({ mesh, stage }) });
  assert.deepEqual(input, before); assert.equal(r.system.layout.elements, 1);
  assert.equal(r.mesh.cells.length, 108); assert.equal(r.mesh.quality.valid, true);
  assert.equal(r.diagnostics.gridSmoothing.converged, true, r.diagnostics.gridSmoothing.reason);
  assert.equal(r.diagnostics.gridSmoothing.farfieldBoundary, 'giles-indexed-y');
  assert.deepEqual(r.diagnostics.gridSmoothing.stationCoordinate.xi, r.sourceInitialization.xi);
  assert.equal(r.mesh.initialization.flowSolved, false);
  assert.equal(r.mesh.initialization.iset.originalFortranExecuted, false);
  assert.ok(notifications.some(n => n.stage === 'smoothing'));
  assert.equal(notifications[0].stage, 'initial');
  const restored = r.system.decode(r.initial).nodes;
  restored.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    assert.ok(Math.hypot(p.x - r.nodes[g][i][j].x, p.y - r.nodes[g][i][j].y) < 2e-13);
  })));
  // This is a topology check independent of the source's logical arrays:
  // one closed hole, no repeated corners or disconnected cut banks.
  const edges = new Map();
  for (const cell of r.mesh.cells) {
    assert.equal(new Set(cell).size, 4);
    for (let k = 0; k < 4; k++) { const a = cell[k], b = cell[(k + 1) % 4], key = a < b ? `${a}:${b}` : `${b}:${a}`;
      edges.set(key, (edges.get(key) ?? 0) + 1); }
  }
  assert.equal(r.mesh.vertices.length - edges.size + r.mesh.cells.length, 0);
  assert.ok([...edges.values()].every(n => n === 1 || n === 2));
  // Enter the existing target-Mach cell equations without a Newton update.
  const gas = initializeStreamtubeDensities(r.system, r.initial), evaluated = r.system.evaluate(gas);
  assert.ok([...evaluated.residual].every(Number.isFinite));
  assert.equal(r.system.conditions.mach, .2);
});

test('ISET reports an exhausted smoother and preserves the source seed without claiming a flow solution', () => {
  const r = createIsetStreamtubeGrid(fixture(), { smoothing: { maxSweeps: 0 } });
  assert.equal(r.diagnostics.gridSmoothing.converged, false);
  assert.equal(r.diagnostics.gridSmoothing.retainedOriginal, true);
  assert.equal(r.mesh.initialization.flowSolved, false);
  const raw = createIsetStreamtubeGrid(fixture(), { smoothing: false });
  assert.deepEqual(r.nodes, raw.nodes);
  assert.equal(raw.diagnostics.gridSmoothing, undefined);
  assert.throws(() => createIsetStreamtubeGrid({ ...fixture(), surfaceStations: { upper: [0, .5, 1], lower: [0, .5, 1] } }), /at least five/);
  assert.throws(() => createIsetStreamtubeGrid(fixture(), { smoothing: [] }), /SLOR controls/);
});
