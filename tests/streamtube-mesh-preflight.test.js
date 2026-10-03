import test from 'node:test';
import assert from 'node:assert/strict';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { createInitialStreamtubeTopology } from '../src/geometry/streamtube-topology.js';
import { createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { prepareStreamtubeMesh, streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { streamtubeMeshConnectivity } from '../src/geometry/streamtube-mesh-connectivity.js';
import { initializeStreamtubeDensities } from '../src/euler/streamtube-initial-state.js';

test('current coarse main/flap preflight preserves failed builds and repairs both geometry and gas-section feasibility', () => {
  const input = { elements: [{ points: naca4('2412', 80) },
    { points: transform(naca4('0012', 80), { chord: .3, x: .94, y: -.08, angle: -15 }) }], alpha: 4 };
  const prepared = createPanelStreamtubeGrid(createInitialStreamtubeTopology(input, { tubes: 11 }));
  const raw = streamtubeMeshSnapshot(prepared), chart = prepared.system.geometryChart(), before = structuredClone(prepared.nodes);
  assert.equal(raw.quality.valid, false); assert.equal(raw.quality.invalidCells.length, 1);
  const rawEvents = [], unchecked = prepareStreamtubeMesh(prepared, {
    gridRepair: false, onMesh: (mesh, stage) => rawEvents.push({ mesh, stage }),
  });
  assert.deepEqual(rawEvents.map(e => e.stage), ['initial']);
  assert.equal(unchecked.diagnostics.gridRepair, undefined);
  assert.equal(unchecked.diagnostics.gridSmoothing, undefined);
  assert.deepEqual(unchecked.mesh.vertices, raw.vertices);
  assert.deepEqual(prepared.nodes, before);
  assert.deepEqual(prepared.system.geometryChart(), chart);
  const failed = prepareStreamtubeMesh(prepared, { gridRepair: { maxSweeps: 0 } });
  assert.equal(failed.mesh.quality.valid, false); assert.equal(failed.diagnostics.gridRepair.converged, false);
  assert.deepEqual(prepared.nodes, before); assert.deepEqual(prepared.system.geometryChart(), chart);
  assert.deepEqual(failed.mesh.vertices, raw.vertices);

  const events = [], r = prepareStreamtubeMesh(prepared, { onMesh: (mesh, stage) => events.push({ mesh, stage }) });
  assert.deepEqual(events.map(e => e.stage), ['initial', 'initial']);
  assert.equal(events[0].mesh.quality.valid, false); assert.equal(events[1].mesh.quality.valid, true);
  assert.equal(r.mesh.initialization.flowSolved, false); assert.equal(r.diagnostics.gridRepair.converged, true);
  assert.deepEqual(r.mesh.cells, raw.cells);
  const original = streamtubeMeshConnectivity(prepared.system.layout, before);
  original.fixed.forEach(i => assert.deepEqual(r.mesh.vertices[i], original.vertices[i]));
  assert.ok(r.diagnostics.gridRepair.maxDisplacement < .001);
  assert.ok(r.diagnostics.tracedMaxStreamfunctionDrift < 1e-7);
  assert.ok(Number.isFinite(r.diagnostics.maxStreamfunctionDrift));
  assert.ok(r.diagnostics.maxStreamfunctionDrift > r.diagnostics.tracedMaxStreamfunctionDrift, 'repaired nodes must not be presented as unchanged panel traces');
  for (const c of r.mesh.cells) {
    const p = c.map(i => r.mesh.vertices[i]);
    for (let i = 0; i < 4; i++) {
      const a = p[i], b = p[(i + 1) % 4], d = p[(i + 2) % 4];
      assert.ok((b.x - a.x) * (d.y - a.y) > (b.y - a.y) * (d.x - a.x));
    }
  }
  const gasInitial = initializeStreamtubeDensities(r.system, r.initial), gas = r.system.evaluate(gasInitial);
  assert.ok(gas.diagnostics.maxMach < .6); assert.ok(gas.diagnostics.residual > 1e-3, 'initialization must not be relabeled as a converged flow');
  assert.ok(gas.diagnostics.maxStagnationPressureError < 2e-13);
  const repairedChart = r.system.geometryChart(), unchanged = prepareStreamtubeMesh(r);
  assert.deepEqual(unchanged.nodes, r.nodes); assert.deepEqual(unchanged.initial, r.initial);
  assert.deepEqual(r.system.geometryChart(), repairedChart);
});
