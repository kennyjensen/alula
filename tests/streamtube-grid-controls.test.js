import test from 'node:test';
import assert from 'node:assert/strict';
import { createInitialStreamtubeTopology } from '../src/geometry/streamtube-topology.js';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { solveStreamtubeAssembly } from '../src/euler/streamtube-result.js';
import { streamtubeAuditRegions } from '../src/geometry/streamtube-grid-audit.js';

test('separate upper/lower/gap counts preserve mass profiles and chord scaling changes the requested surface allocation', () => {
  const input = { elements: [{ points: naca4('0012', 80) }, { points: transform(naca4('0012', 80), { chord: .3, x: .94, y: -.08, angle: -15 }) }] };
  const controls = { surfaceIntervals: 32, lowerTubes: 5, gapTubes: 7, upperTubes: 11, surfaceChordExponent: 1, inletIntervals: 16, outletIntervals: 32 };
  const r = createInitialStreamtubeTopology(input, controls);
  assert.deepEqual(r.weights.map(row => row.length), [5, 7, 11]);
  assert.deepEqual(r.gridSpacing.requestedSurfaceIntervalsByElement.map(p => p.intervals).sort((a, b) => b - a), [32, 10]);
  for (const [g, n] of [[0, 5], [2, 11]]) {
    const sum = r.weights[g].reduce((s, w) => s + w, 0), growth = 3 ** (7 / n);
    r.weights[g].forEach((w, j) => {
      const power = g === 0 ? n - 1 - j : j;
      const exact = growth ** power * (growth - 1) / (3 ** 7 - 1);
      assert.ok(Math.abs(w / sum - exact) < 2e-15);
    });
  }
  assert.equal(Math.min(...r.bodies.map(b => b.leadingIndex)), 16);
  assert.equal(r.outerLower.length - 1 - Math.max(...r.bodies.map(b => b.trailingIndex)), 32);
  const rotated = { elements: input.elements.map(e => ({ points: transform(e.points, { chord: 2, angle: 7, x: 3, y: 1 }) })) };
  assert.deepEqual(createInitialStreamtubeTopology(rotated, controls).gridSpacing.requestedSurfaceIntervalsByElement,
    r.gridSpacing.requestedSurfaceIntervalsByElement);
  for (const extra of [{ upperTubes: 2 }, { gapTubes: 32 }, { surfaceChordExponent: -1 }, { surfaceChordExponent: 2 }])
    assert.throws(() => createInitialStreamtubeTopology(input, { ...controls, ...extra }), /controls/);
});

test('exhausted public SLOR initialization retains the original mesh and records its rejected candidate without entering Euler', () => {
  const meshes = [];
  assert.throws(() => solveStreamtubeAssembly({ elements: [{ points: naca4('0012', 80) }], alpha: 0,
    gridIntervals: 8, gridEllipticSmoothing: { maxSweeps: 0 } }, {
    onMesh: mesh => meshes.push(mesh), onIteration: () => assert.fail('An unrelaxed grid must not enter Euler'),
  }), /Elliptic grid initialization failed: sweep limit/);
  assert.ok(meshes.length > 1);
  assert.equal(meshes.at(-1).initialization.gridSmoothing.converged, false);
  assert.equal(meshes.at(-1).initialization.flowSolved, false);
  assert.equal(meshes.at(-1).initialization.gridSmoothing.retainedOriginal, true);
  assert.ok(meshes.at(-1).initialization.gridSmoothing.rejectedMesh.vertices.length > 0);
  assert.deepEqual(meshes.at(-1).vertices, meshes[0].vertices);
  assert.deepEqual(meshes.at(-1).cells, meshes[0].cells);
  assert.ok(meshes.at(-1).vertices.length > 0);
});

test('public curvature/SLOR controls build the actual quad mesh and publish geometry before and during smoothing', t => {
  const events = [], r = solveStreamtubeAssembly({ elements: [{ points: naca4('0012', 80) }], alpha: 0,
    gridIntervals: 16, gridLowerTubes: 7, gridUpperTubes: 9, gridInletIntervals: 16, gridOutletIntervals: 32,
    gridSurfaceSpacing: 'curvature', gridCurvatureSpacing: { exponent: .5, leadingSpacingRatio: .2, trailingSpacingRatio: .4 },
    gridEllipticSmoothing: true, gridCrosslinePlacement: 'potential' }, { meshOnly: true, onMesh: (mesh, stage) => events.push({ mesh, stage }),
    onIteration: () => assert.fail('Build mesh must not run Euler iterations') });
  assert.equal(r.status, 'mesh-ready'); assert.equal(events[0].stage, 'initial');
  assert.ok(events.some(e => e.stage === 'smoothing')); assert.equal(events.at(-1).stage, 'initial');
  assert.equal(r.mesh.quality.valid, true); assert.equal(r.mesh.initialization.flowSolved, false);
  assert.deepEqual(r.mesh.initialization.tubes, [8, 10]);
  assert.equal(r.mesh.initialization.gridSmoothing.converged, true);
  // Check every region boundary, including shared internal cuts, against
  // the actual base mesh published before smoothing. Interior nodes must move.
  const original = events[0].mesh, boundary = new Set();
  for (const { indices } of streamtubeAuditRegions(original))
    indices.forEach((row, i) => row.forEach((id, j) => {
      if (!i || i === indices.length - 1 || !j || j === row.length - 1) boundary.add(id);
    }));
  for (const { mesh } of events.filter(e => !e.mesh.initialization.gridRefinement)) {
    assert.deepEqual(mesh.cells, original.cells);
    for (const id of boundary) assert.deepEqual(mesh.vertices[id], original.vertices[id]);
  }
  const smoothed = events.filter(e => !e.mesh.initialization.gridRefinement).at(-1).mesh;
  assert.ok(smoothed.vertices.some((p, id) => !boundary.has(id)
    && (p.x !== original.vertices[id].x || p.y !== original.vertices[id].y)));
  const banks = mesh => streamtubeAuditRegions(mesh).map(r => r.nodes.map(row => [row[0], row.at(-1)]));
  assert.deepEqual(banks(r.mesh), banks(smoothed), 'Wall subdivision must retain every material bank');
  assert.equal(r.mesh.initialization.surfaceDistributions.length, 2);
  for (const d of r.mesh.initialization.surfaceDistributions) {
    assert.ok(Math.abs(d.preMatching.actualLeadingSpacingRatio - .2) < 1e-8);
    assert.ok(Math.abs(d.preMatching.actualTrailingSpacingRatio - .4) < 1e-8);
    assert.ok(d.actualLeadingSpacingRatio > 0 && d.actualTrailingSpacingRatio > 0);
    assert.ok(d.count >= d.preMatching.count);
  }
  assert.ok(r.mesh.cells.every(c => c.length === 4));
  t.diagnostic(JSON.stringify({ cells: r.mesh.cells.length, quality: r.mesh.quality,
    regions: r.mesh.initialization.gridSmoothing.regions.map(q => q.history.at(-1)) }));
});
