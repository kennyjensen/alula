import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { streamtubeAuditRegions, measureStreamtubeGridSpacing, auditStreamtubeMesh } from '../src/geometry/streamtube-grid-audit.js';

function meshFor(xs) {
  const vertices = xs.flatMap(x => [0, .25, .5, .75, 1].map(y => ({ x, y }))), cells = [];
  for (let i = 0; i < xs.length - 1; i++) for (let j = 0; j < 4; j++) cells.push([5 * i + j, 5 * (i + 1) + j, 5 * (i + 1) + j + 1, 5 * i + j + 1]);
  return { topology: 'intrinsic-quadrilateral-streamtubes', vertices, cells,
    initialization: { streamwiseSegments: xs.length - 1, tubes: [4], massFlows: [[.25, .25, .25, .25]], flowSolved: false } };
}

test('recorded shared Xi metadata reaches every independent harmonic audit region', () => {
  const xs = [0, .1, .3, .6, 1], mesh = meshFor(xs), offset = mesh.vertices.length;
  const upperVertices = mesh.vertices.map(p => ({ x: p.x, y: p.y + 1 }));
  const upperCells = mesh.cells.map(cell => cell.map(i => i + offset));
  mesh.vertices.push(...upperVertices); mesh.cells.push(...upperCells);
  mesh.initialization.tubes.push(4); mesh.initialization.massFlows.push([.25, .25, .25, .25]);
  mesh.initialization.gridSmoothing = { stationCoordinate: { source: 'prescribed test Xi', xi: xs.slice() } };
  const saved = structuredClone(mesh), regions = streamtubeAuditRegions(mesh);
  assert.equal(regions.length, 2);
  for (const region of regions) {
    assert.deepEqual(region.streamwiseCoordinates, xs);
    assert.notEqual(region.streamwiseCoordinates, mesh.initialization.gridSmoothing.stationCoordinate.xi);
  }
  const result = auditStreamtubeMesh(JSON.parse(JSON.stringify(mesh)));
  for (const harmonic of result.harmonic) {
    assert.equal(harmonic.passed, true);
    assert.ok(harmonic.levels.at(-1).maximum.crosslineIntervals < 1e-12);
    assert.deepEqual(harmonic.levels.at(-1).streamwiseCoordinates, xs);
  }
  assert.equal(result.accepted, false);
  assert.deepEqual(mesh, saved);
  mesh.initialization.gridSmoothing.stationCoordinate.xi[1] = .12;
  assert.equal(regions[0].streamwiseCoordinates[1], .1);
});

test('recorded Xi metadata fails explicitly when malformed and is optional for old meshes', () => {
  const mesh = meshFor([0, .25, .5, .75, 1]);
  assert.equal(streamtubeAuditRegions(mesh)[0].streamwiseCoordinates, undefined);
  for (const xi of [[0, .2, 1], [0, .25, .25, .75, 1], [0, .25, NaN, .75, 1], [.01, .25, .5, .75, 1], null]) {
    mesh.initialization.gridSmoothing = { stationCoordinate: { xi } };
    assert.throws(() => streamtubeAuditRegions(mesh), /Recorded streamwise coordinates/);
  }
});

test('recorded horizontal farfield conditions reach the reference without imposing false Dirichlet Xi values', () => {
  const xs = [0, .1, .4, .7, 1], mesh = meshFor(xs);
  const boundaryConditions = { lower: 'giles-vertical', upper: 'giles-vertical' };
  mesh.initialization.gridSmoothing = { regions: [{ coordinateEquations: { boundaryConditions } }] };
  const saved = structuredClone(mesh), regions = streamtubeAuditRegions(mesh);
  assert.deepEqual(regions[0].boundaryConditions, boundaryConditions);
  assert.notEqual(regions[0].boundaryConditions, boundaryConditions);
  const spacing = measureStreamtubeGridSpacing(regions);
  assert.equal(spacing.spacingLocations.filter(p => p.fixedBoundary).length, 0);
  const report = auditStreamtubeMesh(JSON.parse(JSON.stringify(mesh)));
  // On this rectangle, xi=x and eta=y solve the mixed physical problem
  // exactly. Assigned uniform Xi therefore differ by .15 at x=.1, including
  // the free boundary. Dividing by the local Xi interval .25 gives .6.
  for (const level of report.harmonic[0].levels) {
    assert.deepEqual(level.boundaryConditions, boundaryConditions);
    assert.ok(level.unknownsByCoordinate.xi > level.unknownsByCoordinate.eta);
    assert.ok(Math.abs(level.maximum.crosslineIntervals - .6) < 1e-12);
    assert.ok(level.maximum.tubeIntervals < 1e-12);
  }
  assert.equal(report.accepted, false);
  assert.deepEqual(mesh, saved);
  boundaryConditions.lower = 'fixed';
  assert.equal(regions[0].boundaryConditions.lower, 'giles-vertical');
  const old = meshFor(xs);
  assert.deepEqual(streamtubeAuditRegions(old)[0].boundaryConditions, { lower: 'fixed', upper: 'fixed' });
  for (const bad of [[], [{ coordinateEquations: { boundaryConditions: null } }],
    [{ coordinateEquations: { boundaryConditions: { lower: 'unknown' } } }]]) {
    old.initialization.gridSmoothing = { regions: bad };
    assert.throws(() => streamtubeAuditRegions(old), /boundary condition/);
  }
});

test('serialized curved-farfield geometry reaches independent reference refinement', () => {
  const mesh = meshFor([0, .25, .5, .75, 1]);
  for (const point of mesh.vertices) point.y += .05 * point.x ** 2;
  const boundaryCurves = Object.fromEntries([['lower', 0], ['upper', 4]].map(([side, j]) => [side, {
    points: Array.from({ length: 5 }, (_, i) => ({ ...mesh.vertices[5 * i + j] })),
    slopes: Array.from({ length: 5 }, (_, i) => .1 * i / 4),
  }]));
  mesh.initialization.gridSmoothing = { regions: [{ coordinateEquations: {
    boundaryConditions: { lower: 'normal-curve', upper: 'normal-curve' }, boundaryCurves,
  } }] };
  const original = structuredClone(mesh), regions = streamtubeAuditRegions(mesh);
  const report = auditStreamtubeMesh(JSON.parse(JSON.stringify(mesh)), { refinements: [1, 2], maxUnknowns: 100 });
  assert.equal(report.accepted, false);
  for (const level of report.harmonic[0].levels) {
    for (const side of ['lower', 'upper']) {
      assert.deepEqual(level.boundaryCurves[side].points, boundaryCurves[side].points);
      assert.deepEqual(level.boundaryCurves[side].slopes, boundaryCurves[side].slopes);
    }
    assert.match(level.geometryRefinement, /fixed cubic Hermite/);
  }
  assert.deepEqual(mesh, original);
  boundaryCurves.upper.slopes[1] = .7;
  assert.equal(regions[0].boundaryCurves.upper.slopes[1], .025);
});
test('grid audit reconstructs exported cells and detects boundary-imposed spacing jumps without declaring physical validity', () => {
  const mesh = meshFor([0, .1, .4, .7, 1]), copy = structuredClone(mesh), regions = streamtubeAuditRegions(mesh);
  const spacing = measureStreamtubeGridSpacing(regions);
  assert.ok(Math.abs(spacing.maximumAdjacentRatio - 3) < 1e-14); assert.equal(spacing.minimumCornerSine, 1);
  assert.equal(spacing.spacingLocations.length, 5); assert.ok(spacing.flaggedCells.length > 0);
  assert.equal(spacing.spacingLocations.filter(p => p.fixedBoundary).length, 2);
  const r = auditStreamtubeMesh(mesh);
  assert.equal(r.accepted, false); assert.equal(r.status, 'grid checks need attention');
  // A straight rectangular channel has exact streamlines even with a bad
  // cross-line distribution: distinguish geometry screening from physics.
  assert.ok(r.harmonic[0].levels.at(-1).maximum.tubeIntervals < 1e-12);
  assert.ok(r.harmonic[0].levels.at(-1).maximum.crosslineIntervals > .05);
  assert.equal(r.fieldErrors.length, mesh.cells.length); assert.deepEqual(mesh, copy);
  const exact = auditStreamtubeMesh(meshFor([0, .25, .5, .75, 1]));
  assert.equal(exact.harmonic[0].passed, true); assert.equal(exact.spacing.flaggedCells.length, 0); assert.equal(exact.accepted, false);
  const bad = structuredClone(mesh); bad.cells[1][0] = 0;
  assert.throws(() => streamtubeAuditRegions(bad), /Disconnected/);
  assert.throws(() => auditStreamtubeMesh({ ...mesh, initialization: { ...mesh.initialization, flowSolved: true } }), /compressible/);
});

test('independent reference rejects the retained positive, SLOR-converged default grid', t => {
  const fixture = JSON.parse(readFileSync(new URL('./fixtures/streamtube-default-slor-audit.json', import.meta.url)));
  const r = auditStreamtubeMesh(fixture.mesh);
  assert.equal(r.accepted, false); assert.equal(r.harmonic.length, 3);
  assert.ok(r.spacing.maximumAdjacentRatio > 9); assert.ok(r.spacing.minimumCornerSine > 0);
  assert.ok(r.harmonic[0].levels.at(-1).maximum.tubeIntervals > .9);
  assert.ok(r.harmonic.every(h => !h.passed)); assert.ok(r.harmonic.every(h => !h.referenceResolved));
  assert.ok(r.harmonic.every(h => h.levels.every(l => Object.values(l.linear).every(v => v.relativeResidual < 1e-10))));
  t.diagnostic(JSON.stringify({ flagged: r.spacing.flaggedCells.length, maximumAdjacentRatio: r.spacing.maximumAdjacentRatio,
    mismatch: r.harmonic.map(h => h.levels.at(-1).maximum) }));
});
