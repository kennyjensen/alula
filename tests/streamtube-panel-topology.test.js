import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { naca4 } from '../src/geometry/airfoil.js';
import { createInitialStreamtubeTopology } from '../src/geometry/streamtube-topology.js';
import { preparePanelStreamtubeTopology, validatePanelTopologySeed } from '../src/euler/streamtube-panel-topology.js';
import { createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { streamfunctionAt, potentialDifference } from '../src/inviscid/streamfunction.js';
import { velocityAt } from '../src/inviscid/linear-vortex.js';

const input = JSON.parse(fs.readFileSync('docs/solver-reliability/gui-defaults-current/30p30n/input.json'));
let prepared;
const recovery = () => prepared ??= preparePanelStreamtubeTopology(input, { surfaceIntervals: 16, tubes: 7 });

test('existing monotone topology returns exactly with no alternate panel solution', () => {
  const single = { elements: [{ points: naca4('0012', 40) }], alpha: 4, mach: .2 };
  assert.deepEqual(preparePanelStreamtubeTopology(single, {}), { topology: createInitialStreamtubeTopology(single) });
});

test('30P30N reentrant branches retain exact geometry and use ordered physical inlet streamfunction cuts', () => {
  const before = JSON.stringify(input);
  assert.throws(() => createInitialStreamtubeTopology(input), error => error.code === 'STREAMTUBE_NONMONOTONE_BODY' && error.diagnostics.side === 'lower');
  const { topology, panelSolution: panel, diagnostics } = recovery();
  assert.deepEqual(topology.bodies.map(b => b.element), [2, 1, 0]);
  assert.equal(topology.primaryBody, 1);
  topology.bodies.forEach(body => {
    assert.deepEqual(body.points, input.elements[body.element].points);
    assert.deepEqual(body.trailingEdge, input.elements[body.element].trailingEdge);
    assert.ok(body.surfaceFractions.every((f, i, row) => !i || f > row[i - 1]));
  });
  assert.ok(topology.captureLevels.every((p, i, row) => !i || p > row[i - 1]));
  diagnostics.inletCuts.forEach((cut, i) => {
    assert.ok(Math.abs(streamfunctionAt(cut, panel.field) - diagnostics.bodyStreamfunctions[i]) < 1e-11);
    assert.ok(velocityAt(cut, panel.field).u > 0);
    assert.deepEqual(topology.cutPaths[i][0], { x: cut.x, y: cut.y });
  });
  assert.equal(diagnostics.panelSolves, 1); assert.equal(diagnostics.panelStrengthsChanged, false);
  assert.equal(diagnostics.solidGeometryChanged, false); assert.equal(JSON.stringify(input), before);
});

test('reindexing preserves the complete shared field and its physical potential differences', () => {
  const { panelSolution: panel } = recovery();
  const old = JSON.parse(fs.readFileSync('docs/solver-reliability/gui-defaults-current/30p30n/recovery-panel-source.json')).panel;
  old.field.gamma = Float64Array.from(Object.values(old.field.gamma));
  assert.deepEqual(Array.from(panel.field.gamma), Array.from(old.field.gamma));
  // The saved JSON field cannot retain signed zeros. Compare its exact JSON
  // representation; leave all live field values, including -0, untouched.
  for (const key of ['panels', 'basePanels']) assert.equal(JSON.stringify(panel.field[key].map(({ element, ...p }) => p)),
    JSON.stringify(old.field[key].map(({ element, ...p }) => p)));
  const points = [{ x: -2, y: -.4 }, { x: -1, y: .7 }, { x: 2, y: .5 }];
  for (const point of points) {
    assert.equal(streamfunctionAt(point, panel.field), streamfunctionAt(point, old.field));
    assert.deepEqual(velocityAt(point, panel.field), velocityAt(point, old.field));
  }
  assert.equal(potentialDifference(points[0], points[1], panel.field), potentialDifference(points[0], points[1], old.field));
});

test('reentrant topology cannot bypass potential-chart and trusted-field guards', () => {
  const { topology, panelSolution } = recovery();
  assert.throws(() => preparePanelStreamtubeTopology(input, {}, { crosslinePlacement: 'supplied-x' }), /nonmonotone body/);
  assert.throws(() => createPanelStreamtubeGrid(topology, { crosslinePlacement: 'supplied-x', panelSolution }), /shared panel potential chart/);
  assert.throws(() => createPanelStreamtubeGrid(topology, { crosslinePlacement: 'potential' }), /shared panel potential chart/);
  assert.equal(validatePanelTopologySeed(panelSolution, topology), panelSolution);
  assert.throws(() => validatePanelTopologySeed({ ...panelSolution }, topology), /differs/);
  assert.throws(() => validatePanelTopologySeed(panelSolution, { ...topology, alpha: topology.alpha + .1 }), /differs/);
  const changed = structuredClone(topology); changed.bodies[0].points[1].x += 1e-5;
  assert.throws(() => validatePanelTopologySeed(panelSolution, changed), /differs/);
  assert.throws(() => createInitialStreamtubeTopology(input, { parametricLowerBranches: 'yes' }), /controls/);
});
