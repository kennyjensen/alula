// SPDX-License-Identifier: GPL-2.0-or-later
// A domain-transfer regression on a small coupled root, not a physical
// airfoil-accuracy claim. Production equations and Newton remain unchanged.
import test from 'node:test';
import assert from 'node:assert/strict';
import { extendCoupledLaterally } from '../scripts/validation/coupled-lateral-extension.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { adjustStreamtubeInlets } from '../src/geometry/streamtube-grid-maintenance.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directStreamtubeVolumeGeometry } from './oracles/streamtube-control-volume-geometry.js';

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const create = ({ restart: f }) => createCoupledStreamtubeBody(f.input, {
  ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL,
});
const maximum = values => Math.max(0, ...Array.from(values, Math.abs));
const near = (a, b, tolerance = 2e-12) => {
  if (typeof a === 'number') {
    assert.ok(Number.isFinite(b) && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
  } else if (a && typeof a === 'object') {
    assert.deepEqual(Object.keys(a), Object.keys(b));
    for (const key of Object.keys(a)) near(a[key], b[key], tolerance);
  } else assert.deepEqual(a, b);
};
const convex = nodes => {
  const geometry = directStreamtubeVolumeGeometry(nodes);
  assert.equal(geometry.valid, true);
  assert.deepEqual(geometry.concavePrimal, []);
};

test('lateral extension preserves a nonzero capture root and reconverges the complete two-element system', t => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'isentropic', wakeGeometry: 'independent-banks', wakeOutlet: 'banks' };
  const inputBefore = structuredClone(input);
  const controls = { reynolds: 1e5, edgeMatching: 'section-velocity', maxIterations: 12,
    tolerance: 1e-10, stepAcceptance: 'admissible' };
  const root = solveCoupledStreamtubeIses(input, controls);
  assert.equal(root.converged, true, root.reason);
  assert.deepEqual(input, inputBefore);
  const checkpoint = serialize(root.checkpoint), frozen = structuredClone(checkpoint);
  const parent = create(checkpoint), before = parent.evaluate(parent.initial), old = parent.euler.layout;
  assert.deepEqual(before.families, checkpoint.families);
  const capture = old.globals.capture.find(col => col !== null);
  assert.ok(Math.abs(parent.initial[capture]) > 1e-6, 'Exercise a material nonzero normalized capture coordinate.');

  const lowerWidths = [.2, .3], upperWidths = [.4];
  const result = extendCoupledLaterally(checkpoint, { lowerWidths, upperWidths });
  const child = result.system, after = child.evaluate(child.initial), layout = child.euler.layout;
  const offset = g => g === 0 ? 2 : 0;
  assert.deepEqual(layout.tubes, [4, 2, 3]);
  assert.equal(layout.nx, old.nx);
  assert.deepEqual(result.checkpoint.continuation, checkpoint.continuation);
  assert.deepEqual(child.bl.stations, parent.bl.stations);
  assert.deepEqual(child.bl.snapshotActive(), parent.bl.snapshotActive());
  assert.deepEqual(child.bl.trips, parent.bl.trips);
  assert.deepEqual(child.euler.fractions, parent.euler.fractions);
  assert.deepEqual(child.conditions, parent.conditions);
  for (const [key, value] of Object.entries(parent.euler.conditions))
    if (key !== 'massScale') assert.deepEqual(child.euler.conditions[key], value);
  near(after.outer.allocation.totalMass, before.outer.allocation.totalMass + .9);
  near(after.outer.captured, before.outer.captured.map((v, i) => i === 0 ? v - .5 : i === old.elements + 1 ? v + .4 : v));
  for (const col of layout.globals.capture) if (col !== null) assert.equal(child.initial[col], 0);
  assert.deepEqual(after.outer.stagnation, before.outer.stagnation);
  assert.deepEqual(after.outer.strengths, before.outer.strengths);
  near(after.outer.allocation.groups[0].slice(0, 2).map(tube => tube.massFlow), [.3, .2]);
  near(after.outer.allocation.groups.at(-1).at(-1).massFlow, .4);

  // Audit actual cells and both node representations, independently of the
  // transfer helper's diagnostics and normalized capture implementation.
  for (const field of ['nodes', 'undisplacedNodes']) before.outer[field].forEach((grid, g) => grid.forEach((row, i) =>
    row.forEach((p, j) => near(after.outer[field][g][i][j + offset(g)], p))));
  before.outer.allocation.groups.forEach((group, g) => group.forEach((tube, j) => {
    const jj = j + offset(g);
    near(after.outer.allocation.groups[g][jj].massFlow, tube.massFlow);
    for (let i = 0; i < old.nx; i++) {
      assert.equal(child.initial[layout.densityIndex(i, g, jj)], parent.initial[old.densityIndex(i, g, j)]);
      near(after.outer.sections[i][g][jj], before.outer.sections[i][g][j]);
    }
    for (let i = 0; i < old.nx - 1; i++) near(after.outer.cells[i][g][jj], before.outer.cells[i][g][j]);
  }));
  assert.deepEqual(child.initial.slice(child.ne), parent.initial.slice(parent.ne));
  near(after.layers.states, before.layers.states, 1e-14);
  near(after.layers.residual, before.layers.residual);
  near(after.residual.slice(child.ne), before.residual.slice(parent.ne));
  near(after.edges, before.edges);
  near(after.outer.outletBankTangency, before.outer.outletBankTangency);
  for (const body of old.bodies.keys()) for (const side of ['upper', 'lower'])
    for (let i = 1; i < old.nx; i++) near(after.outer.bodyPressure(body, i, side), before.outer.bodyPressure(body, i, side));

  const parentInlet = adjustStreamtubeInlets(before.outer.nodes, old.bodies, checkpoint.continuation.fractions);
  const childInlet = adjustStreamtubeInlets(after.outer.nodes, layout.bodies, result.checkpoint.continuation.fractions);
  for (let b = 0; b < old.elements; b++) for (let i = 0; i <= old.bodies[b].leadingIndex; i++) {
    near(childInlet.nodes[b][i].at(-1), parentInlet.nodes[b][i].at(-1));
    near(childInlet.nodes[b + 1][i][0], parentInlet.nodes[b + 1][i][0]);
  }
  convex(after.outer.nodes);

  const serialized = serialize(result.checkpoint), serializedBefore = structuredClone(serialized), replay = create(serialized);
  assert.deepEqual(replay.evaluate(replay.initial).residual, after.residual);
  const zero = solveCoupledStreamtubeIses(undefined, { ...controls, resume: serialized, maxIterations: 0 });
  assert.equal(zero.initialRedistribution.resumed, true);
  assert.deepEqual(zero.initialRedistribution.passages, []);
  assert.deepEqual(zero.families, after.families);
  assert.deepEqual(zero.flow.nodes, after.outer.nodes);
  assert.deepEqual(serialize(zero.checkpoint.continuation), checkpoint.continuation);

  const solved = solveCoupledStreamtubeIses(undefined, { ...controls, resume: serialized });
  assert.equal(solved.converged, true, solved.reason);
  assert.equal(solved.boundaryLayer.surfaces.length, 4);
  assert.equal(solved.boundaryLayer.wakes.length, 2);
  assert.ok(maximum(Object.values(solved.families)) < 1e-10);
  assert.ok(solved.linearDiagnostics.maxRelativeResidual < 1e-10);
  convex(solved.flow.nodes);
  const final = create(serialize(solved.checkpoint)), x = final.initial, jacobian = final.jacobian(x), h = 2e-6;
  const errors = [];
  for (const family of ['euler', 'bl', 'both']) {
    const d = x.map((q, i) => (family === 'euler' && i >= final.ne || family === 'bl' && i < final.ne) ? 0
      : (i < final.ne ? .001 : Math.max(.01, Math.abs(q))) * Math.sin(i * .43 + .2));
    const exact = sparseProduct(jacobian, d);
    const plus = final.residual(x.map((q, i) => q + h * d[i])), minus = final.residual(x.map((q, i) => q - h * d[i]));
    const error = maximum(exact.map((q, i) => {
      const fd = (plus[i] - minus[i]) / (2 * h);
      return (fd - q) / Math.max(1, Math.abs(fd), Math.abs(q));
    }));
    assert.ok(error < 5e-6, `${family}: ${error}`); errors.push({ family, error });
  }
  final.evaluate(x);
  assert.deepEqual(serialized, serializedBefore);
  assert.deepEqual(checkpoint, frozen);
  assert.deepEqual(serialize(root.checkpoint), frozen);
  assert.deepEqual(parent.evaluate(parent.initial).residual, before.residual);
  t.diagnostic(JSON.stringify({ parentUnknowns: parent.n, unknowns: child.n,
    parentCaptureCoordinate: parent.initial[capture], parentIterations: root.history.length - 1,
    iterations: solved.history.length - 1, families: solved.families, errors, transfer: result.diagnostics }));
});
