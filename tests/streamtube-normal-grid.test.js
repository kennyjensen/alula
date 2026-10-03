import test from 'node:test';
import assert from 'node:assert/strict';
import { initializePanelStreamtubeBody, createPanelStreamtubeGrid } from '../src/euler/streamtube-body-initializer.js';
import { solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { directBodyConservation } from './oracles/streamtube-body.js';
import { stagnationPressureError } from './oracles/streamtube-entropy.js';
import { sparseProduct } from '../src/numerics/sparse.js';

const fixture = () => ({ ...intrinsicBodyFixture({ elements: 2, alpha: 2, bodySegments: 16, tubes: 5, tubeGrowth: 3,
  surfaceSpacing: 'cosine' }), streamwiseMode: 'isentropic' });
const controls = { crosslinePlacement: 'potential', outerCrosslineSpread: .15, normalSpacing: 'stagnation' };
const max = a => Math.max(...Array.from(a, Math.abs));

test('normal-spacing fit controls all four physical LE cells and the two-body Euler solution preserves mass, enthalpy and total pressure', () => {
  const input = fixture(), before = structuredClone(input), grid = initializePanelStreamtubeBody(input, controls), { system, initial } = grid;
  assert.deepEqual(input, before); assert.equal(streamtubeMeshSnapshot(grid).quality.valid, true);
  const targets = grid.diagnostics.normalAllocation.groups.flatMap(g => Object.values(g.targets));
  assert.equal(targets.length, 4);
  targets.forEach(t => assert.ok(Math.abs(t.actualAspect / 2.5 - 1) < 1.1e-4));
  const decoded = system.decode(initial);
  for (let g = 0; g < grid.input.weights.length; g++) {
    const expected = grid.input.captureLevels[g + 1] - grid.input.captureLevels[g];
    assert.ok(Math.abs(decoded.allocation.groups[g].reduce((sum, t) => sum + t.massFlow, 0) - expected) < 1e-14);
  }
  // Exercise normal-position, density and global capture/stagnation columns
  // together on the new distribution, independently of its analytic assembly.
  const direction = initial.map((_, i) => .1 * Math.sin(i + .4)), h = 1e-7;
  const matrix = system.jacobian(initial, { sparse: true }), product = sparseProduct(matrix, direction);
  const a = system.residual(initial.map((v, i) => v + h * direction[i])), b = system.residual(initial.map((v, i) => v - h * direction[i]));
  let error = 0;
  product.forEach((v, i) => { error = Math.max(error, Math.abs(v - (a[i] - b[i]) / (2 * h)) / Math.max(1, Math.abs(v))); });
  assert.ok(error < 2e-5, `Fitted-grid Jv error ${error}`);
  const r = solveStreamtubeBody(system, { initial, tolerance: 1e-10, maxIterations: 30 });
  assert.equal(r.converged, true, r.reason); assert.ok(r.linearDiagnostics.maxRelativeResidual <= 1e-10);
  const conservation = directBodyConservation(r, grid.input.bodies, system.conditions);
  for (const k of [0, 3]) assert.ok(Math.abs(conservation.balance[k]) < 2e-9);
  assert.ok(max(conservation.cutTraction) < 2e-9);
  for (const block of conservation.blocks) {
    assert.ok(max(block.internalCancellation) < 2e-9);
    for (const k of [0, 3]) assert.ok(block.maxLocal[k] < 2e-9);
  }
  const pressure = stagnationPressureError(r.sections.map(row => row.flat()), {
    gamma: system.conditions.gamma, referencePressure: system.conditions.pInf, freestreamMach: grid.input.mach });
  assert.ok(pressure.maxRelativeError < 2e-10);
  assert.equal(streamtubeMeshSnapshot({ system, nodes: r.nodes }).quality.valid, true);
  assert.match(r.forceStatus, /Unvalidated/);
});

test('fitted normal allocation preserves units and is insensitive to tighter tracing and smaller seeds', () => {
  const input = fixture(), a = createPanelStreamtubeGrid(input, controls);
  const b = createPanelStreamtubeGrid(input, { ...controls, relativeTolerance: 2e-10, seedDistance: 1e-4 });
  const scale = points => points.map(p => ({ x: 2 * p.x, y: 2 * p.y }));
  const doubled = { ...input, bodies: input.bodies.map(body => ({ ...body, points: scale(body.points), stagnationParameter: 2 * body.stagnationParameter })),
    outerLower: scale(input.outerLower), outerUpper: scale(input.outerUpper), cutPaths: input.cutPaths.map(scale) };
  const c = createPanelStreamtubeGrid(doubled, controls);
  assert.equal(a.system.layout.n, b.system.layout.n); assert.equal(a.system.layout.n, c.system.layout.n);
  let difference = 0;
  for (let g = 0; g < a.nodes.length; g++) for (let i = 0; i < a.nodes[g].length; i++) for (let j = 0; j < a.nodes[g][i].length; j++) {
    const p = a.nodes[g][i][j], q = b.nodes[g][i][j], doubled = c.nodes[g][i][j];
    difference = Math.max(difference, Math.hypot(p.x - q.x, p.y - q.y));
    assert.ok(Math.hypot(2 * p.x - doubled.x, 2 * p.y - doubled.y) < 1e-7);
  }
  assert.ok(difference < 2e-7, `Normal-spacing trace refinement moves the grid by ${difference}`);
  assert.throws(() => createPanelStreamtubeGrid(input, { ...controls, stagnationAspectRatio: 0 }), /controls/);
});

test('automatic normal matching preserves wider original tubes and the physical grid exactly', () => {
  const input = fixture(), before = structuredClone(input);
  const original = createPanelStreamtubeGrid(input, { ...controls, normalSpacing: 'supplied' });
  const automatic = createPanelStreamtubeGrid(input, { ...controls, normalSpacing: 'automatic', stagnationAspectRatio: .01 });
  assert.ok(automatic.diagnostics.normalAllocation.groups.every(g => g.retained));
  assert.deepEqual(automatic.input.weights, original.input.weights);
  assert.deepEqual(automatic.nodes, original.nodes);
  assert.deepEqual(input, before);
});
