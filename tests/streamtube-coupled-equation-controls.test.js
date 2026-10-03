import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { streamtubeEquationControls } from '../src/euler/streamtube-equation-selection.js';
import { coupledAssemblyConditions, solveCoupledStreamtubeAssembly } from '../src/euler/streamtube-coupled-assembly.js';
import { coupledMachPlan } from '../src/euler/tests/streamtube-coupled-mach-assembly.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';

test('explicit equations retain Reynolds/trip normalization and select compressible BL gas before grid work', () => {
  const input = { elements: [{}, {}], referenceChord: 2, reynolds: 3e6,
    materialTrips: [[.2, .3], [.4, .5]], transitionMode: 'automatic' };
  const bodies = [{ element: 1 }, { element: 0 }];
  const baseline = coupledAssemblyConditions(input, bodies, 1);
  for (const eulerIsmom of [1, 2, 3, 4]) {
    const next = coupledAssemblyConditions({ ...input, eulerIsmom }, bodies, 1);
    assert.equal(next.options.blThermodynamics, 'historical-common-isentrope');
    delete next.options.blThermodynamics;
    assert.deepEqual(next, baseline);
  }
  for (const eulerIsmom of [0, 5, '2', 2.5, null, NaN]) {
    assert.throws(() => solveCoupledStreamtubeAssembly({ ...input, eulerIsmom }), /ISMOM/);
    assert.throws(() => coupledMachPlan({ ...input, flowModel: 'streamtube-grid', quadBoundaryLayers: true, eulerIsmom }), /ISMOM/);
  }
});

// Structural cache/routing fixture only. Planning must reject changing the
// operator before reconstructing a flow; this fixture is not a numerical root.
function planningParent(ismom) {
  const sourceCase = { flowModel: 'streamtube-grid', quadBoundaryLayers: true,
    eulerIsmom: ismom, mach: .2, referenceChord: 1, reynolds: 1e6, elements: [{}] };
  const bodies = [{ element: 0 }], settings = coupledAssemblyConditions(sourceCase, bodies, 1);
  const equations = streamtubeEquationControls(ismom), nodes = [[[{ x: 0, y: 0 }]]];
  const families = { euler: 0, boundaryLayer: 0, edgeMatching: 0 };
  return { model: 'research-streamtube-euler-bl', converged: true, mesh: { quality: { valid: true } },
    sourceCase, ...settings.normalization, families,
    solverSettings: { tolerance: 1e-10, ...equations },
    checkpoint: { version: 1, families, continuation: {}, restart: {
      input: { bodies, mach: .2, ...equations }, options: settings.options,
      initialEuler: { x: [0], nodes, undisplacedNodes: structuredClone(nodes) }, initialBL: [.03, 1, 2, 1] } } };
}

test('Mach routing keeps explicit ISMOM and rejects a cached operator from another selection', () => {
  for (const ismom of [1, 2, 3, 4]) {
    const parent = planningParent(ismom), target = { ...parent.sourceCase, mach: .74 };
    const before = structuredClone(parent);
    const plan = coupledMachPlan(target, parent);
    assert.equal(plan.route, 'warm-hybrid');
    assert.equal(plan.targetMach, .74);
    assert.deepEqual(parent, before);
    assert.equal(coupledMachPlan(target).route, 'cold-baseline');
    assert.throws(() => coupledMachPlan(target, undefined, { epsilonP: 2e-5 }), /public dissipation controls/);
    assert.throws(() => coupledMachPlan(target, undefined, {
      upwind: { mucon: 1, mcrit: .95, boundary: { kind: 'unfiltered-first-two' } },
    }), /public dissipation controls/);
    assert.throws(() => coupledMachPlan({ ...target, eulerIsmom: ismom % 4 + 1 }, parent), /stale/);
    const mismatched = structuredClone(parent);
    mismatched.checkpoint.restart.input.hybrid.ismom = ismom % 4 + 1;
    assert.throws(() => coupledMachPlan(target, mismatched), /parent must use/);
    delete mismatched.checkpoint.restart.input.hybrid.ismom;
    assert.throws(() => coupledMachPlan(target, mismatched), /parent must use/);
  }
});

test('all four coupled operators accept the same supplied physical BL state without a profile march', () => {
  const saved = JSON.parse(fs.readFileSync(new URL('fixtures/streamtube-isentropic-coupled-root.json', import.meta.url)));
  const before = JSON.stringify(saved);
  let baseline;
  for (const ismom of [1, 2, 3, 4]) {
    const input = { ...saved.input, ...streamtubeEquationControls(ismom) };
    const system = createCoupledStreamtubeBody(input, { ...saved.options, edgeMatching: 'section-velocity',
      blThermodynamics: 'historical-common-isentrope', initialEuler: saved.initialEuler, initialBL: saved.initialBL });
    assert.equal(system.initialization.suppliedBL, true);
    assert.equal(system.euler.conditions.hybrid.ismom, ismom);
    assert.deepEqual(Array.from(system.initial.subarray(system.ne)), saved.initialBL);
    const value = system.evaluate(system.initial);
    assert.ok(value.residual.every(Number.isFinite));
    const physical = { nodes: value.outer.nodes, layers: value.layers.states,
      residual: Array.from(value.residual.subarray(system.ne)) };
    if (baseline) assert.deepEqual(physical, baseline);
    else baseline = physical;
  }
  assert.equal(JSON.stringify(saved), before);
});
