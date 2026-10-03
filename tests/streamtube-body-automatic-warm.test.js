// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { solveStreamtubeBodyAutomatic } from '../src/euler/tests/streamtube-body-automatic.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const tolerance = 1e-10;
const budgets = { tolerance, directMaxIterations: 0, stageMaxIterations: 15,
  initialFractionStep: .5, maxStages: 6 };
const upwind = { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } };
const sources = new Map();
function source(mode = 'momentum') {
  if (!sources.has(mode)) {
    const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2, mach: .2 }),
      flowModel: 'compressible', streamwiseMode: mode, upwind: structuredClone(upwind),
      ...(mode === 'hybrid' ? { hybrid: { epsilonP: .001 } } : {}) };
    const root = solveStreamtubeBodyAutomatic(input, { tolerance, directMaxIterations: 15 });
    assert.equal(root.converged, true, root.reason);
    assert.equal(root.system.layout.n, 113);
    assert.equal(streamtubeMeshSnapshot({ system: root.system, nodes: root.nodes }).quality.valid, true);
    sources.set(mode, { root, seed: { input, initialEuler: { x: Array.from(root.x), nodes: structuredClone(root.nodes) } } });
  }
  const saved = sources.get(mode);
  return { root: saved.root, initialFlow: structuredClone(saved.seed), input: structuredClone(saved.seed.input) };
}
function assertPhysicalIdentity(a, b) {
  assert.deepEqual(a.nodes, b.nodes);
  assert.deepEqual(a.captured, b.captured);
  assert.deepEqual(a.stagnation, b.stagnation);
  assert.deepEqual(a.strengths, b.strengths);
  assert.deepEqual(a.allocation, b.allocation);
  for (let i = 0; i < a.sections.length; i++) for (let g = 0; g < a.sections[i].length; g++)
    for (let j = 0; j < a.sections[i][g].length; j++) {
      assert.equal(a.sections[i][g][j].rho, b.sections[i][g][j].rho);
      assert.equal(a.sections[i][g][j].q, b.sections[i][g][j].q);
    }
}

test('same-Mach automatic full-flow entry preserves the accepted physical state with zero Newton updates', () => {
  const { root, initialFlow, input } = source(), before = structuredClone({ input, initialFlow });
  const result = solveStreamtubeBodyAutomatic(input, { ...budgets, stageMaxIterations: 0, initialFlow });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.continuation.reachedTarget, true);
  assert.equal(result.continuation.currentMach, input.mach);
  assert.equal(result.continuation.attempts.length, 1);
  assert.equal(result.continuation.attempts[0].label, 'target-warm');
  assert.equal(result.continuation.attempts[0].iterations, 0);
  assert.equal(result.linearDiagnostics.solves, 0);
  assert.equal(result.continuation.initialFlow.coldInitializationSkipped, true);
  assert.equal(result.continuation.initialFlow.completePhysicalState, true);
  assert.equal(result.continuation.densityReinitializedOnWarmRestart, false);
  assert.match(result.continuation.attempts[0].startup.densityInitialization, /no isentropic inversion/);
  assertPhysicalIdentity(root, result);
  assert.deepEqual(result.residual, root.residual);
  assert.deepEqual({ input, initialFlow }, before);
  assert.equal(result.physicalAcceptance, false);
});

test('failed warm increments halve from the last accepted root and cannot label its small residual as target convergence', () => {
  const { root, initialFlow, input } = source(); input.mach = .3;
  const result = solveStreamtubeBodyAutomatic(input, { ...budgets, initialFlow,
    stageMaxIterations: 0, maxSubdivisions: 2, maxStages: 8 });
  assert.equal(result.converged, false);
  assert.equal(result.continuation.reachedTarget, false);
  assert.equal(result.continuation.currentMach, .2);
  assert.equal(result.continuation.targetMach, .3);
  assert.equal(result.system.conditions.mach, .2);
  assert.ok(result.diagnostics.residual <= tolerance);
  assert.match(result.reason, /stopped at Mach 0.2/);
  const attempts = result.continuation.attempts;
  assert.deepEqual(attempts.map(a => a.fraction), [1, .5, .25, .125]);
  assert.ok(attempts.every(a => !a.converged && a.iterations === 0 && a.gridValid));
  assert.ok(attempts.every(a => a.startup.sourceMach === .2 && a.startup.physicalDensityPreserved));
  assertPhysicalIdentity(root, result);
});

test('warm fallback reaches the exact target with copied caller data and detached callbacks', () => {
  const { initialFlow, input } = source('hybrid'); input.mach = .3;
  const reference = structuredClone({ input, initialFlow }), stages = []; let meshes = 0;
  const result = solveStreamtubeBodyAutomatic(input, { ...budgets, initialFlow,
    onStage: event => {
      stages.push({ ...event }); event.mach = NaN; event.fraction = -1;
      input.mach = NaN; input.hybrid.epsilonP = NaN; input.upwind.mcrit = NaN;
      initialFlow.input.mach = NaN; initialFlow.input.hybrid.epsilonP = NaN;
      initialFlow.initialEuler.x.fill(NaN);
      initialFlow.initialEuler.nodes[0][0][0].x = NaN;
    },
    onIteration: event => { event.residual = NaN; event.mach = NaN; },
    onMesh: event => {
      meshes++;
      for (const key of ['system', 'flow', 'nodes']) assert.equal(Object.hasOwn(event, key), false);
      assert.equal(event.mesh.quality.valid, true);
      event.mesh.vertices.forEach(p => { p.x = NaN; p.y = NaN; });
      event.mesh.cells.forEach(c => c.fill(-1));
      event.mesh.quality.valid = false; event.iteration.residual = NaN;
    },
  });
  assert.equal(result.converged, true, result.reason);
  assert.equal(result.continuation.reachedTarget, true);
  assert.equal(result.continuation.currentMach, reference.input.mach);
  assert.equal(result.system.conditions.mach, reference.input.mach);
  assert.equal(result.system.conditions.streamwiseMode, 'hybrid');
  assert.deepEqual(result.system.conditions.hybrid, reference.input.hybrid);
  assert.deepEqual(result.system.conditions.upwind, reference.input.upwind);
  assert.ok(meshes > 0);
  assert.ok(result.diagnostics.residual <= tolerance);
  const attempts = result.continuation.attempts;
  assert.deepEqual(attempts.map(a => a.mach), [.3, .25, .3]);
  assert.equal(attempts[0].converged, false);
  assert.ok(attempts.slice(1).every(a => a.converged && a.gridValid));
  assert.ok(attempts.every(a => a.startup.physicalMassPreserved && a.startup.physicalDensityPreserved
    && a.history.every(h => Number.isFinite(h.residual))));
  assert.deepEqual(attempts.map(({ label, mach, fraction }) => ({ label, mach, fraction, targetMach: .3 })), stages);
  assert.equal(streamtubeMeshSnapshot({ system: result.system, nodes: result.nodes }).quality.valid, true);
});

test('warm stage limit returns the last accepted intermediate operating point explicitly unconverged', () => {
  const { initialFlow, input } = source(); input.mach = .3;
  const result = solveStreamtubeBodyAutomatic(input, { ...budgets, initialFlow, maxStages: 2 });
  assert.equal(result.converged, false);
  assert.equal(result.continuation.reachedTarget, false);
  assert.equal(result.continuation.currentMach, .25);
  assert.equal(result.continuation.targetMach, .3);
  assert.ok(result.diagnostics.residual <= tolerance);
  assert.ok(result.continuation.attempts.at(-1).converged);
  assert.match(result.reason, /stage limit/);
});

test('complete-flow entry rejects invalid sources and changed equations before attempting continuation', () => {
  const { initialFlow, input } = source('hybrid'); let stages = 0;
  const run = (caseInput = input, seed = initialFlow, extra = {}) => solveStreamtubeBodyAutomatic(caseInput,
    { ...budgets, initialFlow: seed, onStage: () => { stages++; }, ...extra });
  assert.throws(() => run(input, initialFlow, { initialGeometry: {} }), /mutually exclusive/);
  assert.throws(() => run(input, initialFlow, { initialEuler: {} }), /mutually exclusive/);
  for (const invalid of [null, {}, { input, initialEuler: { x: [], nodes: initialFlow.initialEuler.nodes } },
    { ...initialFlow, initialEuler: { ...initialFlow.initialEuler, x: [NaN] } }])
    assert.throws(() => run(input, invalid), /initialFlow|finite encoded source state/);
  const changed = structuredClone(initialFlow); changed.initialEuler.x[0] += .001;
  assert.throws(() => run(input, changed), /source must already satisfy the residual tolerance/);
  assert.throws(() => run({ ...input, hybrid: { epsilonP: .002 } }), /changed hybrid/);
  assert.throws(() => run({ ...input, upwind: { ...upwind, mcrit: .98 } }), /changed upwind/);
  const momentum = { ...input, streamwiseMode: 'momentum' }; delete momentum.hybrid;
  assert.throws(() => run(momentum), /changed streamwiseMode/);
  const folded = structuredClone(initialFlow);
  // An interior free streamline node is independent of the solid contour.
  // Move it beyond the next streamwise station: positive residual tolerance
  // is immaterial because a folded source is rejected before its evaluation.
  const nodes = folded.initialEuler.nodes;
  nodes[0][2][1].x = nodes[0][4][1].x;
  const checkSystem = createStreamtubeBodySystem(folded.input);
  assert.equal(streamtubeMeshSnapshot({ system: checkSystem, nodes }).quality.valid, false);
  assert.throws(() => run(input, folded), /valid convex quadrilateral grid/);
  assert.equal(stages, 0);
});

test('warm continuation propagates cancellation without trying another operating point', () => {
  for (const callback of ['onStage', 'onIteration', 'onMesh']) {
    const { initialFlow, input } = source(); input.mach = .3;
    const sentinel = new Error(`cancel warm ${callback}`); let stages = 0, calls = 0;
    const options = { ...budgets, directMaxIterations: 15, initialFlow, onStage: () => { stages++; } };
    options[callback] = () => { calls++; throw sentinel; };
    assert.throws(() => solveStreamtubeBodyAutomatic(input, options), e => e === sentinel);
    assert.equal(calls, 1);
    assert.equal(stages, callback === 'onStage' ? 0 : 1);
  }
});

test('omitting initialFlow retains the archived cold wrapper results exactly', async () => {
  const metadata = JSON.parse(readFileSync(new URL('../docs/hybrid-flow/body-automatic-warm-before.json', import.meta.url)));
  const bytes = readFileSync(new URL(`../${metadata.archive}`, import.meta.url));
  assert.equal(createHash('sha256').update(bytes).digest('hex'), metadata.sha256);
  const archived = bytes.toString().replace(/from '(\.\/[^']+)'/g, (_, path) => {
    const original = new URL(path, new URL('../src/euler/streamtube-body-automatic.js', import.meta.url));
    const resolved = existsSync(original) ? original
      : new URL(path, new URL('../src/euler/tests/streamtube-body-automatic.js', import.meta.url));
    return `from '${resolved.href}'`;
  });
  const previous = (await import(`data:text/javascript;base64,${Buffer.from(archived).toString('base64')}`)).solveStreamtubeBodyAutomatic;
  const input = intrinsicBodyFixture({ bodySegments: 4, tubes: 2, mach: .2 });
  const options = { ...budgets, stageMaxIterations: 0, maxStages: 2 };
  const old = previous(input, options), current = solveStreamtubeBodyAutomatic(input, options);
  const { system: oldSystem, ...oldData } = old, { system: currentSystem, ...currentData } = current;
  assert.deepEqual(currentData, oldData);
  assert.deepEqual(currentSystem.conditions, oldSystem.conditions);
  assert.equal(Object.hasOwn(current.continuation, 'initialFlow'), false);
});
