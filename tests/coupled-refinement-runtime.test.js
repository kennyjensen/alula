// Browser preparation/zero-update replay against the retained paired seed.
// This performs no Newton updates, Jacobian assembly or linear factorization.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { coupledRefinementPlan, solveCoupledStreamtubeRefinement } from '../src/euler/streamtube-coupled-refinement-assembly.js';
import { checkpointDataEqual } from '../src/euler/streamtube-nested-checkpoint.js';
import { quadCoupledResultForDisplay } from '../src/ui/quad-coupled-result.js';
import { createQuadMeshProgress } from '../src/ui/quad-mesh-progress.js';

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const read = path => JSON.parse(fs.readFileSync(new URL(path, import.meta.url), 'utf8'));
const coarse = read('../docs/current-multielement-automatic-16x9-slor-browser.json');
const validated = read('../docs/current-multielement-paired-refined-seed.json');
const input = coarse.input, parent = { ...coarse.result, sourceCase: structuredClone(input) };

test('refinement plan enforces a converged unchanged case and the fixed node budget before flow evaluation', () => {
  const plan = coupledRefinementPlan(input, parent);
  assert.equal(plan.streamwiseFactor, 2);
  assert.equal(plan.normalInterpolation, 'streamfunction-quadratic');
  assert.deepEqual(plan.normalSubdivisions, validated.controls.normalSubdivisions);
  assert.deepEqual(plan.parent, { nx: 139, tubes: [12, 15, 12] });
  assert.deepEqual(plan.child, { nx: 278, tubes: [15, 21, 15] });
  assert.equal(plan.nodeCount, 15066);
  assert.equal(plan.maxNodes, 50000);
  assert.throws(() => coupledRefinementPlan(input, { ...parent, converged: false }), /converged coupled/);
  assert.throws(() => coupledRefinementPlan(input, { ...parent, mesh: { quality: { valid: false } } }), /valid grid/);
  assert.throws(() => coupledRefinementPlan(input, { ...parent, sourceCase: undefined }), /stale/);
  assert.throws(() => coupledRefinementPlan({ ...input, mach: .25 }, parent), /stale/);
  assert.throws(() => coupledRefinementPlan({ ...input, gridIntervals: 32 }, parent), /stale/);
  assert.throws(() => coupledRefinementPlan(input, { ...parent, checkpoint: undefined }), /complete continuation/);
  assert.throws(() => coupledRefinementPlan(input, { ...parent, families: { ...parent.families, boundaryLayer: 1 } }), /residuals differ/);
  const families = { ...parent.families, boundaryLayer: .001 };
  assert.throws(() => coupledRefinementPlan(input, { ...parent, families, checkpoint: { ...parent.checkpoint, families } }), /convergence tolerance/);
  assert.throws(() => coupledRefinementPlan(input, parent, { maxNodes: plan.nodeCount - 1 }), /node budget/);
  assert.throws(() => coupledRefinementPlan(input, parent, { maxNodes: 50001 }), /node budget/);
  const oversized = { ...parent, checkpoint: { ...parent.checkpoint,
    restart: { ...parent.checkpoint.restart, input: { ...parent.checkpoint.restart.input, outerLower: Array(1001) } } } };
  assert.throws(() => coupledRefinementPlan(input, oversized), /exceeding the 50,000 node budget/);
});

test('browser structural equality preserves strict numeric, typed-array and property distinctions', () => {
  assert.equal(checkpointDataEqual({ a: [1, 2], b: 3 }, { b: 3, a: [1, 2] }), true);
  assert.equal(checkpointDataEqual({ a: undefined }, {}), false);
  assert.equal(checkpointDataEqual([0], [-0]), false);
  assert.equal(checkpointDataEqual([], Array(1)), false);
  assert.equal(checkpointDataEqual(new Float64Array([1, 2]), new Float64Array([1, 2])), true);
  assert.equal(checkpointDataEqual(new Float64Array([1, 2]), [1, 2]), false);
  assert.equal(checkpointDataEqual({ a: [1, 2] }, { a: [1, 2.000000000000001] }), false);
});

test('runtime paired refinement exactly reproduces the validated seed, preserves history and publishes before zero-update resume', t => {
  assert.equal(coarse.passed, true); assert.equal(validated.passed, true);
  const before = serialize(parent), events = [], progress = createQuadMeshProgress(input.referenceChord);
  // A label from an earlier transfer must not survive this default linear one.
  // Only result metadata changes in this control; the accepted state is exact.
  const labeledParent = { ...parent, solverSettings: { ...parent.solverSettings, streamwiseInterpolation: 'surface-pchip' } };
  let latestMesh, firstMesh, prepared, latestCheckpoint;
  const raw = solveCoupledStreamtubeRefinement(input, labeledParent, { maxIterations: 0,
    onPrepared: value => { prepared = value; },
    onCheckpoint: value => { latestCheckpoint = value; },
    onStage: stage => { events.push({ kind: 'stage', stage: stage.stage }); progress.stage(stage); },
    onMesh: (mesh, stage, state) => {
      events.push({ kind: 'mesh', stage }); latestMesh = progress.mesh(mesh); firstMesh ??= latestMesh;
      if (state) {
        assert.ok(prepared.system.bl.stations.length > 0);
        assert.deepEqual(serialize(latestCheckpoint), validated.checkpoint);
        assert.deepEqual(state.flow.nodes, latestCheckpoint.restart.initialEuler.nodes);
      }
    },
    onIteration: iteration => events.push({ kind: 'iteration', iteration: iteration.iteration }),
  });
  assert.deepEqual(serialize(raw.checkpoint), validated.checkpoint);
  assert.equal(Object.hasOwn(raw.solverSettings, 'streamwiseInterpolation'), false);
  assert.equal(labeledParent.solverSettings.streamwiseInterpolation, 'surface-pchip');
  assert.deepEqual(serialize(raw.refinement.transfer), validated.transfer);
  assert.deepEqual(serialize(raw.restart), validated.restart);
  assert.deepEqual(raw.families, validated.result.families);
  assert.deepEqual(raw.refinement.parent, { nx: 139, tubes: [12, 15, 12], unknowns: 12184, cells: 5421 });
  assert.deepEqual(raw.refinement.child, { nx: 278, tubes: [15, 21, 15], unknowns: 30999, cells: 14178 });
  assert.equal(raw.refinement.level, 1);
  assert.equal(raw.refinement.transfer.initialSMOVERepeated, false);
  assert.equal(raw.refinement.transfer.maintenanceHistoryPreserved, true);
  assert.equal(raw.initialRedistribution.resumed, true);
  assert.deepEqual(raw.initialRedistribution.passages, []);
  assert.equal(raw.linearDiagnostics.solves, 0);
  assert.equal(raw.history.length, 1);
  assert.equal(raw.converged, false);
  assert.equal(raw.boundaryLayer.surfaces.length, 4);
  assert.equal(raw.boundaryLayer.wakes.length, 2);
  assert.deepEqual(events.filter(event => event.kind === 'stage').map(event => event.stage), ['coupled-refinement', 'coupled']);
  assert.ok(events.findIndex(event => event.kind === 'mesh') < events.findIndex(event => event.kind === 'iteration'));
  assert.equal(firstMesh.cells.length, 14178);
  assert.equal(firstMesh.quality.valid, true);
  assert.equal(firstMesh.flow.iteration, 0);
  assert.deepEqual(firstMesh.vertices, latestMesh.vertices);
  const display = quadCoupledResultForDisplay(raw, input, latestMesh);
  assert.ok(['cl', 'cd', 'cm'].every(key => Number.isFinite(display[key])));
  assert.equal(display.coefficientStatus, 'unconverged');
  assert.equal(display.physicalAcceptance, false);
  assert.deepEqual(display.sourceCase, input);
  assert.equal(display.initialization.attempts[0].kind, 'paired-nested-refinement');
  assert.deepEqual(serialize(parent), before);
  t.diagnostic(JSON.stringify({ unknowns: raw.x.length, cells: display.mesh.cells.length,
    exactValidatedSeed: true, exactValidatedTransfer: true, initialSMOVERepeated: false,
    linearSolves: raw.linearDiagnostics.solves, newtonUpdates: raw.history.length - 1,
    meshMessages: events.filter(event => event.kind === 'mesh').length }));
});
