import test from 'node:test';
import assert from 'node:assert/strict';
import { createCoupledShockBroadeningDriver } from '../scripts/validation/coupled-shock-broadening-driver.js';

const clone = value => structuredClone(value);
const checkpoint = () => ({ version: 1, families: { euler: 3, boundaryLayer: 0, edgeMatching: 0 },
  restart: { input: { upwind: { mcrit: .99, mucon: 1 } }, options: { transitionState: [0] },
    initialEuler: { x: [0, 4], nodes: [[1, 2]], undisplacedNodes: [[1, 2]] }, initialBL: [0] },
  continuation: { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing',
    blUpdate: 'xfoil', fractions: [[0, 1]], lastRedistributedStagnation: [0], preferredOrdering: 'amd', pivotTolerance: .001 } });
const envelope = cp => ({ version: 1, kind: 'mses-shock-broadening', checkpoint: cp,
  targetMcrit: .99, previousDensityChange: .2, acceptedUpdates: 10 });

function fixture({ residual = cp => Math.max(0, 3 - cp.restart.initialBL[0]), failRebase,
  reject = false, convex = true, mutateResult } = {}) {
  const seen = { rebases: [], solves: [] };
  const threshold = ({ targetMcrit, densityChange }) => targetMcrit - Math.min(.2, densityChange);
  const refresh = cp => {
    const r = residual(cp); cp.families = { euler: r, boundaryLayer: 0, edgeMatching: 0 }; return r;
  };
  const run = createCoupledShockBroadeningDriver({ threshold,
    densityChange: (a, b, n) => { assert.equal(n, 1); return Math.abs(Math.expm1(b[0] - a[0])); },
    rebase: (source, mcrit) => {
      seen.rebases.push(mcrit);
      if (failRebase?.(source, mcrit)) throw Object.assign(new Error('manufactured interface-pressure failure'), { code: 'DOMAIN' });
      const cp = clone(source); cp.restart.input.upwind.mcrit = mcrit; const r = refresh(cp);
      return { checkpoint: cp, residual: [r], diagnostics: { densityCount: 1, sourceMcrit: source.restart.input.upwind.mcrit,
        targetMcrit: mcrit, physicalStatePreserved: true } };
    },
    solve: (input, options) => {
      assert.equal(input, undefined); assert.equal(options.maxBacktracks, 12);
      assert.equal(options.iterationGeometry, options.resume.continuation.iterationGeometry);
      assert.equal(options.blUpdate, 'xfoil');
      seen.solves.push({ limit: options.maxIterations, mcrit: options.resume.restart.input.upwind.mcrit });
      const cp = clone(options.resume), history = [{ iteration: 0, residual: refresh(cp) }], linearDiagnostics = { solves: 0 };
      options.onCheckpoint(cp, { history, linearDiagnostics });
      if (options.maxIterations && history[0].residual > options.tolerance) {
        linearDiagnostics.solves = 1;
        if (!reject) {
          cp.restart.initialEuler.x[0] += Math.log1p(.1); cp.restart.initialBL[0]++;
          cp.restart.options.transitionState[0]++; cp.continuation.lastRedistributedStagnation[0] += .01;
          history.push({ iteration: 1, step: .25, residual: refresh(cp) });
          options.onCheckpoint(cp, { history, linearDiagnostics });
        }
      }
      const r = refresh(cp);
      const valid = typeof convex === 'function' ? convex(cp) : convex;
      const result = { checkpoint: cp, history, residual: [r], families: cp.families, linearDiagnostics,
        x: [...cp.restart.initialEuler.x, ...cp.restart.initialBL],
        flow: { nodes: clone(cp.restart.initialEuler.nodes), undisplacedNodes: clone(cp.restart.initialEuler.undisplacedNodes) },
        boundaryLayer: { transitionState: clone(cp.restart.options.transitionState) },
        mesh: { quality: { valid } }, converged: r <= options.tolerance && valid,
        reason: reject ? 'manufactured ordinary rejection' : !valid ? 'Invalid final displacement grid' : 'iteration limit',
        ...(reject ? { lastRejectedStep: { stage: 'admissibility', message: 'original rejection detail' } } : {}) };
      mutateResult?.(result, options); return result;
    } });
  return { run, seen };
}

test('accepted updates, applied density history and serialized resume exclude repeated iteration zero', () => {
  const source = checkpoint(), before = clone(source), a = fixture(), records = [], iterations = [];
  const full = a.run(source, { maxIterations: 3, onRecord: (kind, value) => {
    records.push({ kind, value: clone(value) });
    // Observers receive isolated objects, including the full CP.
    value.controller.checkpoint.restart.initialEuler.x[0] = 123;
  }, onIteration: entry => iterations.push(entry) });
  assert.equal(full.converged, true); assert.equal(full.targetRestored, true); assert.equal(full.targetChecked, true);
  assert.equal(full.operations.acceptedUpdates, 3); assert.equal(full.linearDiagnostics.solves, 3);
  assert.deepEqual(iterations.map(e => e.iteration), [1, 2, 3]);
  assert.equal(records.filter(r => r.kind === 'update').length, 3);
  assert.ok(Math.abs(full.controller.previousDensityChange - .1) < 1e-15);
  assert.ok(a.seen.solves.filter(s => s.limit === 1).slice(1).every(s => Math.abs(s.mcrit - .89) < 1e-15));
  const b = fixture(), first = b.run(source, { maxIterations: 1 });
  const continued = b.run(JSON.parse(JSON.stringify(first.controller)), { maxIterations: 2 });
  assert.equal(continued.converged, true); assert.deepEqual(continued.controller, full.controller);
  assert.equal(continued.linearDiagnostics.solves, 2); assert.equal(continued.linearDiagnostics.cumulativeSolves, 3);
  assert.deepEqual(source, before);
});

test('a broadened root restores target before Newton; a nonconvex target root cannot pass or loop', () => {
  const residual = cp => cp.restart.input.upwind.mcrit < .99 ? 0 : Math.max(0, 1 - cp.restart.initialBL[0]);
  const { run, seen } = fixture({ residual });
  const result = run(envelope(checkpoint()), { maxIterations: 1 });
  assert.equal(result.converged, true);
  assert.deepEqual(seen.solves.filter(s => s.limit === 1), [{ limit: 1, mcrit: .99 }]);
  assert.deepEqual(seen.rebases.slice(0, 2), [.79, .99]);
  const invalid = fixture({ residual: () => 0, convex: false });
  const failed = invalid.run(checkpoint(), { maxIterations: 0 });
  assert.equal(failed.converged, false); assert.equal(failed.targetChecked, true);
  assert.equal(failed.reason, 'Invalid final displacement grid');
  assert.deepEqual(invalid.seen.solves, [{ limit: 0, mcrit: .99 }]);
});

test('failed rebasing and failed Newton retain the previous complete phase/chart/density history', () => {
  const source = checkpoint(), before = clone(source);
  const a = fixture({ failRebase: (_, mcrit) => mcrit < .99 });
  const result = a.run(source, { maxIterations: 3 });
  assert.equal(result.converged, false); assert.equal(result.error.code, 'DOMAIN');
  assert.equal(result.error.stage, 'threshold rebase'); assert.equal(result.operations.acceptedUpdates, 1);
  assert.equal(result.controller.checkpoint.restart.initialBL[0], 1);
  assert.equal(result.controller.checkpoint.restart.options.transitionState[0], 1);
  assert.equal(result.controller.checkpoint.continuation.lastRedistributedStagnation[0], .01);
  assert.equal(result.controller.checkpoint.restart.input.upwind.mcrit, .99);
  assert.equal(result.linearDiagnostics.solves, 1); assert.deepEqual(source, before);
  const b = fixture({ reject: true }), wrapped = envelope(checkpoint()), rejected = b.run(wrapped, { maxIterations: 4 });
  assert.equal(rejected.operations.acceptedUpdates, 0); assert.equal(rejected.linearDiagnostics.solves, 1);
  assert.equal(rejected.controller.previousDensityChange, .2); assert.equal(rejected.controller.acceptedUpdates, 10);
  assert.deepEqual(rejected.controller.checkpoint.restart.initialEuler, wrapped.checkpoint.restart.initialEuler);
  assert.equal(rejected.error.lastRejectedStep.message, 'original rejection detail');
  assert.equal(b.seen.solves.length, 1);
});

test('observer cancellation retains an already committed complete update without another solve', () => {
  for (const cancellation of [null, undefined, 'stop requested', Object.freeze(new Error('stop requested'))]) {
    const { run, seen } = fixture(); let retained;
    let caught = false;
    try {
      run(checkpoint(), { maxIterations: 4, onRecord: (kind, value) => {
        if (kind === 'update') retained = value.controller;
      }, onIteration: () => { throw cancellation; } });
    } catch (error) { caught = true; assert.equal(error, cancellation); }
    assert.equal(caught, true); assert.equal(retained.acceptedUpdates, 1);
    assert.equal(retained.checkpoint.restart.initialBL[0], 1); assert.equal(retained.linearSolves, 1);
    assert.equal(seen.solves.length, 1);
  }
});

test('restoration failures retain target next action, and resumed tolerance/policy are explicit', () => {
  const residual = cp => cp.restart.input.upwind.mcrit < .99 ? 0 : 1;
  const failed = fixture({ residual, failRebase: (_, mcrit) => mcrit === .99 }).run(envelope(checkpoint()),
    { maxIterations: 1, tolerance: 1e-8 });
  assert.equal(failed.converged, false); assert.equal(failed.controller.nextAction, 'target-check');
  assert.equal(failed.controller.checkpoint.restart.input.upwind.mcrit, .79);
  const next = fixture({ residual }); next.run(failed.controller, { maxIterations: 0 });
  assert.equal(next.seen.rebases[0], .99);
  assert.equal(failed.controller.tolerance, 1e-8);
  assert.throws(() => next.run(failed.controller, { tolerance: 1e-10 }), /controls do not match/);
  assert.throws(() => next.run({ ...failed.controller, policy: {} }), /controls do not match/);
});

test('sampled intermediate reports remain resumable, and stale final physical output cannot certify a root', () => {
  const sampled = fixture({ convex: cp => cp.restart.initialBL[0] >= 3 });
  assert.equal(sampled.run(checkpoint(), { maxIterations: 3 }).converged, true);
  for (const alter of [r => { r.x[0] += 1; }, r => { r.flow.nodes[0][0] += 1; },
    r => { r.boundaryLayer.transitionState[0] += 1; }, r => { r.residual[0] = 1e-12; },
    r => { r.checkpoint.continuation.fractions[0][1] = .9; }]) {
    const { run } = fixture({ residual: () => 0, mutateResult: alter });
    const result = run(checkpoint(), { maxIterations: 0 });
    assert.equal(result.converged, false); assert.equal(result.targetChecked, false);
    assert.match(result.reason, /does not exactly match/);
  }
});
