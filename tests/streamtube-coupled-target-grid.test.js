import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { targetMachGridPlan, solveCoupledTargetGrid } from '../src/euler/tests/streamtube-coupled-target-grid.js';
const input = { flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .74, alpha: 2.68,
  gridIntervals: 64, gridTubes: 7, transitionMode: 'automatic', materialTrips: [[1, 1]], elements: [{}] };

test('target-Mach grid ordering preserves low grids, explicit warm/prepared cases and nonterminal trips', () => {
  assert.equal(targetMachGridPlan(input).coarseIntervals, 32);
  assert.equal(targetMachGridPlan({ ...input, gridIntervals: 32 }).coarseIntervals, 16);
  for (const patch of [{ gridIntervals: 16 }, { mach: .2 }, { transitionMode: 'fixed-trip' }, { materialTrips: [[.5, 1]] }])
    assert.equal(targetMachGridPlan({ ...input, ...patch }), null);
  for (const options of [{ enabled: false }, { parentResult: {} }, { preparedEuler: {} }])
    assert.equal(targetMachGridPlan(input, options), null);
  assert.throws(() => targetMachGridPlan(input, { enabled: null }));
});

test('coarse failure stays provisional and coarse accepted events cannot be cached as fine roots', () => {
  const before = structuredClone(input), events = [];
  const result = solveCoupledTargetGrid(input, { plan: targetMachGridPlan(input),
    onCheckpoint: (_cp, details) => events.push(details),
    solveCoarse: (coarse, observers) => {
      assert.equal(coarse.gridIntervals, 32); assert.equal(coarse.mach, .74);
      observers.onCheckpoint({}, { kind: 'accepted', reachedTarget: true, mach: .72 });
      return { converged: false, stateConverged: true, mach: .72, reason: 'stage limit', sourceCase: coarse, mesh: { quality: { valid: true } } };
    } });
  assert.equal(result.converged, false); assert.equal(result.gridSequence.reachedTarget, false);
  assert.equal(result.stateConverged, true, 'retain the lower-Mach root without calling the target solved');
  assert.equal(result.sourceCase.gridIntervals, 32); assert.equal(result.requestedCase.gridIntervals, 64);
  assert.equal(events[0].kind, 'coarse-accepted'); assert.equal(events[0].reachedTarget, false);
  assert.deepEqual(input, before);
});

// Explicit adapter stubs exercise ordering/metadata; they do not claim a PDE root.
test('successful target-Mach transfer retains tube counts and requires the fine root', async () => {
  const order = [], normalization = { referenceChord: 1, solverLength: 1, referenceReynolds: 1e6, kernelReynolds: 1e6 };
  const opts = { reynolds: 1e6, ncrit: 4, edgeMatching: 'section-velocity', tripFractions: [[1, 1]] };
  const checkpoint = { restart: { input: { mach: .74, alpha: 2.68, weights: [Array(10).fill(1), Array(10).fill(1)] }, options: opts } };
  const root = { converged: true, mach: .74, mesh: { quality: { valid: true }, cells: [] }, checkpoint,
    sourceCase: { ...input, gridIntervals: 32 }, ...normalization, initialization: {}, solverSettings: {} };
  let fineConverged = true;
  globalThis.__targetGridTest = {
    checkpointDataEqual: (a, b) => JSON.stringify(a) === JSON.stringify(b),
    coupledAssemblyConditions: () => ({ options: opts, normalization }),
    createCoupledStreamtubeBody: () => ({ bl: { hasFiniteBase: false } }),
    prepareCoupledGridSequence: args => {
      order.push('refine'); assert.equal(args.sourceResult.converged, true);
      assert.equal(args.input.mach, .74); assert.deepEqual(args.input.weights.map(r => r.length), [10, 10]);
      assert.equal(args.sourceGridIntervals, 32); assert.equal(args.requestedGridIntervals, 64);
      return { system: {}, transfer: { nominalGridIntervals: 64 } };
    },
    solveCoupledGridLevel: () => { order.push('fine-solve'); return { ...root, converged: fineConverged,
      reason: fineConverged ? 'residual' : 'iteration limit', solverInput: { bodies: [{}] }, linearDiagnostics: { solves: 3 } }; },
  };
  try {
    const code = fs.readFileSync('src/euler/tests/streamtube-coupled-target-grid.js', 'utf8')
      .replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis.__targetGridTest;`);
    const { solveCoupledTargetGrid: solve } = await import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
    const options = { plan: targetMachGridPlan(input), solveCoarse: () => { order.push('coarse-target'); return root; } };
    const result = solve(input, options);
    assert.deepEqual(order, ['coarse-target', 'refine', 'fine-solve']);
    assert.equal(result.converged, true); assert.equal(result.gridSequence.reachedTarget, true);
    let returned = false, caught = false;
    try { solve(input, { ...options, onStage: () => { throw undefined; } }); returned = true; }
    catch (error) { caught = true; assert.equal(error, undefined); }
    assert.equal(returned, false); assert.equal(caught, true, 'observer cancellation must escape unchanged');
    fineConverged = false;
    const failed = solve(input, options);
    assert.equal(failed.converged, false); assert.equal(failed.gridSequence.reachedTarget, false);
    assert.match(failed.reason, /did not converge/);
    const bad = structuredClone(root); bad.checkpoint.restart.input.mach = .2;
    assert.throws(() => solve(input, { ...options, parentResult: bad }), /checkpoint Mach/);
  } finally { delete globalThis.__targetGridTest; }
});

test('stalled operating point refines at actual conditions before resuming; failed fine solves retain the coarse root', async () => {
  const normalization = { referenceChord: 1, solverLength: 1, referenceReynolds: 1e6, kernelReynolds: 1e6 };
  const opts = { reynolds: 1e6, ncrit: 4, edgeMatching: 'section-velocity', tripFractions: [[1, 1]] };
  const families = { euler: 1e-12, boundaryLayer: 1e-12, edgeMatching: 1e-12 };
  const cp = { families, restart: { input: { mach: .72, alpha: 2.1, weights: [Array(10).fill(1), Array(10).fill(1)] }, options: opts } };
  const root = { converged: false, stateConverged: true, mach: .72, actualAlpha: 2.1,
    mesh: { quality: { valid: true }, cells: [] }, checkpoint: cp, families,
    sourceCase: { ...input, mach: .72, alpha: 2.1, gridIntervals: 32 },
    ...normalization, initialization: {}, solverSettings: {} };
  const order = [], events = [];
  let failFine = false, failResume = false;
  globalThis.__retainedGridTest = {
    checkpointDataEqual: (a, b) => JSON.stringify(a) === JSON.stringify(b),
    coupledAssemblyConditions: () => ({ options: opts, normalization }),
    createCoupledStreamtubeBody: () => ({ bl: { hasFiniteBase: false } }),
    prepareCoupledGridSequence: args => {
      order.push('refine'); assert.equal(args.input.mach, .72); assert.equal(args.input.alpha, 2.1);
      assert.equal(args.sourceResult.converged, true);
      return { system: {}, transfer: { nominalGridIntervals: 64 } };
    },
    solveCoupledGridLevel: () => {
      order.push('fine-solve'); return { ...root, converged: !failFine, reason: 'iteration limit',
        solverInput: { bodies: [{}] }, linearDiagnostics: { solves: 3 } };
    },
  };
  try {
    const code = fs.readFileSync('src/euler/tests/streamtube-coupled-target-grid.js', 'utf8')
      .replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}} = globalThis.__retainedGridTest;`);
    const { solveCoupledTargetGrid: solve } = await import('data:text/javascript;base64,' + Buffer.from(code).toString('base64'));
    const options = { plan: targetMachGridPlan(input), solveCoarse: () => structuredClone(root),
      onStage: e => events.push(e),
      resumeOperatingPoint: fine => {
        order.push('continue'); assert.equal(fine.sourceCase.gridIntervals, 64);
        assert.equal(fine.mach, .72); assert.equal(fine.actualAlpha, 2.1);
        if (failResume) return { ...fine, converged: false, stateConverged: true, reason: 'stage limit' };
        const checkpoint = structuredClone(fine.checkpoint);
        checkpoint.restart.input.mach = .74; checkpoint.restart.input.alpha = 2.68;
        return { ...fine, checkpoint, converged: true, mach: .74, actualAlpha: 2.68,
          sourceCase: { ...input } };
      } };
    const solved = solve(input, options);
    assert.equal(solved.converged, true);
    assert.deepEqual(order, ['refine', 'fine-solve', 'continue']);
    assert.equal(events[0].gridStrategy, 'refine-retained-operating-point');
    const resumed = solve(input, { ...options, parentResult: root,
      solveCoarse: () => { throw new Error('must retain warm source'); } });
    assert.equal(resumed.converged, true);
    const warm64 = { ...root, sourceCase: { ...root.sourceCase, gridIntervals: 64 } };
    assert.equal(targetMachGridPlan({ ...input, gridIntervals: 128 }, { parentResult: warm64 }).coarseIntervals, 64);
    assert.equal(targetMachGridPlan({ ...input, gridIntervals: 128, reynolds: 9 }, { parentResult: warm64 }), null);
    failResume = true;
    const stalled = solve(input, options);
    assert.equal(stalled.converged, false); assert.equal(stalled.gridSequence.reachedTarget, false);
    assert.equal(stalled.gridSequence.actualGridIntervals, 64);
    assert.equal(stalled.actualAlpha, 2.1);
    failFine = true; order.length = 0;
    const failed = solve(input, options);
    assert.deepEqual(order, ['refine', 'fine-solve']);
    assert.equal(failed.sourceCase.gridIntervals, 32);
    assert.deepEqual(failed.checkpoint, root.checkpoint);
    assert.equal(failed.stateConverged, true);
  } finally { delete globalThis.__retainedGridTest; }
});
