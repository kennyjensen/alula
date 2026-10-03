import { coupledConvergenceSatisfied } from '../src/euler/streamtube-coupled-convergence.js';
import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { planCoupledGridSequence, prepareCoupledGridSequence } from '../src/euler/streamtube-coupled-grid-sequence.js';

const input = (weights, nx = 8) => ({ bodies: Array.from({ length: weights.length - 1 }, (_, element) =>
  ({ element, leadingIndex: 2, trailingIndex: 5 })), weights, outerLower: Array.from({ length: nx + 1 }, (_, x) => ({ x, y: -2 })) });

test('the measured NLR dimensions produce the qualified nested grid without hardcoded body or passage counts', () => {
  // Actual saved coarse mass weights; the planner sees only dimensions and
  // weights, not an airfoil name, Mach, Reynolds number or case identifier.
  const weights = [[729, 243, 81, 27, 9, 3, .25, .25, .25, .25],
    [.0556302334890786, .0556302334890786, .0556302334890786, .0556302334890786,
      .6234898018587335, .9009688679024191, 1, .9009688679024191, .6234898018587339,
      .05563023348907852, .05563023348907852, .05563023348907852, .05563023348907852],
    [.25, .25, .25, .25, 3, 9, 27, 81, 243, 729]];
  const args = { sourceInput: input(weights, 130), input: input([Array(12).fill(1), Array(15).fill(1), Array(12).fill(1)], 231),
    sourceGridIntervals: 16, requestedGridIntervals: 32 }, before = structuredClone(args);
  const p = planCoupledGridSequence(args);
  assert.equal(p.streamwiseFactor, 2); assert.deepEqual(p.child, { nx: 260, tubes: [12, 15, 12] });
  assert.equal(p.child.nx * p.child.tubes.reduce((a, b) => a + b, 0), 10140);
  assert.equal(p.nodeCount, 10962); assert.equal(p.requestedGridRetained, false);
  assert.deepEqual(p.normalSubdivisions[0], [3, 1, 1, 1, 1, 1, 1, 1, 1, 1]);
  assert.deepEqual(p.normalSubdivisions[2], [1, 1, 1, 1, 1, 1, 1, 1, 1, 3]);
  assert.deepEqual(args, before);
});

test('generic sequencing rounds nominal resolution upward, never merges tubes, and partitions each parent mass', () => {
  const sourceInput = input([[10, 3, 1], [1, 1, 1, 1], [2, 9, 4]]);
  const requested = input([[1, 1], Array(9).fill(1), Array(8).fill(1)]);
  const p = planCoupledGridSequence({ sourceInput, input: requested, sourceGridIntervals: 16, requestedGridIntervals: 47 });
  assert.equal(p.streamwiseFactor, 3); assert.equal(p.nominalGridIntervals, 48);
  assert.deepEqual(p.child.tubes, [3, 9, 8]); assert.deepEqual(p.normalSubdivisions[0], [1, 1, 1]);
  sourceInput.weights.forEach((weights, g) => weights.forEach((mass, j) => {
    const n = p.normalSubdivisions[g][j]; assert.ok(n >= 1 && n <= 4);
    assert.ok(Math.abs(Array(n).fill(mass / n).reduce((a, b) => a + b, 0) - mass) < 1e-14);
  }));
  const equal = p.normalSubdivisions[1]; assert.deepEqual(equal, [3, 2, 2, 2], 'equal mass ties are deterministic');
});

test('unsupported multilevel jumps, normal jumps, bad dimensions, no-op and node-budget excess fail before numerical construction', () => {
  const sourceInput = input([[2, 1], [1, 2]]), fine = input([Array(3).fill(1), Array(3).fill(1)]);
  const args = { sourceInput, input: fine, sourceGridIntervals: 16, requestedGridIntervals: 32 };
  assert.throws(() => planCoupledGridSequence({ ...args, requestedGridIntervals: 128 }), /converged intermediate level/);
  assert.throws(() => planCoupledGridSequence({ ...args, input: input([Array(9).fill(1), Array(3).fill(1)]) }), /parent tube/);
  assert.throws(() => planCoupledGridSequence(args, { maxNodes: 10 }), /node budget/);
  assert.throws(() => planCoupledGridSequence(args, { maxNodes: 50001 }), /node budget/);
  assert.throws(() => planCoupledGridSequence({ ...args, input: sourceInput, requestedGridIntervals: 16 }), /does not refine/);
  assert.throws(() => planCoupledGridSequence({ ...args, sourceGridIntervals: NaN }), /explicit valid/);
  const broken = structuredClone(sourceInput); broken.weights[0][0] = -1;
  assert.throws(() => planCoupledGridSequence({ ...args, sourceInput: broken }), /passage dimensions/);
});

// Orchestration-contract tests replace numerical constructors/refinement.
// The real local finite-base refiner and exact NLR source/child replay are
// validated separately; no synthetic state here is claimed to solve a PDE.
let serial = 0;
async function harness({ partial = false } = {}) {
  const calls = { constructs: 0, replays: 0, refinements: 0, profiles: 0 }, flags = {};
  const physical = { reynolds: 2.7e6, ncrit: 4, edgeMatching: 'section-velocity', transitionMode: 'automatic' };
  const options = { ...physical, tripFractions: [[1, 1], [1, 1]], transitionState: [2, 3, 4, 5] };
  const body = element => ({ element, leadingIndex: 2, trailingIndex: 5,
    trailingEdge: { kind: 'finite-base', upperIndex: 0, lowerIndex: 3 },
    points: [{ x: element, y: .01 }, { x: element - 1, y: 0 }, { x: element, y: -.01 }, { x: element, y: 0 }, { x: element, y: .01 }] });
  const a = { ...input([[10, 2], [1, 3], [2, 10]]), bodies: [body(1), body(0)],
    alpha: 6, mach: .185, wakeGeometry: 'independent-banks', wakeOutlet: 'banks',
    flowModel: 'compressible', streamwiseMode: 'isentropic', geometryDomain: 'positive-simple' };
  const b = structuredClone(a); b.weights = [[1, 1, 1], [1, 1, 1], [1, 1, 1]];
  const families = partial ? { euler: .01, boundaryLayer: .02, edgeMatching: .03 }
    : { euler: 1e-13, boundaryLayer: 2e-13, edgeMatching: 3e-13 };
  const residual = Float64Array.from(Object.values(families)), state = Float64Array.of(.2, .4, .03, 1, 2, 1);
  const nodes = [[[{ x: 0, y: 0 }, { x: 1, y: 0 }], [{ x: 0, y: 1 }, { x: 1, y: 1 }]]];
  const euler = i => ({ conditions: { mach: i.mach, alpha: i.alpha, flowModel: i.flowModel, streamwiseMode: i.streamwiseMode },
    layout: { bodies: i.bodies, independentWakeBanks: i.wakeGeometry === 'independent-banks' },
    decode: () => ({ nodes: structuredClone(nodes) }) });
  const system = (i, phase = options.transitionState, refined = false) => ({ initial: state.slice(), ne: 2,
    conditions: { ...physical }, euler: euler(i),
    bl: { transitionMode: 'automatic', trips: [[1, 1], [1, 1]], snapshotActive: () => phase.slice() },
    evaluate: () => { calls.replays++; return { residual: flags.badReplay ? Float64Array.of(1, 2, 3) : residual.slice(),
      families: refined ? { euler: .1, boundaryLayer: .2, edgeMatching: .3 } : { ...families }, outer: { nodes: structuredClone(nodes) } }; },
    admissible: () => !flags.invalid });
  const sourceSystem = system(a), saved = { input: a, options, initialEuler: { x: state.slice(0, 2),
    nodes: structuredClone(nodes), undisplacedNodes: structuredClone(nodes) }, initialBL: state.slice(2) };
  const sourceResult = { converged: true, residual, families, mesh: { quality: { valid: true } }, checkpoint: { restart: saved } };
  if (partial) Object.assign(sourceResult, { converged: false, reason: 'iteration limit', status: 'unconverged',
    lastRejectedStep: null, initialRedistribution: { accepted: true }, x: state.slice(),
    history: [{ iteration: 20, residual: .03, ...families }], boundaryLayer: { transitionState: options.transitionState.slice() },
    checkpoint: { version: 1, restart: saved, families: { ...families }, continuation: {
      fractions: [[0, .5, 1], [0, .5, 1]], lastRedistributedStagnation: [.5, .5], preferredOrdering: 'amd',
      pivotTolerance: .001, iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing',
      blUpdate: 'xfoil', projectionGeometry: 'boundary-increment' } } });
  const stubs = { ...shearPolicy, coupledConvergenceSatisfied,
    prepareRefinedCoupledGridProfile(prepared) { calls.profiles++; return { ...prepared, nativeProfilePrepared: true }; },
    createStreamtubeBodySystem(i) { calls.constructs++; return euler(i); },
    initialStreamtubeDisplacement: () => ({ surfaces: [], wakes: [] }),
    createCoupledStreamtubeBody: i => system(i),
    streamtubeMeshSnapshot: () => ({ quality: { valid: !flags.invalidMesh }, cells: Array(24).fill({}), initialization: {} }),
    refineCoupledStreamtubeBody(i, parent, plan) {
      calls.refinements++; calls.plan = structuredClone(plan);
      const child = structuredClone(i); child.gridSpacing = {};
      child.bodies.forEach(v => { v.leadingIndex *= plan.streamwiseFactor; v.trailingIndex *= plan.streamwiseFactor; });
      return { input: child, system: system(child, [4, 7, 8, 10], true), options: { ...options },
        initialEuler: saved.initialEuler, initialBL: saved.initialBL, diagnostics: { finiteBaseTransfer: { present: true } } };
    },
  };
  const key = `__gridSequence${++serial}`; globalThis[key] = stubs;
  const text = fs.readFileSync(new URL('../src/euler/streamtube-coupled-grid-sequence.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}}=globalThis.${key};`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(text).toString('base64'));
  return { ...module, calls, flags, args: { sourceSystem, sourceResult, input: b,
    options: { ...options }, initialEuler: saved.initialEuler, sourceGridIntervals: 16, requestedGridIntervals: 32 },
    release: () => delete globalThis[key] };
}

test('native profile preparation is explicit and the omitted or interpolate predictor leaves the original seed unchanged', async () => {
  const h = await harness();
  try {
    const omitted = h.prepareCoupledGridSequence(h.args);
    const interpolate = h.prepareCoupledGridSequence(h.args, { blPredictor: 'interpolate' });
    assert.equal(h.calls.profiles, 0);
    assert.deepEqual(omitted.initialEuler, interpolate.initialEuler); assert.deepEqual(omitted.initialBL, interpolate.initialBL);
    assert.deepEqual(omitted.transfer, interpolate.transfer);
    const native = h.prepareCoupledGridSequence(h.args, { blPredictor: 'xfoil-mrchdu' });
    assert.equal(h.calls.profiles, 1); assert.equal(native.nativeProfilePrepared, true);
  } finally { h.release(); }
  const invalid = await harness();
  try {
    for (const blPredictor of ['unknown', null, false])
      assert.throws(() => invalid.prepareCoupledGridSequence(invalid.args, { blPredictor }), /BL predictor/);
    assert.equal(invalid.calls.constructs, 0); assert.equal(invalid.calls.refinements, 0); assert.equal(invalid.calls.profiles, 0);
  } finally { invalid.release(); }
});

test('accepted orchestration transfers complete state with updated phases, unchanged physics and explicit replacement grid', async () => {
  const h = await harness(), initial = h.args.sourceSystem.initial.slice(), requested = structuredClone(h.args.input);
  const p = h.prepareCoupledGridSequence(h.args);
  assert.equal(h.calls.refinements, 1); assert.equal(h.calls.plan.streamwiseFactor, 2);
  assert.equal(p.initialization.thicknessFactor, 1); assert.equal(p.initialization.flowSolved, false);
  assert.equal(p.transfer.requestedGridRetained, false); assert.equal(p.transfer.initialGuessOnly, true);
  assert.deepEqual(p.options.transitionState, [4, 7, 8, 10]);
  assert.deepEqual(p.transfer.sourceTransitionState, [2, 3, 4, 5]);
  assert.deepEqual(p.input.bodies.map(v => v.trailingEdge), h.args.input.bodies.map(v => v.trailingEdge));
  assert.deepEqual(p.input.gridSpacing.surfaceIntervalsByElement, [{ element: 1, intervals: 6 }, { element: 0, intervals: 6 }]);
  assert.deepEqual(h.args.input, requested); assert.deepEqual(h.args.sourceSystem.initial, initial);
  assert.equal(p.transfer.operations.globalLinearSolves, 0); assert.equal(p.transfer.operations.nativeBLInitializations, 0);
  h.release();
});

test('complete-state sequencing needs target geometry and conditions, not an unused target Euler gas state', async () => {
  const h = await harness();
  try {
    delete h.args.initialEuler;
    const result = h.prepareCoupledGridSequence(h.args);
    assert.equal(result.transfer.requestedGridRetained, false);
    assert.equal(result.transfer.thicknessFactor, 1);
    assert.equal(h.calls.refinements, 1);
    assert.equal(h.calls.constructs, 1, 'Only target condition/geometry metadata is constructed');
    assert.equal(result.transfer.operations.nativeBLInitializations, 0);
  } finally { h.release(); }
});

test('source convergence, replay, target physical controls, material identities, trips and final admissibility remain mandatory', async () => {
  const changes = [
    [h => { h.args.sourceResult.converged = false; }, /converged accepted/],
    [h => { h.args.sourceResult.families.euler = .01; }, /converged accepted/],
    [h => { h.args.input.mach = .2; }, /physical\/model input mach/],
    [h => { h.args.options.reynolds *= 2; }, /option reynolds/],
    [h => { h.args.options.ncrit = 9; }, /option ncrit/],
    [h => { h.args.options.hkFloorLinearization = 'native'; }, /Hk-floor linearization policy/],
    [h => { h.args.sourceSystem.conditions.hkFloorLinearization = 'native'; }, /Hk-floor linearization policy/],
    [h => { h.args.options.tripFractions = [[.9, 1], [1, 1]]; }, /material trips/],
    [h => { h.args.input.bodies[0].trailingEdge.upperIndex = 1; }, /trailing-edge topology/],
    [h => { h.args.input.bodies.reverse(); }, /passage order/],
    [h => { h.args.input.bodies[0].points[0].y += .01; }, /solid contour/],
    [h => { h.args.sourceSystem.initial[0] += .1; }, /source checkpoint/],
    [h => { h.flags.badReplay = true; }, /replay exactly/],
    [h => { h.flags.invalidMesh = true; }, /source mesh/],
    [h => { h.flags.invalid = true; }, /replay exactly/],
  ];
  for (const [change, expression] of changes) {
    const h = await harness(); change(h);
    assert.throws(() => h.prepareCoupledGridSequence(h.args), expression);
    assert.equal(h.calls.refinements, 0); h.release();
  }
});

test('unaccepted source and oversize planning reject before constructing any target system', () => {
  assert.throws(() => prepareCoupledGridSequence({ sourceResult: { converged: false } }), /converged accepted/);
});

test('partial complete-state refinement is opt-in and explicitly remains an unconverged initial guess', async () => {
  const h = await harness({ partial: true }), before = structuredClone(h.args.sourceResult);
  try {
    assert.throws(() => h.prepareCoupledGridSequence(h.args), /converged accepted/);
    assert.equal(h.calls.constructs, 0);
    const p = h.prepareCoupledGridSequence(h.args, { allowPartialSource: true });
    assert.equal(p.transfer.method, 'partial-coupled-grid-sequence');
    assert.equal(p.transfer.sourceConverged, false); assert.equal(p.transfer.sourceReason, 'iteration limit');
    assert.equal(p.transfer.sourceIteration, 20); assert.deepEqual(p.transfer.sourceFamilies, before.families);
    assert.equal(p.transfer.initialGuessOnly, true); assert.equal(p.transfer.equationsChanged, false);
    assert.equal(p.transfer.requestedGridRetained, false); assert.equal(p.transfer.maintenanceHistoryInherited, false);
    assert.equal(p.initialization.flowSolved, false); assert.match(p.initialization.method, /unconverged accepted/);
    assert.equal(p.transfer.thicknessFactor, 1); assert.equal(p.transfer.operations.globalNewtonUpdates, 0);
    assert.equal(p.transfer.operations.globalLinearSolves, 0);
    assert.deepEqual(h.args.sourceResult, before); assert.equal(h.calls.refinements, 1);
  } finally { h.release(); }
});

test('partial source accepts each saved ordering policy without changing the complete source or inheriting factor history', async () => {
  for (const policy of [undefined, 'auto', 'station', 'station-auto', 'aligned-auto']) {
    const h = await harness({ partial: true });
    try {
      const c = h.args.sourceResult.checkpoint.continuation;
      if (policy !== undefined) c.linearOrdering = policy;
      if (policy === 'station-auto') c.stationFallback = true;
      const before = structuredClone(h.args.sourceResult);
      const prepared = h.prepareCoupledGridSequence(h.args, { allowPartialSource: true });
      assert.equal(prepared.transfer.method, 'partial-coupled-grid-sequence');
      assert.equal(prepared.transfer.maintenanceHistoryInherited, false);
      assert.equal(Object.hasOwn(prepared.options, 'linearOrdering'), false);
      assert.equal(Object.hasOwn(prepared.options, 'stationFallback'), false);
      assert.deepEqual(h.args.sourceResult, before); assert.equal(h.calls.refinements, 1);
    } finally { h.release(); }
  }
  for (const policy of ['unknown', false, 1]) {
    const h = await harness({ partial: true });
    try {
      h.args.sourceResult.checkpoint.continuation.linearOrdering = policy;
      assert.throws(() => h.prepareCoupledGridSequence(h.args, { allowPartialSource: true }), /accepted|converged/);
      assert.equal(h.calls.constructs, 0); assert.equal(h.calls.refinements, 0);
    } finally { h.release(); }
  }
  for (const [policy, fallback] of [['station-auto', null], ['station-auto', 'true'], ['station-auto', 1],
    ['station', true], ['auto', true], [undefined, true]]) {
    const h = await harness({ partial: true });
    try {
      const c = h.args.sourceResult.checkpoint.continuation;
      if (policy !== undefined) c.linearOrdering = policy;
      c.stationFallback = fallback;
      assert.throws(() => h.prepareCoupledGridSequence(h.args, { allowPartialSource: true }), /accepted|converged/);
      assert.equal(h.calls.constructs, 0); assert.equal(h.calls.refinements, 0);
    } finally { h.release(); }
  }
});

test('opting into partial refinement leaves the default converged-source result schema unchanged', async () => {
  const h = await harness();
  try {
    const exact = h.prepareCoupledGridSequence(h.args), optIn = h.prepareCoupledGridSequence(h.args, { allowPartialSource: true });
    for (const key of ['input', 'options', 'initialEuler', 'initialBL', 'mesh', 'transfer', 'initialization'])
      assert.deepEqual(optIn[key], exact[key]);
    assert.equal(Object.hasOwn(optIn.transfer, 'sourceConverged'), false);
    for (const value of [null, 1, 'true'])
      assert.throws(() => h.prepareCoupledGridSequence(h.args, { allowPartialSource: value }), /partial coupled grid-sequence policy/);
  } finally { h.release(); }
});

test('partial refinement excludes cancellations, rejected steps and incomplete or stale accepted checkpoints', async () => {
  const changes = [
    h => { h.args.sourceResult.reason = 'cancelled'; },
    h => { h.args.sourceResult.reason = 'Coupled ISES update rejected during admissibility'; },
    h => { h.args.sourceResult.lastRejectedStep = { stage: 'admissibility' }; },
    h => { delete h.args.sourceResult.lastRejectedStep; },
    h => { h.args.sourceResult.initialRedistribution.accepted = false; },
    h => { h.args.sourceResult.checkpoint.version = 2; },
    h => { delete h.args.sourceResult.checkpoint.continuation; },
    h => { h.args.sourceResult.checkpoint.continuation.fractions[0][1] = 1; },
    h => { h.args.sourceResult.checkpoint.continuation.lastRedistributedStagnation[0] = NaN; },
    h => { h.args.sourceResult.checkpoint.restart.initialEuler.x[0] = NaN; },
    h => { h.args.sourceResult.checkpoint.restart.initialEuler.nodes[0][0][0].x = NaN; },
    h => { h.args.sourceResult.checkpoint.restart.initialBL[0] = Infinity; },
    h => { delete h.args.sourceResult.checkpoint.restart.initialEuler.undisplacedNodes; },
    h => { h.args.sourceResult.checkpoint.restart.options.transitionState[0] = -1; },
    h => { h.args.sourceResult.boundaryLayer.transitionState[0]++; },
    h => { h.args.sourceResult.residual[0] = NaN; },
    h => { h.args.sourceResult.checkpoint.families.euler += .01; },
    h => { h.args.sourceResult.history.at(-1).boundaryLayer += .01; },
    h => { h.args.sourceResult.history.at(-1).iteration = 0; },
    h => { h.args.sourceResult.mesh.quality.valid = false; },
  ];
  for (const change of changes) {
    const h = await harness({ partial: true });
    try {
      change(h); assert.throws(() => h.prepareCoupledGridSequence(h.args, { allowPartialSource: true }), /accepted iteration-limit/);
      assert.equal(h.calls.constructs, 0); assert.equal(h.calls.refinements, 0);
    } finally { h.release(); }
  }
});

test('partial refinement still independently checks every source row, phase, physical cell and target condition', async () => {
  const changes = [
    [h => { h.args.sourceResult.x[0] += .01; }, /partial coupled source state/],
    [h => { h.args.sourceResult.checkpoint.restart.initialEuler.nodes[0][0][0].x += .01; }, /partial coupled source state/],
    [h => { h.args.sourceSystem.bl.snapshotActive = () => [2, 3, 4, 6]; }, /partial coupled source state/],
    [h => { h.flags.badReplay = true; }, /replay exactly/],
    [h => { h.flags.invalidMesh = true; }, /source mesh/],
    [h => { h.flags.invalid = true; }, /replay exactly/],
    [h => { h.args.input.mach = .2; }, /physical\/model input mach/],
    [h => { h.args.options.ncrit = 9; }, /option ncrit/],
    [h => { h.args.options.hkFloorLinearization = 'native'; }, /Hk-floor linearization policy/],
    [h => { h.args.options.tripFractions = [[.9, 1], [1, 1]]; }, /material trips/],
    [h => { h.args.input.bodies[0].points[0].x += .01; }, /solid contour/],
  ];
  for (const [change, expression] of changes) {
    const h = await harness({ partial: true });
    try {
      change(h); assert.throws(() => h.prepareCoupledGridSequence(h.args, { allowPartialSource: true }), expression);
      assert.equal(h.calls.refinements, 0);
    } finally { h.release(); }
  }
});

test('private preparation carries explicit shear policy independently of physical options and profile predictor',async()=>{
 for(const coordinate of[undefined,'linear','logarithmic']) for(const predictor of['interpolate','xfoil-mrchdu']){
  const h=await harness();
  h.args.sourceResult.checkpoint.continuation={blUpdate:'xfoil',...(coordinate===undefined?{}:{shearCoordinate:coordinate})};
  const before=structuredClone(h.args.sourceResult);
  try {const prepared=h.prepareCoupledGridSequence(h.args,{blPredictor:predictor});
    assert.equal(prepared.iterationControls.shearCoordinate,coordinate??'linear');
    assert.equal(prepared.options.shearCoordinate,undefined);
    assert.deepEqual(h.args.sourceResult,before);
  }finally{h.release();}
 }
 for(const coordinate of[null,false,{},'bad']){
  const h=await harness();h.args.sourceResult.checkpoint.continuation={blUpdate:'xfoil',shearCoordinate:coordinate};
  try{assert.throws(()=>h.prepareCoupledGridSequence(h.args),/shear-coordinate/);assert.equal(h.calls.refinements,0);}finally{h.release();}
 }
});
