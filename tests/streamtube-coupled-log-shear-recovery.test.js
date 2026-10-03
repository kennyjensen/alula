import test from 'node:test';
import assert from 'node:assert/strict';
import { planCoupledLogarithmicShearRecovery as planRecovery } from '../src/euler/streamtube-coupled-log-shear-recovery.js';

const families = ['euler', 'boundaryLayer', 'edgeMatching'];
const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
const startupPlan = { kind: 'ncrit-transition-startup', sourceNcrit: 4, targetNcrit: 9,
  shearCoordinate: 'linear', sourceShearCoordinate: 'linear', hkFloorLinearization: 'native', sourceIterations: 40 };
const controls = { startupPlan, startupAttempt: 2, maxIterations: 40, tolerance: 1e-10 };
function fixture(maxIterations = 40) {
  const history = Array.from({ length: maxIterations + 1 }, (_, iteration) => ({ iteration,
    euler: .2 + (maxIterations - iteration) * .01,
    boundaryLayer: 2 + (maxIterations - iteration) * .01,
    edgeMatching: .4 + (maxIterations - iteration) * .01,
    step: iteration ? .001 : 0, backtracks: 0, rejections: [] }));
  const f = Object.fromEntries(families.map(k => [k, history.at(-1)[k]]));
  const checkpoint = { version: 1, families: { ...f }, restart: {
    input: { mach: .2, bodies: [{ leadingIndex: 2, trailingIndex: 8 }], hybrid: { ismom: 4 },
      wakeGeometry: 'independent-banks', wakeDisplacementMotion: 'te-center' },
    options: { reynolds: 1e6, ncrit: 4, transitionMode: 'automatic', tripFractions: [[1, 1]],
      transitionState: [2, 3], hkFloorLinearization: 'native' },
    initialEuler: { x: Float64Array.of(0, 0), nodes: structuredClone(nodes), undisplacedNodes: structuredClone(nodes) },
    initialBL: Float64Array.of(.03, 1, 2, 1) }, continuation: {
      fractions: [[0, .5, 1]], lastRedistributedStagnation: [.2], iterationGeometry: 'ises-sampled',
      stepAcceptance: 'admissible', stagnationLimiter: 'listing', blUpdate: 'xfoil', projectionGeometry: 'boundary-increment',
      linearOrdering: 'station-auto', stationFallback: true, preferredOrdering: 'amd', pivotTolerance: .001 } };
  return { converged: false, reason: 'iteration limit', initialRedistribution: { accepted: true }, families: f,
    mesh: { quality: { valid: true, invalidCells: [], minArea: .01, minCornerSine: .5 } }, history, checkpoint,
    conditions: { mach: .2, reynolds: 1e6, ncrit: 4, transitionMode: 'automatic', hkFloorLinearization: 'native' },
    x: Float64Array.of(0, 0, .03, 1, 2, 1), residual: Float64Array.of(f.euler, -f.boundaryLayer, f.edgeMatching, 0, 0, 0),
    flow: { nodes: structuredClone(nodes), undisplacedNodes: structuredClone(nodes) },
    linearDiagnostics: { solves: maxIterations, refinements: 0, pivotRecoveries: 0, maxRelativeResidual: 1e-14,
      iterations: Array.from({ length: maxIterations }, (_, i) => ({ iteration: i + 1 })) }, lastRejectedStep: null };
}

test('a complete stalled linear second-start source transfers only its cloned policy marker', () => {
  const result = fixture(), before = structuredClone(result), p = planRecovery(result, controls);
  assert(p); assert.equal(p.kind, 'coupled-logarithmic-shear-recovery');
  assert.equal(p.originalIterations, 40); assert.equal(p.additionalIterations, 40); assert.equal(p.maximumTotalIterations, 80);
  assert.equal(p.sourceShearCoordinate, 'linear'); assert.equal(p.shearCoordinate, 'logarithmic');
  const expected = structuredClone(result.checkpoint); expected.continuation.shearCoordinate = 'logarithmic';
  assert.deepEqual(p.resume, expected); assert.deepEqual(result, before);
  p.resume.restart.initialBL[0] = 17; p.resume.restart.initialEuler.nodes[0][0][0].x = 99;
  assert.deepEqual(result, before, 'The transferred continuation must not own the source arrays.');
});

test('detected early stagnation still permits the existing logarithmic recovery', () => {
  const result = fixture(3);
  result.reason = 'Coupled progress stalled: negligible-state-change.';
  result.progressControl = { stopReason: result.reason, recoveries: 0 };
  const plan = planRecovery(result, controls);
  assert(plan); assert.equal(plan.originalIterations, 3); assert.equal(plan.maximumTotalIterations, 43);
  assert.equal(plan.resume.continuation.shearCoordinate, 'logarithmic');
  result.progressControl.stopReason = 'unrelated'; assert.equal(planRecovery(result, controls), null);
});

test('the extra chunk respects public caps and no already-log checkpoint can receive it twice', () => {
  const budgets = [];
  for (const maxIterations of [1, 5, 20, 40, 80]) {
    const p = planRecovery(fixture(maxIterations), { ...controls, maxIterations,
      startupPlan: { ...startupPlan, sourceIterations: maxIterations } });
    assert(p, `cap ${maxIterations}`);
    assert.equal(p.additionalIterations, Math.min(40, maxIterations));
    assert.equal(p.maximumTotalIterations, maxIterations + Math.min(40, maxIterations));
    budgets.push({ maxIterations, additionalIterations: p.additionalIterations });
  }
  const repeated = fixture(); repeated.checkpoint.continuation.shearCoordinate = 'logarithmic';
  assert.equal(planRecovery(repeated, controls), null);
  for (const maxIterations of [0, -1, 1.5, Infinity, NaN]) assert.equal(planRecovery(fixture(), { ...controls, maxIterations }), null);
});

test('successful linear roots, first attempts and unrelated startup routes remain unchanged', () => {
  const success = fixture(); success.converged = true;
  const before = structuredClone(success); assert.equal(planRecovery(success, controls), null); assert.deepEqual(success, before);
  for (const patch of [{ startupAttempt: 1 }, { startupAttempt: 3 }, { startupPlan: null },
    { startupPlan: { ...startupPlan, kind: 'other' } }, { startupPlan: { ...startupPlan, sourceNcrit: 5 } },
    { startupPlan: { ...startupPlan, shearCoordinate: 'logarithmic' } },
    { startupPlan: { ...startupPlan, targetNcrit: Infinity } }, { tolerance: 0 }, { tolerance: NaN }])
    assert.equal(planRecovery(fixture(), { ...controls, ...patch }), null);
  const legacy = fixture(), explicit = fixture(); explicit.checkpoint.continuation.shearCoordinate = 'linear';
  assert.deepEqual(planRecovery(legacy, controls).resume, planRecovery(explicit, controls).resume);
  const invalidPolicy = fixture(); invalidPolicy.checkpoint.continuation.shearCoordinate = 'invalid';
  assert.throws(() => planRecovery(invalidPolicy, controls), /shear-coordinate policy/);
});

test('termination, family, history and complete checkpoint mismatches forbid transfer', () => {
  const mutations = [
    ['non-iteration failure', r => { r.reason = 'admissibility'; }],
    ['last rejection', r => { r.lastRejectedStep = { stage: 'admissibility' }; }],
    ['missing checkpoint', r => { delete r.checkpoint; }],
    ['checkpoint version', r => { r.checkpoint.version = 2; }],
    ['invalid grid flag', r => { r.mesh.quality.valid = false; }],
    ['invalid cell', r => { r.mesh.quality.invalidCells.push(0); }],
    ['no accepted initial chart', r => { r.initialRedistribution.accepted = false; }],
    ['short history', r => { r.history.pop(); }],
    ['renumbered history', r => { r.history[20].iteration = 21; }],
    ['history-family mismatch', r => { r.history.at(-1).boundaryLayer += 1; }],
    ['checkpoint-family mismatch', r => { r.checkpoint.families.euler += .1; }],
    ['nonfinite family', r => { r.families.euler = Infinity; }],
    ['already-small residual', r => { for (const k of families) r.families[k] = 1e-13; }],
    ['Euler vector mismatch', r => { r.x[0] += .1; }],
    ['BL vector mismatch', r => { r.x[2] += .1; }],
    ['residual vector length', r => { r.residual = r.residual.slice(1); }],
    ['residual-family maximum mismatch', r => { r.residual[1] = -4; }],
    ['physical nodes mismatch', r => { r.flow.nodes[0][0][0].x += .1; }],
    ['undisplaced nodes mismatch', r => { r.flow.undisplacedNodes[0][0][0].y += .1; }],
    ['wrong Ncrit', r => { r.checkpoint.restart.options.ncrit = 5; }],
    ['wrong BL policy', r => { r.checkpoint.continuation.blUpdate = 'giles'; }],
    ['missing transition phase', r => { delete r.checkpoint.restart.options.transitionState; }],
    ['noninteger transition phase', r => { r.checkpoint.restart.options.transitionState[0] = .2; }],
    ['missing inlet fractions', r => { delete r.checkpoint.continuation.fractions; }],
    ['nonmonotone inlet fractions', r => { r.checkpoint.continuation.fractions[0] = [0, .7, .6]; }],
    ['missing last redistribution', r => { delete r.checkpoint.continuation.lastRedistributedStagnation; }],
    ['nonfinite last redistribution', r => { r.checkpoint.continuation.lastRedistributedStagnation[0] = NaN; }],
  ];
  for (const [name, mutate] of mutations) {
    const r = fixture(); mutate(r); const before = structuredClone(r);
    assert.equal(planRecovery(r, controls), null, name); assert.deepEqual(r, before, name + ' source changed');
  }
});
