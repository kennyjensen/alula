// Zero-step driver contracts with declared numerical stubs. No Jacobian,
// factorization, global Newton update, physical flow or mesh calculation runs.
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';

const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url), 'utf8');
const clone = structuredClone;
const families = { euler: .1, boundaryLayer: .2, edgeMatching: .3 };
const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
function checkpoint(coordinate, blUpdate = 'xfoil') {
  return { version: 1, families: clone(families), restart: {
    input: { mach: .2, bodies: [{ leadingIndex: 1 }], wakeGeometry: 'independent-banks' },
    options: { ncrit: 4, reynolds: 1e6, transitionMode: 'automatic', tripFractions: [[1, 1]], transitionState: [2, 3],
      hkFloorLinearization: 'native' },
    initialEuler: { x: [1], nodes: clone(nodes), undisplacedNodes: clone(nodes) }, initialBL: [.1, .2, .3] },
    continuation: { fractions: [[0, 1]], lastRedistributedStagnation: [.2], preferredOrdering: 'amd', pivotTolerance: .001,
      iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing',
      ...(blUpdate === undefined ? {} : { blUpdate }), projectionGeometry: blUpdate === 'xfoil' ? 'boundary-increment' : 'fixed',
      ...(coordinate === undefined ? {} : { shearCoordinate: coordinate }) } };
}
const controls = { maxIterations: 0, tolerance: 1e-10, iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing' };
let serial = 0;
async function harness() {
  const calls = { constructs: [], evaluations: 0, forbidden: [] }, key = `__shearCheckpoint${++serial}`;
  const forbidden = name => () => { calls.forbidden.push(name); throw new Error(`Forbidden numerical call: ${name}`); };
  const stubs = {
    createCoupledStreamtubeBody(input, options) {
      calls.constructs.push({ input: clone(input), options: clone(options) });
      const initial = Float64Array.from([...(options.initialEuler?.x ?? [1]), ...(options.initialBL ?? [.1, .2, .3])]);
      return { initial, ne: 1,
        euler: { layout: { densityCount: 1, bodies: input.bodies, rows: [{ kind: 'streamwise' }] },
          conditions: { pressureScale: 1 },
          setDisplacement() {}, adoptGeometry: x => x.slice() },
        bl: { snapshotActive: () => [2, 3], thicknesses: () => ({}), surfaces: [{}, {}],
          updateActive: () => ({ changed: false, changes: [] }) },
        admissible: () => true,
        admissibleValue(state) { return this.evaluate(state); },
        evaluate() { calls.evaluations++; return { residual: Float64Array.of(.1, .2, .3, 0), families: clone(families),
          outer: { nodes: clone(nodes), undisplacedNodes: clone(nodes), stagnation: [.2], diagnostics: { maxMach: .2 } } }; },
        jacobian: forbidden('Jacobian'),
      };
    },
    coupledStreamtubeResult: (system, state, details) => ({ ...details, x: Array.from(state) }),
    createIterationProgress: () => ({ seed() {} }),
    physicalIterationVector: () => Float64Array.of(0),
    residualMerit: () => 0,
    captureStreamtubeInletFractions: () => [[0, 1]],
    redistributeStreamtubeTangentially: nodes => ({ nodes, solution: { pairs: [], relativeResidual: 0 }, maxDisplacement: 0 }),
    assertConvexStreamtubeGrid() {},
  };
  for (const name of ['coupledStreamtubeTripEvents', 'proposeCoupledDensityNewton', 'proposeCoupledXfoilBLUpdate',
    'proposeLogarithmicShearUpdate', 'respondToCoupledProjectionGeometry', 'prepareCoupledTransitionProfileTrial',
    'solveSparseDirect', 'createStreamtubeStationOrdering', 'solveCoupledLinearSystem', 'adjustStreamtubeInlets',
    'dekinkStreamtubeInteriors', 'prepareConvexWakeGridUpdate', 'requireCoupledResidualDecrease']) stubs[name] = forbidden(name);
  globalThis[key] = stubs;
  // Both aliased proposal imports are distinct stubs; never import solver numerics.
  const text = source.replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => {
    const bindings = names.split(',').map(name => name.trim().split(/\s+as\s+/).at(-1)).join(', ');
    return `const { ${bindings} } = globalThis.${key};`;
  });
  try { return { ...(await import('data:text/javascript;base64,' + Buffer.from(text + '\n//' + serial).toString('base64'))), calls }; }
  finally { delete globalThis[key]; }
}

test('fresh zero-step solves keep ordinary and native-Hk startup linear unless log is explicit', async () => {
  for (const coordinate of [undefined, 'linear', 'logarithmic']) for (const hkFloorLinearization of [undefined, 'native']) {
    const h = await harness(), cp = checkpoint(coordinate), published = [];
    const options = { ...controls, ...cp.restart.options, hkFloorLinearization, initialEuler: cp.restart.initialEuler,
      initialBL: cp.restart.initialBL, blUpdate: 'xfoil', projectionGeometry: 'boundary-increment',
      ...(coordinate === undefined ? {} : { shearCoordinate: coordinate }), onCheckpoint: cp => published.push(cp) };
    const before = clone({ input: cp.restart.input, options: Object.fromEntries(Object.entries(options).filter(([, v]) => typeof v !== 'function')) });
    const result = h.solveCoupledStreamtubeIses(cp.restart.input, options);
    const expected = coordinate === 'logarithmic' ? 'logarithmic' : undefined;
    assert.equal(result.shearCoordinate, expected);
    assert.equal(result.checkpoint.continuation.shearCoordinate, expected);
    assert.equal(result.checkpoint.continuation.eliminateAmplification, true);
    assert.equal(result.initialRedistribution.accepted, true);
    assert.equal(published.length, 1); assert.equal(published[0].continuation.shearCoordinate, expected);
    assert.deepEqual(h.calls.forbidden, []); assert.equal(result.linearDiagnostics.solves, 0);
    assert.equal(result.history.length, 1); assert.equal(result.history[0].iteration, 0);
    for (const call of h.calls.constructs) assert.equal(Object.hasOwn(call.options, 'shearCoordinate'), false);
    assert.deepEqual({ input: cp.restart.input, options: Object.fromEntries(Object.entries(options).filter(([, v]) => typeof v !== 'function')) }, before);
  }
});

test('zero-step checkpoint resume preserves log markers while omitted and explicit linear remain legacy-compatible', async () => {
  for (const coordinate of [undefined, 'linear', 'logarithmic']) {
    const h = await harness(), cp = checkpoint(coordinate), before = clone(cp), published = [];
    const result = h.solveCoupledStreamtubeIses(undefined, { ...controls, resume: cp, onCheckpoint: cp => published.push(cp) });
    assert.equal(result.checkpoint.continuation.shearCoordinate, coordinate === 'logarithmic' ? coordinate : undefined);
    assert.equal(result.shearCoordinate, coordinate === 'logarithmic' ? coordinate : undefined);
    assert.equal(result.checkpoint.restart.options.hkFloorLinearization, 'native');
    assert.equal(result.checkpoint.restart.options.shearCoordinate, undefined);
    assert.equal(result.checkpoint.continuation.eliminateAmplification, undefined);
    assert.deepEqual(result.x, [...cp.restart.initialEuler.x, ...cp.restart.initialBL]);
    assert.deepEqual(result.checkpoint.restart.initialEuler.nodes, cp.restart.initialEuler.nodes);
    assert.deepEqual(result.checkpoint.restart.initialEuler.undisplacedNodes, cp.restart.initialEuler.undisplacedNodes);
    assert.equal(published.length, 1); assert.equal(h.calls.constructs.length, 1);
    assert.deepEqual(cp, before); assert.deepEqual(h.calls.forbidden, []);
  }
});

test('resume coordinate disagreement rejects before constructing a system', async () => {
  for (const [saved, requested] of [[undefined, 'logarithmic'], ['linear', 'logarithmic'], ['logarithmic', 'linear'],
    ['logarithmic', null], [undefined, null]]) {
    const h = await harness(), cp = checkpoint(saved), before = clone(cp);
    assert.throws(() => h.solveCoupledStreamtubeIses(undefined, { ...controls, resume: cp, shearCoordinate: requested }), /shear-coordinate controls do not match/);
    assert.equal(h.calls.constructs.length, 0); assert.deepEqual(h.calls.forbidden, []); assert.deepEqual(cp, before);
  }
});

test('malformed coordinate markers and logarithmic Giles updates reject before constructing a system', async () => {
  for (const coordinate of [null, false, 0, '', 'wrong', {}, []]) for (const resumed of [false, true]) {
    const h = await harness(), cp = checkpoint(coordinate);
    assert.throws(() => h.solveCoupledStreamtubeIses(resumed ? undefined : cp.restart.input,
      resumed ? { ...controls, resume: cp } : { ...controls, shearCoordinate: coordinate, blUpdate: 'xfoil' }), /Invalid coupled ISES iteration controls/);
    assert.equal(h.calls.constructs.length, 0); assert.deepEqual(h.calls.forbidden, []);
  }
  for (const resumed of [false, true]) {
    const h = await harness(), cp = checkpoint('logarithmic', 'giles');
    assert.throws(() => h.solveCoupledStreamtubeIses(resumed ? undefined : cp.restart.input,
      resumed ? { ...controls, resume: cp } : { ...controls, shearCoordinate: 'logarithmic' }), /Invalid coupled ISES iteration controls/);
    assert.equal(h.calls.constructs.length, 0); assert.deepEqual(h.calls.forbidden, []);
  }
});

test('zero-step log resumes still enforce exact residual replay before publishing a checkpoint', async () => {
  const h = await harness(), cp = checkpoint('logarithmic'), before = clone(cp);
  cp.families.euler += .01;
  let published = 0;
  assert.throws(() => h.solveCoupledStreamtubeIses(undefined, { ...controls, resume: cp, onCheckpoint: () => published++ }), /residual does not replay exactly/);
  assert.equal(published, 0); assert.equal(h.calls.constructs.length, 1); assert.deepEqual(h.calls.forbidden, []);
  assert.deepEqual(cp.restart, before.restart); assert.deepEqual(cp.continuation, before.continuation);
});
