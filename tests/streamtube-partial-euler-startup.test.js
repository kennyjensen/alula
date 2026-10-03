// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { initializeCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-initializer.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { selectCoupledEulerStartup } from '../src/euler/streamtube-coupled-startup.js';
import { initialStreamtubeDisplacement } from '../src/euler/streamtube-geometry.js';
import { limitStreamtubeGridStep, interpolateStreamtubeGridNodes } from '../src/geometry/streamtube-convex-step.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { finiteBaseBodyFixture } from './fixtures/finite-base-body.js';

const plain = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const seeds = new Map();
function fixture(finite = false) {
  if (!seeds.has(finite)) {
    const data = finite ? { ...finiteBaseBodyFixture({ bodySegments: 4, tubes: 2 }), mach: .03 }
      : intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
    data.bodies.forEach((b, element) => { b.element = element; });
    const source = createStreamtubeBodySystem(data);
    const solverInput = finite ? { ...data, displacement: initialStreamtubeDisplacement(source.layout, source.baseGeometry) } : data;
    const flow = solveStreamtubeIses(solverInput, { maxIterations: 0, iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible' });
    assert.equal(flow.initialRedistribution.accepted, true, flow.reason);
    assert.equal(flow.finalQuality.valid, true);
    const system = createStreamtubeBodySystem(flow.solverInput);
    const mesh = streamtubeMeshSnapshot({ system, nodes: flow.nodes });
    // Only the public stop envelope is manufactured. The packed state,
    // coordinates, gas equations, finite-base transfer and BL seed are real.
    // No 20-step Euler solve is needed to test this routing decision.
    flow.history = [...flow.history, { iteration: 20, residual: flow.diagnostics.residual, step: .1 }];
    const precursor = { status: 'unconverged', mach: data.mach, solverInput: flow.solverInput, flow, mesh,
      diagnostics: { iterations: 20, equationResidual: flow.diagnostics.residual, reason: 'iteration limit',
        solverStopReason: 'iteration limit', cells: mesh.cells.length } };
    const input = { elements: data.bodies.map(b => ({ points: b.points })), alpha: data.alpha, mach: data.mach,
      reynolds: 2.7e6, ncrit: 4, transitionMode: 'automatic', materialTrips: data.bodies.map(() => [1, 1]), gridIntervals: 4, gridTubes: 2 };
    seeds.set(finite, { input, precursor });
  }
  return structuredClone(seeds.get(finite));
}

// Mock only expensive orchestration endpoints: the Euler result is the
// saved real tiny seed above; the coupled driver evaluates zero updates.
// Crucially, ordinary BL initialization and all its domain checks run.
async function assembly(precursor, { initialize = initializeCoupledStreamtubeBody, failEuler } = {}) {
  const url = new URL('../src/euler/streamtube-coupled-assembly.js', import.meta.url), key = `partial-start-${Math.random()}`;
  const eulerCalls = [], initializationCalls = [], coupledCalls = [];
  globalThis[key] = {
    solveStreamtubeAssembly(input, options) {
      eulerCalls.push({ input: structuredClone(input), options });
      if (failEuler) throw failEuler;
      return precursor;
    },
    initializeCoupledStreamtubeBody(input, options, controls) {
      initializationCalls.push({ input: structuredClone(input), options: structuredClone({ ...options, initialEdgeVelocity: undefined }), controls });
      return initialize(input, options, controls);
    },
    solveCoupledStreamtubeIses(input, options) {
      coupledCalls.push({ input: structuredClone(input), options: structuredClone({ ...options,
        onIteration: undefined, onMesh: undefined, onCheckpoint: undefined }) });
      return solveCoupledStreamtubeIses(input, { ...options, maxIterations: 0 });
    },
  };
  let code = fs.readFileSync(url, 'utf8');
  for (const name of ['solveStreamtubeAssembly', 'initializeCoupledStreamtubeBody', 'solveCoupledStreamtubeIses']) {
    let matched = 0;
    code = code.replace(/import \{([^}]+)\} from '([^']+)';/g, (statement, names, source) => {
      const imports = names.split(',').map(v => v.trim());
      if (!imports.includes(name)) return statement;
      matched++;
      const retained = imports.filter(v => v !== name);
      return `${retained.length ? `import { ${retained.join(', ')} } from '${source}';\n` : ''}`
        + `const { ${name} } = globalThis[${JSON.stringify(key)}];`;
    });
    assert.equal(matched, 1, `Replace only the expensive ${name} seam and retain other real imports.`);
  }
  code = code.replace(/from '(\.[^']+)'/g, (_, relative) => `from '${new URL(relative, url).href}'`);
  try {
    const module = await import(`data:text/javascript;base64,${Buffer.from(code).toString('base64')}`);
    return { solve: module.solveCoupledStreamtubeAssembly, eulerCalls, initializationCalls, coupledCalls };
  } finally { delete globalThis[key]; }
}
const controls = { maxIterations: 1, maxStartupAttempts: 1, coarseStartup: false, tolerance: 1e-10 };

test('direct app startup shares a 40-update Euler budget across orders and respects explicit budgets', async () => {
  for (const [extra, expectedBudget, firstOrder] of [[{}, 20, false], [{ direct: true }, 40, true],
    [{ direct: true, eulerMaxIterations: 12 }, 12, true]]) {
    const { input, precursor } = fixture(), stop = new Error('inspect precursor controls');
    const d = await assembly(precursor, { initialize() { throw stop; } });
    assert.throws(() => d.solve(input, { ...controls, ...extra }), error => error === stop);
    assert.equal(d.eulerCalls.length, 1);
    assert.equal(d.eulerCalls[0].options.maxIterations, expectedBudget);
    assert.equal(d.eulerCalls[0].options.firstOrderStartup, firstOrder);
    assert.deepEqual(d.eulerCalls[0].input, input);
  }
});

test('partial adaptive Euler transfers its state without making temporary MCRIT the viscous target', async () => {
  const { input, precursor } = fixture();
  const data = { ...precursor.solverInput, streamwiseMode: 'momentum',
    upwind: { mucon: 1, mcrit: .85, boundary: { kind: 'unfiltered-first-two' } } };
  const flow = solveStreamtubeIses(data, { maxIterations: 0, iterationGeometry: 'ises-sampled',
    stepAcceptance: 'armijo', adaptiveMcrit: true, targetMcrit: .99 });
  flow.history.push({ iteration: 20, residual: flow.diagnostics.residual });
  precursor.flow = flow; precursor.solverInput = flow.solverInput;
  precursor.diagnostics.equationResidual = flow.diagnostics.residual;
  const before = plain(precursor), stop = new Error('inspect requested-law handoff');
  const driver = await assembly(precursor, { initialize(data, options) {
    assert.equal(data.upwind.mcrit, .99);
    assert.equal(data.upwind.mucon, 1);
    const source = createStreamtubeBodySystem(precursor.solverInput);
    assert.deepEqual(plain(options.initialEuler.x.slice(0, source.layout.densityCount)),
      plain(precursor.flow.x.slice(0, source.layout.densityCount)));
    assert.deepEqual(options.initialEuler.nodes, precursor.flow.nodes);
    throw stop;
  } });
  assert.throws(() => driver.solve(input, controls), error => error === stop);
  assert.deepEqual(plain(precursor), before);
});

test('an accepted iteration-limit Euler state enters ordinary automatic BL initialization without claiming convergence', async () => {
  const { input, precursor } = fixture(), original = plain(precursor), d = await assembly(precursor), stages = [], checkpoints = [];
  const result = d.solve(input, { ...controls, onStage: value => stages.push(value), onCheckpoint: value => checkpoints.push(value) });
  assert.equal(d.initializationCalls.length, 1);
  assert.equal(d.coupledCalls.length, 1);
  assert.equal(result.converged, false);
  assert.equal(result.physicalAcceptance, false);
  assert.deepEqual(result.initialization.euler, { iterations: 20, residual: precursor.diagnostics.equationResidual,
    cells: precursor.diagnostics.cells, gridSmoothing: undefined, gasInitialization: undefined,
    startup: 'partial-inviscid', converged: false, reason: 'iteration limit' });
  assert.equal(stages.find(s => s.stage === 'boundary-layer-initialization').eulerPreparation.converged, false);
  assert.equal(checkpoints[0].result.status, 'unconverged');
  const init = d.initializationCalls[0], n = createStreamtubeBodySystem(precursor.solverInput).layout.densityCount;
  assert.deepEqual(plain(init.options.initialEuler.x.slice(0, n)), plain(precursor.flow.x.slice(0, n)));
  assert.equal(init.options.ncrit, 4);
  assert.equal(init.options.transitionMode, 'automatic');
  assert.deepEqual(init.options.tripFractions, [[1, 1], [1, 1]]);
  assert.equal(d.coupledCalls[0].options.maxIterations, 1);
  assert.equal(d.coupledCalls[0].options.tolerance, controls.tolerance);
  assert.equal(d.coupledCalls[0].options.blUpdate, 'xfoil');
  assert.equal(d.coupledCalls[0].options.projectionGeometry, 'boundary-increment');
  assert.deepEqual(plain(precursor), original);
});

test('a partial first-order precursor preserves its gas state and restores the requested viscous order', async () => {
  const { input, precursor } = fixture();
  const data = { ...precursor.solverInput, streamwiseMode: 'momentum',
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const flow = solveStreamtubeIses(data, { maxIterations: 0, iterationGeometry: 'ises-sampled',
    stepAcceptance: 'armijo', firstOrderStartup: true, adaptiveMcrit: true });
  assert.equal(flow.solverInput.upwind.mucon, -1);
  flow.history.push({ iteration: 20, residual: flow.diagnostics.residual });
  precursor.flow = flow; precursor.solverInput = flow.solverInput;
  precursor.diagnostics.equationResidual = flow.diagnostics.residual;
  const before = plain(precursor), stop = new Error('inspect requested-order handoff');
  const driver = await assembly(precursor, { initialize(data, options) {
    assert.equal(data.upwind.mcrit, .99);
    assert.equal(data.upwind.mucon, 1);
    const source = createStreamtubeBodySystem(precursor.solverInput);
    assert.deepEqual(plain(options.initialEuler.x.slice(0, source.layout.densityCount)),
      plain(precursor.flow.x.slice(0, source.layout.densityCount)));
    assert.deepEqual(options.initialEuler.nodes, precursor.flow.nodes);
    throw stop;
  } });
  assert.throws(() => driver.solve(input, controls), error => error === stop);
  assert.deepEqual(plain(precursor), before);
});

test('the fully converged precursor path keeps its existing result and stage schemas', async () => {
  const { input, precursor } = fixture();
  const partialDriver = await assembly(structuredClone(precursor));
  const partialResult = plain(partialDriver.solve(input, controls));
  precursor.status = 'research-converged'; precursor.flow.converged = true; precursor.flow.residualConverged = true;
  precursor.flow.reason = 'residual'; precursor.diagnostics.reason = 'residual'; precursor.diagnostics.solverStopReason = 'residual';
  const d = await assembly(precursor), stages = [];
  const result = d.solve(input, { ...controls, onStage: value => stages.push(value) });
  assert.deepEqual(Object.keys(result.initialization.euler), ['iterations', 'residual', 'cells', 'gridSmoothing', 'gasInitialization']);
  assert(stages.every(stage => !('eulerPreparation' in stage)));
  assert.equal(d.initializationCalls.length, 1);
  for (const key of ['startup', 'converged', 'reason']) delete partialResult.initialization.euler[key];
  assert.deepEqual(plain(result), partialResult, 'Handoff classification must not change the physical initializer, driver or result fields.');
});

test('partial finite-base transfer preserves the prescribed gap and Euler variables before BL initialization', async () => {
  const { input, precursor } = fixture(true), stop = new Error('Stop after finite-base transfer.');
  const d = await assembly(precursor, { initialize(data, options) {
    const base = createStreamtubeBodySystem(data);
    const system = createStreamtubeBodySystem({ ...data,
      displacement: initialStreamtubeDisplacement(base.layout, base.baseGeometry) });
    const x = system.adoptGeometry(options.initialEuler.x, options.initialEuler.nodes), value = system.evaluate(x);
    assert.equal(system.layout.independentWakeBanks, true);
    assert(system.baseGeometry[0].width > 0);
    assert(system.displacement.wakes[0].some(d => d > 0));
    const n = system.layout.densityCount;
    assert.deepEqual(plain(x.slice(0, n)), plain(precursor.flow.x.slice(0, n)));
    const source = createStreamtubeBodySystem(precursor.solverInput);
    for (const name of Object.keys(source.layout.globals)) {
      const from = source.layout.globals[name], to = system.layout.globals[name];
      if (Array.isArray(from)) from.forEach((column, k) => {
        if (column !== null) assert.equal(x[to[k]], precursor.flow.x[column]);
      });
      else assert.equal(x[to], precursor.flow.x[from]);
    }
    assert(value.nodes.every(group => group.every(row => row.every(p => Number.isFinite(p.x) && Number.isFinite(p.y)))));
    throw stop;
  } });
  assert.throws(() => d.solve(input, controls), error => error === stop && error.stage === 'boundary-layer-initialization');
  assert.equal(d.initializationCalls.length, 1);
  assert.equal(d.coupledCalls.length, 0);
});

test('partial startup publishes and transfers the earlier complete seed while retaining the terminal Euler evidence', async () => {
  const { input, precursor } = fixture();
  const early = solveStreamtubeIses(precursor.solverInput, { maxIterations: 0,
    iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', retainBestCheckpoint: true });
  precursor.flow.bestCheckpoint = early.bestCheckpoint;
  precursor.flow.diagnostics.residual = 2 * early.diagnostics.residual;
  precursor.flow.finalQuality.minCornerSine = .5 * early.finalQuality.minCornerSine;
  precursor.diagnostics.equationResidual = precursor.flow.diagnostics.residual;
  const selected = selectCoupledEulerStartup(precursor), checkpoints = [], meshes = [];
  assert.ok(selected);
  const before = plain(precursor), d = await assembly(precursor);
  const result = d.solve(input, { ...controls,
    onCheckpoint: checkpoint => checkpoints.push(checkpoint), onMesh: mesh => meshes.push(mesh) });
  assert.equal(d.eulerCalls[0].options.retainBestCheckpoint, true);
  assert.deepEqual(plain(precursor), before);
  assert.deepEqual(plain(checkpoints[0].result), before);
  assert.equal(checkpoints[1].stage, 'euler-startup-selection');
  assert.equal(result.initialization.euler.iterations, 20);
  assert.equal(result.initialization.euler.selectedState.selectedIteration, 0);
  assert.equal(result.initialization.euler.selectedState.selectedResidual, early.diagnostics.residual);
  assert.equal(meshes[0].initialization.eulerStartup.selectedIteration, 0);
  const columns = createStreamtubeBodySystem(precursor.solverInput).layout.densityCount;
  assert.deepEqual(plain(d.initializationCalls[0].options.initialEuler.x.slice(0, columns)), plain(early.x.slice(0, columns)));
  assert.equal(result.converged, false);
});

test('rejected, invalid, unstarted or malformed precursor envelopes never trigger partial BL startup', async () => {
  const mutations = [
    p => { p.flow.lastRejectedStep = { stage: 'admissibility' }; },
    p => { p.flow.reason = 'Rejected Newton step'; },
    p => { p.diagnostics.reason = 'Final grid is invalid'; },
    p => { p.mesh.quality.valid = false; },
    p => { p.flow.finalQuality.valid = false; },
    p => { p.mesh.quality.invalidCells = [0]; },
    p => { p.flow.initialRedistribution.accepted = false; },
    p => { p.diagnostics.iterations = 0; p.flow.history = p.flow.history.slice(0, 1); },
    p => { p.flow.x[0] = NaN; },
    p => { p.flow.residual[0] = Infinity; },
    p => { p.diagnostics.equationResidual = NaN; },
    p => { p.flow.residual = []; },
  ];
  for (const mutate of mutations) {
    const { input, precursor } = fixture(); mutate(precursor);
    const d = await assembly(precursor);
    assert.throws(() => d.solve(input, controls), error => error.code === 'coupled-euler-precursor');
    assert.equal(d.initializationCalls.length, 0); assert.equal(d.coupledCalls.length, 0);
  }
  const { input, precursor } = fixture(), d = await assembly(precursor);
  assert.throws(() => d.solve(input, { ...controls, maxIterations: 0 }), error => error.code === 'coupled-euler-precursor');
  assert.equal(d.initializationCalls.length, 0);
});

test('partial handoff rechecks actual gas and geometry instead of trusting valid metadata', async () => {
  for (const kind of ['gas', 'grid', 'nodes']) {
    const { input, precursor } = fixture();
    if (kind === 'gas') precursor.flow.x[0] = 1000;
    else if (kind === 'nodes') precursor.flow.nodes = [];
    else {
      const before = precursor.flow.nodes, proposed = structuredClone(before); proposed[0][1][1].y += 1;
      const limit = limitStreamtubeGridStep(before, proposed);
      precursor.flow.nodes = interpolateStreamtubeGridNodes(before, proposed, limit.limiter.boundaryStep * 1.001);
    }
    const d = await assembly(precursor);
    assert.throws(() => d.solve(input, controls), error => error.stage === 'euler');
    assert.equal(d.initializationCalls.length, 0); assert.equal(d.coupledCalls.length, 0);
  }
});

test('cancellation and ordinary BL initialization failures remain failures rather than startup recovery', async () => {
  const { input, precursor } = fixture();
  for (const site of ['euler', 'checkpoint', 'stage', 'bl']) {
    const stop = new Error(`Canceled at ${site}`), d = await assembly(precursor, {
      ...(site === 'euler' ? { failEuler: stop } : {}), ...(site === 'bl' ? { initialize() { throw stop; } } : {}) });
    assert.throws(() => d.solve(input, { ...controls,
      onCheckpoint: () => { if (site === 'checkpoint') throw stop; },
      onStage: value => { if (site === 'stage' && value.stage === 'boundary-layer-initialization') throw stop; },
    }), error => error === stop);
    assert.equal(d.initializationCalls.length, site === 'bl' ? 1 : 0);
    assert.equal(d.coupledCalls.length, 0);
  }
});
