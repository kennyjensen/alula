import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { planCoupledLogarithmicShearRecovery } from '../src/euler/streamtube-coupled-log-shear-recovery.js';
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { streamtubeEquationControls } from '../src/euler/streamtube-equation-selection.js';
import { coupledStartupExtensionPlan } from '../src/euler/streamtube-coupled-startup.js';

// Exercise the real finalizer/driver routing with controlled numerical
// dependencies. No Euler, BL, geometry, Jacobian or linear solve is run.
const archive = '../docs/solver-reliability/nlr-user-current-regression/remaining-budget-draft/';
const read = path => fs.readFileSync(new URL(path, import.meta.url), 'utf8');
let sequence = 0;
async function withStubs(source, stubs) {
  const key = `__coupledStopReason${++sequence}`;
  globalThis[key] = { ...shearPolicy, planCoupledLogarithmicShearRecovery, ...stubs };
  const code = source.replace(/import \{([^}]+)\} from '[^']+';/g,
    (_, names) => `const {${names.replace(/\s+as\s+/g, ': ')}} = globalThis[${JSON.stringify(key)}];`);
  try { return await import(`data:text/javascript;base64,${Buffer.from(code + `\n//# sourceURL=${key}.js`).toString('base64')}`); }
  finally { delete globalThis[key]; }
}
const nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
const families = residual => ({ euler: residual, boundaryLayer: residual / 2, edgeMatching: residual / 4 });
const mesh = valid => ({ quality: { valid, concaveCells: valid ? 0 : 1 }, initialization: {}, cells: [[0, 1, 2, 3]] });

async function finalizer(source = read('../src/euler/streamtube-coupled.js')) {
  const calls = [];
  const module = await withStubs(source, {
    solveStreamtubeBody(_system, options) { calls.push(options); return { nodes, undisplacedNodes: nodes }; },
    streamtubeMeshSnapshot({ system }) { return mesh(system.finalGridValid); },
  });
  function system(residual, valid) {
    return { ne: 1, initial: Float64Array.of(1, .01, .02, .03, .9),
      conditions: { edgeMatching: 'section-velocity' }, initialization: {},
      euler: { finalGridValid: valid },
      bl: { transitionMode: 'automatic', stations: [{ id: 0 }], surfaces: [], wakes: [] },
      evaluate() { return { residual: Float64Array.of(residual), families: families(residual), layers: { states: [{}] } }; } };
  }
  return { result: module.coupledStreamtubeResult, system, calls };
}

test('final grid rejection preserves the numerical stop and never certifies a residual-only root', async () => {
  const h = await finalizer();
  for (const [residual, reason, expected] of [[1, 'iteration limit', 'iteration limit'],
    [0, 'iteration limit', 'residual'], [1, 'Coupled ISES update rejected during admissibility: domain',
      'Coupled ISES update rejected during admissibility: domain']]) {
    const s = h.system(residual, false), r = h.result(s, s.initial, { reason, tolerance: 1e-10 });
    assert.equal(r.reason, 'Invalid final displacement grid');
    assert.equal(r.solverStopReason, expected);
    assert.equal(r.converged, false);
    assert.equal(r.mesh.initialization.flowSolved, false);
    assert.equal(r.mesh.quality.valid, false);
  }
  assert.ok(h.calls.every(c => c.maxIterations === 0));
});

test('valid/converged final result shapes and values remain exactly equal to the archived implementation', async () => {
  const current = await finalizer(), previous = await finalizer(read(archive + 'streamtube-coupled.js.before.txt'));
  for (const [residual, reason] of [[0, 'iteration limit'], [1, 'iteration limit'], [1, 'update rejected']]) {
    const a = current.system(residual, true), b = previous.system(residual, true);
    const r = current.result(a, a.initial, { reason, history: [{ iteration: 0 }] });
    assert.deepEqual(r, previous.result(b, b.initial, { reason, history: [{ iteration: 0 }] }));
    assert.equal(Object.hasOwn(r, 'solverStopReason'), false);
  }
});

test('ISES resume rejects a nonconvex checkpoint before any global solve or coordinate maintenance', async () => {
  const h = await finalizer(), s = h.system(1, false), admissibility = [];
  s.euler.layout = { densityCount: 1 };
  s.bl.snapshotActive = () => [2];
  s.admissibleValue = (_state, options) => {
    admissibility.push(options.requireConvex);
    return options.requireConvex ? null : s.evaluate();
  };
  const evaluate = s.evaluate;
  s.evaluate = () => ({ ...evaluate(), outer: { nodes, undisplacedNodes: nodes, stagnation: [0] } });
  const cp = { version: 1, families: families(1),
    restart: { input: { bodies: [{ leadingIndex: 1 }] }, options: { transitionMode: 'automatic', transitionState: [2] },
      initialEuler: { x: s.initial.slice(0, 1), nodes, undisplacedNodes: nodes }, initialBL: s.initial.slice(1) },
    continuation: { fractions: [[0, 1]], lastRedistributedStagnation: [0], preferredOrdering: 'amd', pivotTolerance: .001,
      iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing', blUpdate: 'xfoil' } };
  const before = structuredClone(cp);
  const driver = await withStubs(read('../src/euler/streamtube-coupled-ises.js'), {
    createCoupledStreamtubeBody: () => s, coupledStreamtubeResult: h.result,
    solveSparseDirect() { throw new Error('No linear solve is permitted in this test.'); },
    redistributeStreamtubeTangentially() { throw new Error('Resume must not repeat initial SMOVE.'); },
  });
  assert.throws(() => driver.solveCoupledStreamtubeIses(undefined, { resume: cp, maxIterations: 0,
    iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible' }), /Inadmissible coupled ISES initial state/);
  assert.deepEqual(admissibility, [true]);
  assert.equal(h.calls.length, 0, 'A rejected initial checkpoint must not reach the flow finalizer.');
  assert.deepEqual(cp, before);
});

async function assembly(config = {}, source = read('../src/euler/streamtube-coupled-assembly.js')) {
  const f = await finalizer(), calls = { solves: [], initializers: 0 }, checkpoints = [];
  const body = { element: 0 };
  const input = { elements: [{ name: 'test', points: [{ x: 1, y: 0 }, { x: 0, y: 0 }] }],
    gridIntervals: 16, gridTubes: 7, mach: .185, reynolds: 2.51e6, ncrit: 9, transitionMode: 'automatic' };
  const system = f.system(1, true);
  system.bl = { ...system.bl, hasFiniteBase: false, trips: [[1, 1]], snapshotActive: () => [2] };
  system.evaluate = () => ({ families: families(1), outer: { nodes, undisplacedNodes: nodes } });
  const stubs = {
    streamtubeEquationControls, coupledStartupExtensionPlan,
    solveStreamtubeAssembly(c) { return { status: 'research-converged', mach: c.mach,
      solverInput: { mach: c.mach, bodies: [body], streamwiseMode: 'isentropic' },
      flow: { x: Float64Array.of(1), nodes }, mesh: mesh(true), diagnostics: { iterations: 1, equationResidual: 1e-12, cells: 1 } }; },
    createStreamtubeBodySystem() { return { layout: { bodies: [body] }, conditions: { lengthScale: 1 },
      adoptGeometry: x => x.slice(), decode: () => ({ nodes, undisplacedNodes: nodes }) }; },
    initialStreamtubeDisplacement: () => ({}), transferStreamtubeGeometry: (_s, x) => x.slice(),
    initializeCoupledStreamtubeBody() { calls.initializers++; return { system, mesh: mesh(true), initialization: { thicknessFactor: 1 } }; },
    transitionRecoveryPlan: () => null,
    solveCoupledStreamtubeIses(solverInput, options) {
      const index = calls.solves.length;
      calls.solves.push(options);
      const valid = index ? !(config.secondInvalid ?? false) : (config.firstValid ?? false);
      const residual = index ? (config.secondInvalid ? 1 : 0) : (config.firstResidual ?? 1);
      const current = f.system(residual, valid);
      const iterations = options.maxIterations;
      const history = Array.from({ length: iterations + 1 }, (_, iteration) => ({ iteration, residual, step: iteration ? 1 : 0 }));
      const checkpoint = { version: 1, restart: { input: solverInput ?? options.resume.restart.input,
        options: { transitionMode: 'automatic', transitionState: [2] },
        initialEuler: { x: current.initial.slice(0, 1), nodes, undisplacedNodes: nodes }, initialBL: current.initial.slice(1) },
        continuation: { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', blUpdate: 'xfoil' }, families: families(residual) };
      checkpoints.push(checkpoint);
      const r = f.result(current, current.initial, { reason: index ? 'iteration limit' : (config.reason ?? 'iteration limit'),
        history, solverInput: checkpoint.restart.input, blUpdate: 'xfoil',
        ...(config.noCheckpoint ? {} : { checkpoint }),
        linearDiagnostics: { solves: iterations, refinements: 0, pivotRecoveries: 0, maxRelativeResidual: 0,
          iterations: history.slice(1).map(({ iteration }) => ({ iteration })) } });
      history.forEach(h => options.onIteration?.(h));
      return r;
    },
  };
  const module = await withStubs(source, stubs);
  return { solve: options => module.solveCoupledStreamtubeAssembly(input, { maxIterations: 40, maxStartupAttempts: 1, ...options }),
    calls, checkpoints, input };
}

test('an iteration-limited nonconvex endpoint resumes the same checkpoint for only the remaining 20 of 40 updates', async () => {
  const h = await assembly(), before = structuredClone(h.input), observed = [];
  const r = h.solve({ onIteration: h => { if (h.stage === 'coupled') observed.push(h.iteration); } });
  assert.deepEqual(h.calls.solves.map(c => c.maxIterations), [20, 20]);
  assert.equal(h.calls.solves[1].resume, h.checkpoints[0]);
  assert.equal(h.calls.solves[1].initialBL, undefined);
  assert.equal(h.calls.solves[1].initialEuler, undefined);
  assert.equal(h.calls.solves[1].iterationGeometry, 'ises-sampled');
  assert.equal(h.calls.solves[1].stepAcceptance, 'admissible');
  assert.equal(h.calls.initializers, 1);
  assert.deepEqual(observed, Array.from({ length: 41 }, (_, i) => i));
  assert.deepEqual(r.history.map(h => h.iteration), observed);
  assert.equal(r.linearDiagnostics.solves, 40);
  assert.equal(r.initialization.attempts[0].iterations, 40);
  assert.equal(r.converged, true);
  assert.equal(Object.hasOwn(r, 'solverStopReason'), false);
  assert.equal(r.physicalAcceptance, false);
  assert.deepEqual(h.input, before);
});

test('a true rejected update, a residual-only invalid root, a missing checkpoint or no remaining budget is not resumed', async () => {
  for (const [config, options] of [[{ reason: 'Coupled ISES update rejected during material-trip transfer: domain' }, {}],
    [{ firstResidual: 0 }, {}], [{ noCheckpoint: true }, {}], [{}, { maxIterations: 20 }]]) {
    const h = await assembly(config), r = h.solve(options);
    assert.equal(h.calls.solves.length, 1);
    assert.equal(r.converged, false);
    assert.equal(r.reason, 'Invalid final displacement grid');
  }
});

test('using the remaining budget does not weaken the final convexity gate', async () => {
  const h = await assembly({ secondInvalid: true }), r = h.solve();
  assert.deepEqual(h.calls.solves.map(c => c.maxIterations), [20, 20]);
  assert.equal(r.converged, false);
  assert.equal(r.mesh.quality.valid, false);
  assert.equal(r.mesh.initialization.flowSolved, false);
  assert.equal(r.reason, 'Invalid final displacement grid');
  assert.equal(r.solverStopReason, 'iteration limit');
});

test('ordinary valid-grid remaining-budget routing stays exactly equal to the archived assembly', async () => {
  const a = await assembly({ firstValid: true });
  const b = await assembly({ firstValid: true }, read(archive + 'streamtube-coupled-assembly.js.before.txt'));
  const expected = b.solve();
  // Add the separately qualified fresh wake chart and explicit legacy linear/exact policy labels.
  const freshChart = value => {
    if (!value || typeof value !== 'object' || ArrayBuffer.isView(value)) return;
    if (value.solverSettings) {
      value.solverSettings.shearCoordinate ??= 'linear'; value.solverSettings.hkFloorLinearization ??= 'exact';
    }
    if (value.wakeGeometry === 'independent-banks' && value.wakeOutlet === 'banks' && !Object.hasOwn(value, 'stepMethod'))
      value.wakeDisplacementMotion = 'te-center';
    for (const child of Object.values(value)) freshChart(child);
  };
  freshChart(expected);
  assert.deepEqual(a.solve(), expected);
  assert.deepEqual(a.calls.solves.map(c => c.maxIterations), [20, 20]);
  assert.equal(a.calls.initializers, b.calls.initializers);
});
