// SPDX-License-Identifier: GPL-2.0-or-later
// Numerical dependencies are controlled stubs: these test adapter routing,
// publication and legacy compatibility, not convergence of any ISMOM model.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { streamtubeEquationControls } from '../src/euler/streamtube-equation-selection.js';

const current = fs.readFileSync(new URL('../src/euler/streamtube-result.js', import.meta.url), 'utf8');
const baseline = fs.readFileSync(new URL('../docs/euler-ismom/inviscid-before/streamtube-result.js.txt', import.meta.url), 'utf8');
let sequence = 0;
const clone = structuredClone;
async function harness(config = {}, source = current) {
  const calls = [], nodes = [[[{ x: 0, y: 0 }, { x: 0, y: 1 }], [{ x: 1, y: 0 }, { x: 1, y: 1 }]]];
  const mesh = (iteration = 0) => ({ vertices: clone(nodes.flat(2)), cells: [[0, 1, 3, 2]],
    quality: { valid: true }, initialization: {}, iteration: { iteration }, flow: { speed: [1] } });
  const conditions = input => ({ mach: input.mach, alpha: input.alpha ?? 0, flowModel: 'compressible',
    pInf: 1, streamwiseMode: input.streamwiseMode,
    ...(input.hybrid ? { hybrid: clone(input.hybrid), upwind: clone(input.upwind) } : {}) });
  const system = input => ({ conditions: conditions(input), layout: { tubes: [1, 1], elements: 1, n: 4 } });
  const flow = input => ({ nodes: clone(nodes), history: [{ iteration: 0 }, { iteration: 1 }],
    converged: !config.unconverged, residualConverged: !config.unconverged, reason: config.unconverged ? 'iteration limit' : 'residual',
    diagnostics: { residual: config.unconverged ? .1 : 1e-12 }, solverInput: clone(input),
    streamwiseMode: input.streamwiseMode, jacobianBackend: 'stub', linearBackend: 'stub',
    stagnationLimiter: 'stub', iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible' });
  const stubs = {
    streamtubeEquationControls,
    createInitialStreamtubeTopology(input) { calls.push(['topology', clone(input)]); return clone(input); },
    createPanelStreamtubeGrid(input) {
      calls.push(['grid', clone(input)]);
      return { input, system: system(input), initial: [0], nodes: clone(nodes), diagnostics: {}, mesh: mesh() };
    },
    prepareStreamtubeMesh(value) { calls.push(['mesh', clone(value.input)]); return value; },
    prepareCurvedStreamtubePreview() { throw new Error('Unexpected curved preview'); },
    refineStreamtubeBody(input) {
      calls.push(['refine', clone(input)]);
      return { input, system: system(input), initial: [0], initialEuler: { nodes: clone(nodes) }, diagnostics: {} };
    },
    streamtubeMeshSnapshot(value) { return mesh(value.iteration?.iteration ?? 0); },
    initializeStreamtubeStartup(value) {
      calls.push(['startup', clone(value.conditions)]);
      if (config.sonic) throw Object.assign(new Error('capacity'), { code: 'streamtube-sonic-capacity',
        diagnostics: { capacityRatio: 1.1, stagnationDensityFallback: { reason: 'domain rejected' } } });
      return { initial: [0], flow: { nodes: clone(nodes) }, diagnostics: { kind: 'stub', method: config.startupMethod } };
    },
    streamtubeFlowSnapshot() { return { speed: [1] }; },
    compareStreamtubeFlow(value) { return value; },
    solveStreamtubeIses(input, options) {
      calls.push(['solve', clone(input), { maxIterations: options.maxIterations, tolerance: options.tolerance,
        iterationGeometry: options.iterationGeometry, stepAcceptance: options.stepAcceptance,
        ...(options.firstOrderStartup === undefined ? {} : { firstOrderStartup: options.firstOrderStartup }),
        ...(options.stopOnGridStagnation === undefined ? {} : { stopOnGridStagnation: options.stopOnGridStagnation }),
        ...(options.retainCheckpoint === undefined ? {} : { retainCheckpoint: options.retainCheckpoint }) }]);
      options.onIteration?.({ iteration: 1, residual: 1e-12 });
      return flow(input);
    },
    streamtubeForceCoefficients(_flow, gas, options) {
      calls.push(['legacy-force', clone(gas), clone(options)]); return { cl: .1, cd: .2, cm: .3 };
    },
    streamtubeSolidPressureForces({ flow: value, conditions: gas, referenceChord, momentReference }) {
      calls.push(['solid-force', clone(gas), { referenceChord, momentReference }]);
      assert.equal(value.solverInput.hybrid.ismom, gas.hybrid.ismom);
      return { cl: .4, pressureIntegralDrag: .5, cm: .6, physicalAcceptance: false };
    },
  };
  const key = `__inviscidIsmom${++sequence}`; globalThis[key] = stubs;
  const rewritten = source.replace(/import \{([^}]+)\} from '[^']+';/g,
    (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  const loaded = await import('data:text/javascript;base64,' + Buffer.from(rewritten + `\n//# sourceURL=${key}.js`).toString('base64'));
  return { ...loaded, calls, release: () => delete globalThis[key] };
}
const input = () => ({ elements: [{ name: 'control', points: [{ x: 0, y: 0 }, { x: 1, y: 0 }] }],
  mach: .63, alpha: 2, referenceChord: 2, momentReference: { x: .5, y: .1 } });

test('subsonic multielement startup enables stagnation recovery without changing single-element routing', async () => {
  for (const mach of [.01, .2, .63]) for (const count of [1, 2]) {
    const h = await harness(), base = input();
    try {
      h.solveStreamtubeAssembly({ ...base, mach, eulerIsmom: 4,
        elements: Array.from({ length: count }, () => clone(base.elements[0])) });
      const [, solverInput, controls] = h.calls.find(c => c[0] === 'solve');
      assert.equal(solverInput.stagnationMotion, 'walls-only');
      assert.equal(controls.stopOnGridStagnation ?? false, count > 1 && mach <= .3);
    } finally { h.release(); }
  }
});

test('all four explicit inviscid choices retain one operator through preparation, solve and physical-pressure publication', async () => {
  for (const ismom of [1, 2, 3, 4]) {
    const h = await harness({ unconverged: ismom === 3 }), c = { ...input(), eulerIsmom: ismom }, before = clone(c), observations = [];
    try {
      const result = h.solveStreamtubeAssembly(c, { maxIterations: 2, onIteration: e => observations.push(e) });
      for (const kind of ['grid', 'refine', 'solve']) {
        const found = h.calls.filter(c => c[0] === kind); assert.equal(found.length, 1);
        for (const [key, value] of Object.entries(streamtubeEquationControls(ismom))) assert.deepEqual(found[0][1][key], value);
        assert.equal(found[0][1].mach, before.mach);
      }
      assert.equal(h.calls.find(c => c[0] === 'solve')[2].retainCheckpoint, true);
      assert.equal(h.calls.find(c => c[0] === 'solve')[2].stepAcceptance, 'armijo');
      assert.equal(h.calls.filter(c => c[0] === 'startup').length, 1);
      assert.equal(h.calls.filter(c => c[0] === 'legacy-force').length, 0);
      assert.equal(h.calls.filter(c => c[0] === 'solid-force').length, 1);
      assert.deepEqual([result.cl, result.cd, result.cm], [.4, .5, .6]);
      assert.equal(result.coefficients.physicalAcceptance, false);
      assert.equal(result.solverInput.hybrid.ismom, ismom);
      assert.deepEqual(result.solverSettings.hybrid, { epsilonP: 1e-5, ismom });
      assert.deepEqual(result.solverSettings.upwind, streamtubeEquationControls(ismom).upwind);
      assert.equal(result.status, ismom === 3 ? 'unconverged' : 'research-converged');
      assert.ok(!result.limitations.includes('shocks are not supported'));
      assert.equal(observations[0].mach, before.mach); assert.deepEqual(c, before);
    } finally { h.release(); }
  }
});

test('explicit research acceptance remains selectable for an ISMOM solve', async () => {
  const h = await harness();
  try {
    h.solveStreamtubeAssembly({ ...input(), eulerIsmom: 4 }, { stepAcceptance: 'admissible' });
    assert.equal(h.calls.find(c => c[0] === 'solve')[2].stepAcceptance, 'admissible');
  } finally { h.release(); }
});

test('temporary first-order startup is routed only for standalone fallback gas initialization', async () => {
  for (const quadBoundaryLayers of [false, true, undefined])
  for (const startupMethod of ['isentropic', 'stagnation-density']) {
    const h = await harness({ startupMethod });
    try {
      const result = h.solveStreamtubeAssembly({ ...input(), eulerIsmom: 4, quadBoundaryLayers });
      const enabled = quadBoundaryLayers === false && startupMethod === 'stagnation-density';
      assert.equal(h.calls.find(c => c[0] === 'solve')[2].firstOrderStartup ?? false, enabled);
      assert.equal(result.solverSettings.firstOrderStartup ?? false, enabled);
    } finally { h.release(); }
  }
  const h = await harness({ startupMethod: 'stagnation-density' });
  try {
    h.solveStreamtubeAssembly({ ...input(), eulerIsmom: 4, quadBoundaryLayers: false }, { firstOrderStartup: false });
    assert.equal(h.calls.find(c => c[0] === 'solve')[2].firstOrderStartup, undefined);
  } finally { h.release(); }
});

test('omitted ISMOM preserves complete archived adapter output, calls and observations', async () => {
  for (const unconverged of [false, true]) {
    const a = await harness({ unconverged }), b = await harness({ unconverged }, baseline), events = [[], []];
    try {
      const invoke = (h, i) => h.solveStreamtubeAssembly(input(), { onMesh: (m, s) => events[i].push([clone(m), s]),
        onIteration: e => events[i].push(clone(e)) });
      const actual = invoke(a, 0), expected = invoke(b, 1);
      // Startup accounting was added after this archived adapter. Its
      // metadata does not change the legacy solve or physical output.
      assert.equal(actual.diagnostics.totalIterations, expected.diagnostics.iterations);
      assert.deepEqual(actual.diagnostics.startupAttempts, [{ attempt: 1,
        iterations: expected.diagnostics.iterations, mucon: undefined,
        reason: expected.diagnostics.reason, residual: expected.diagnostics.residual }]);
      delete actual.diagnostics.totalIterations;
      delete actual.diagnostics.startupAttempts;
      assert.deepEqual(actual, expected); assert.deepEqual(a.calls, b.calls);
      for (const event of events[0]) if (!Array.isArray(event)) {
        assert.equal(event.startupAttempt, 1);
        delete event.startupAttempt;
      }
      assert.deepEqual(events[0], events[1]);
    } finally { a.release(); b.release(); }
  }
});

test('explicit startup rejection preserves selected mode and never starts a fallback solve', async () => {
  const h = await harness({ sonic: true });
  try {
    assert.throws(() => h.solveStreamtubeAssembly({ ...input(), eulerIsmom: 1 }), e =>
      e.code === 'streamtube-sonic-capacity' && /ISMOM 1/.test(e.message) && !/shocks are not supported/.test(e.message));
    assert.equal(h.calls.filter(c => c[0] === 'grid').length, 1);
    assert.equal(h.calls.filter(c => c[0] === 'solve').length, 0);
    assert.equal(h.calls.find(c => c[0] === 'startup')[1].hybrid.ismom, 1);
  } finally { h.release(); }
});

test('invalid ISMOM fails before mesh work; mesh-only explicit choice performs no gas or solve work', async () => {
  const h = await harness();
  try {
    assert.throws(() => h.solveStreamtubeAssembly({ ...input(), eulerIsmom: '4' }), /ISMOM/); assert.deepEqual(h.calls, []);
    const result = h.solveStreamtubeAssembly({ ...input(), eulerIsmom: 4 }, { meshOnly: true });
    assert.equal(result.status, 'mesh-ready');
    assert.ok(h.calls.every(c => !['startup', 'solve', 'solid-force', 'legacy-force'].includes(c[0])));
    assert.equal(h.calls.find(c => c[0] === 'grid')[1].hybrid.ismom, 4);
  } finally { h.release(); }
});
