import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { coupledCheckpointHkPolicy } from '../src/euler/streamtube-coupled-assembly.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// Public routing only. All solve/preparation seams below use explicit stubs;
// these tests neither manufacture an accepted physical root nor run Newton.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const encode = text => 'data:text/javascript;base64,' + Buffer.from(text).toString('base64');
const archive = 'docs/coupled-projection-grid-response/public-integration/before/';
const read = file => fs.readFileSync(file, 'utf8');
let serial = 0;
// Reuse only existing controlled-fixture definitions, never their tests.
// Replacing import.meta.url preserves their original source-relative paths.
async function definitions(file, names, transform = text => text) {
  const url = pathToFileURL(path.resolve(file));
  let text = read(file).split('\ntest(')[0];
  assert(text.length < read(file).length, 'Fixture/test boundary must be present.');
  text = transform(text).replaceAll('import.meta.url', JSON.stringify(url.href))
    .replace(/from '([^']+)'/g, (_, specifier) => `from '${specifier.startsWith('.') ? new URL(specifier, url).href : specifier}'`);
  return import(encode(text + `\nexport { ${names.join(',')} };`));
}
function replaceOnce(text, before, after) {
  assert.equal(text.split(before).length, 2, `Expected one controlled fixture anchor: ${before}`);
  return text.replace(before, after);
}
const cold = await definitions('tests/streamtube-coarse-startup-routing.test.js', ['harness', 'baseCase', 'withFreshWakeChart'], text => {
  text = replaceOnce(text, 'snapshotActive:()=>[2,3]', 'snapshotActive:()=>Array(input.bodies.length*2).fill(2)');
  text = replaceOnce(text, 'calls.solves.push({grid,maxIterations:options.maxIterations})', `calls.solves.push({grid,maxIterations:options.maxIterations,
    wakeGeometry:input.wakeGeometry,blUpdate:options.blUpdate,projectionGeometry:options.projectionGeometry,
    ...(options.hkProjectionRecovery===undefined?{}:{hkProjectionRecovery:options.hkProjectionRecovery}),
    ...(options.resume?{resume:clone(options.resume)}:{})})`);
  text = replaceOnce(text, "options:{transitionMode:options.transitionMode,tripFractions:options.tripFractions,",
    "options:options.resume?.restart.options??{transitionMode:options.transitionMode,tripFractions:options.tripFractions,");
  text = replaceOnce(text, 'continuation:{}', `continuation:{...clone(options.resume?.continuation??{}),
    ...(options.blUpdate===undefined?{}:{blUpdate:options.blUpdate}),
    ...(options.projectionGeometry===undefined?{}:{projectionGeometry:options.projectionGeometry})}`);
  text = replaceOnce(text, "reason:converged?'residual':'line search'", "blUpdate:checkpoint.continuation.blUpdate,reason:converged?'residual':config.iterationLimit?'iteration limit':'line search'");
  text = replaceOnce(text, 'mesh:m,x,checkpoint,families:checkpoint.families',
    'mesh:m,x,...(config.noCheckpoint?{}:{checkpoint}),projectionGeometry:checkpoint.continuation.projectionGeometry,families:checkpoint.families');
  return text;
});
const mach = await definitions('tests/streamtube-coupled-mach-regions.test.js', ['adapter', 'parentFixture'], text =>
  replaceOnce(text, 'const imports = hooks ? {', `const imports = hooks ? {
    ...(hooks.cold ? {'./streamtube-coupled-assembly.js': dataURL(\`export const solveCoupledStreamtubeAssembly=globalThis['\${key}'].cold;
      export {coupledCheckpointHkPolicy,coupledAssemblyConditions} from '${pathToFileURL(path.resolve('src/euler/streamtube-coupled-assembly.js')).href}';\`)} : {}),`));

async function stubbed(file, hooks) {
  const key = `__projectionRouting${serial++}`;
  globalThis[key] = { ...shearPolicy, coupledCheckpointHkPolicy, ...hooks };
  const text = read(file).replace(/import \{([^}]+)\} from '[^']+';/g,
    (_, names) => `const {${names}}=globalThis[${JSON.stringify(key)}];`);
  try { return await import(encode(text)); } finally { delete globalThis[key]; }
}

test('direct startup keeps the Newton interior grid after the BL thickness correction', async () => {
  const h = await cold.harness(), input = cold.baseCase(), before = structuredClone(input);
  try {
    h.solveCoupledStreamtubeAssembly(input, { direct: true, maxIterations: 3 });
    assert.equal(h.calls.solves.length, 1);
    assert.equal(h.calls.solves[0].blUpdate, 'xfoil');
    assert.equal(h.calls.solves[0].projectionGeometry, 'fixed');
    assert.equal(h.calls.solves[0].wakeGeometry, 'independent-banks');
    assert.equal(h.calls.sequences.length, 0);
    assert.deepEqual(input, before);
  } finally { h.release(); }
});

test('an explicit exact-Hk request disables automatic projection recovery', async () => {
  const h = await cold.harness(), input = { ...cold.baseCase(), coupledNativeHk: false };
  try {
    h.solveCoupledStreamtubeAssembly(input, { direct: true, maxIterations: 3 });
    assert.equal(h.calls.solves[0].hkProjectionRecovery, false);
  } finally { h.release(); }
});

test('new automatic two-body coarse and requested fine solves both enable independent-bank response', async () => {
  const h = await cold.harness(), input = cold.baseCase();
  input.elements.push(structuredClone(input.elements[0])); input.materialTrips.push([1, 1]);
  const before = structuredClone(input);
  try {
    const result = h.solveCoupledStreamtubeAssembly(input, { maxIterations: 3 });
    assert.deepEqual(h.calls.solves.map(({ grid, blUpdate, projectionGeometry, wakeGeometry }) =>
      [grid, blUpdate, projectionGeometry, wakeGeometry]),
    [[16, 'xfoil', 'boundary-increment', 'independent-banks'], [36, 'xfoil', 'boundary-increment', 'independent-banks']]);
    assert.equal(h.calls.sequences.length, 1); assert.equal(result.bodies.length, 2);
    assert.equal(result.solverSettings.projectionGeometry, result.checkpoint.continuation.projectionGeometry);
    assert.equal(result.solverSettings.blUpdate, 'xfoil'); assert.deepEqual(input, before);
  } finally { h.release(); }
});

test('thinner fresh retries enable response while continuation chunks inherit their checkpoint', async () => {
  for (const iterationLimit of [false, true]) {
    const h = await cold.harness({ fineFailed: true, iterationLimit }), input = cold.baseCase();
    try {
      const result = h.solveCoupledStreamtubeAssembly(input, {
        maxIterations: iterationLimit ? 25 : 3, maxStartupAttempts: iterationLimit ? 1 : 2, coarseStartup: false });
      assert.equal(h.calls.solves.length, 2);
      assert.equal(h.calls.solves[0].projectionGeometry, 'boundary-increment');
      if (iterationLimit) {
        assert.equal(h.calls.solves[1].projectionGeometry, undefined);
        assert.equal(h.calls.solves[1].resume.continuation.projectionGeometry, 'boundary-increment');
        assert.deepEqual(h.calls.solves.map(c => c.maxIterations), [20, 5]);
      } else assert.equal(h.calls.solves[1].projectionGeometry, 'boundary-increment');
      assert.equal(result.solverSettings.projectionGeometry, 'boundary-increment');
    } finally { h.release(); }
  }
});

test('fixed-trip fresh path retains archived calls, checkpoint and public schema exactly', async () => {
  const a = await cold.harness(), b = await cold.harness({}, read(archive + 'streamtube-coupled-assembly.js.txt'));
  const input = { ...cold.baseCase(), transitionMode: 'fixed-trip', materialTrips: [[.05, .05]] };
  try {
    const current = a.solveCoupledStreamtubeAssembly(input, { maxIterations: 3 });
    const previous = b.solveCoupledStreamtubeAssembly(input, { maxIterations: 3 });
    assert.deepEqual(current, cold.withFreshWakeChart(previous)); assert.deepEqual(a.calls, b.calls);
    assert.equal(Object.hasOwn(current.solverSettings, 'projectionGeometry'), false);
    assert.equal(Object.hasOwn(current.checkpoint.continuation, 'projectionGeometry'), false);
  } finally { a.release(); b.release(); }
});


test('early cold failure without a checkpoint retains driver policy through public Mach reporting', async () => {
  const h = await cold.harness({ fineFailed: true, noCheckpoint: true }), input = cold.baseCase();
  try {
    const coldResult = h.solveCoupledStreamtubeAssembly(input, { maxIterations: 3, maxStartupAttempts: 1, coarseStartup: false });
    assert.equal(coldResult.converged, false); assert.equal(coldResult.checkpoint, undefined);
    assert.equal(coldResult.projectionGeometry, 'boundary-increment');
    assert.equal(coldResult.solverSettings.projectionGeometry, 'boundary-increment');
    // Body coordinate/normalization metadata comes from the existing pure
    // public-plan fixture; no numerical state is validated in this stub.
    const parent = mach.parentFixture(2, false), source = { ...parent, converged: false,
      checkpoint: undefined, projectionGeometry: coldResult.projectionGeometry,
      solverSettings: { ...parent.solverSettings, ...coldResult.solverSettings },
      conditions: { mach: .2 }, reason: 'Initial redistribution failed (controlled)' };
    const module = await mach.adapter('src/euler/tests/streamtube-coupled-mach-assembly.js', {
      cold: () => structuredClone(source), source: () => { throw Error('No warm replay expected'); },
      automatic: () => { throw Error('No continuation after failed cold state'); },
    });
    const result = module.solveCoupledStreamtubeMach(parent.sourceCase, {});
    assert.equal(result.checkpoint, undefined); assert.equal(result.converged, false);
    assert.equal(result.solverSettings.projectionGeometry, 'boundary-increment');
    assert.equal(result.reason, source.reason);
  } finally { h.release(); }
});

test('transition recovery inherits source policy and leaves legacy absence to the low-level default', async () => {
  for (const policy of [undefined, 'fixed', 'boundary-increment']) {
    const families = { euler: .1, boundaryLayer: .2, edgeMatching: .3 }, calls = [];
    const input = { wakeGeometry: 'independent-banks', mach: .2 };
    const result = { checkpoint: { version: 1, restart: { input, options: { transitionMode: 'automatic' },
      initialEuler: { x: [0] }, initialBL: [0] }, continuation: policy === undefined ? {} : { projectionGeometry: policy } },
      families, mesh: { initialization: { gridSmoothing: { enabled: true } } } };
    const before = structuredClone(result), mapped = { input, options: { transitionMode: 'automatic' },
      initialEuler: { x: [1] }, initialBL: [2], system: { n: 2 }, diagnostics: { mapped: true } };
    const hooks = {
      createCoupledStreamtubeBody: () => ({ n: 1, initial: [0], evaluate: () => ({ families }) }),
      refineCoupledStreamtubeBody: () => mapped,
      solveCoupledStreamtubeIses: (input, options) => {
        calls.push({ input, options }); return { mesh: { initialization: {} }, checkpoint: result.checkpoint };
      }, streamtubeMeshSnapshot: () => { throw Error('No mesh callback expected'); },
    };
    const module = await stubbed('src/euler/streamtube-transition-recovery.js', hooks);
    const value = module.recoverCoupledTransition(result, { plan: { normalFactor: 1 }, maxIterations: 0 });
    assert.equal(calls.length, 1); assert.equal(calls[0].options.projectionGeometry, policy);
    assert.equal(Object.hasOwn(calls[0].options, 'projectionGeometry'), policy !== undefined);
    assert.equal(calls[0].options.blUpdate, 'xfoil'); assert.equal(calls[0].input, mapped.input);
    assert.equal(value.mesh.initialization.gridSmoothing.enabled, true); assert.deepEqual(result, before);
  }
});

function warmHooks(parent, returnedPolicy) {
  const calls = [];
  return { calls, hooks: {
    source: () => ({ system: { euler: { conditions: { lengthScale: parent.solverLength } }, bl: { transitionMode: 'automatic' } } }),
    automatic: (_, options) => {
      calls.push(structuredClone(options.initialCheckpoint));
      const checkpoint = structuredClone(options.initialCheckpoint);
      if (returnedPolicy !== null) {
        if (returnedPolicy === undefined) delete checkpoint.continuation.projectionGeometry;
        else checkpoint.continuation.projectionGeometry = returnedPolicy;
      }
      return { checkpoint, projectionGeometry: 'boundary-increment', solverInput: checkpoint.restart.input, families: checkpoint.families,
        conditions: { mach: checkpoint.restart.input.mach }, mesh: { quality: { valid: true } },
        converged: false, stateConverged: true, continuation: { boundaryLayerReinitializedOnWarmRestart: false } };
    },
  } };
}

test('public warm continuation preserves legacy policy and reports actual checkpoint policy over stale labels', async () => {
  for (const policy of [undefined, 'boundary-increment']) {
    const parent = mach.parentFixture(2, false);
    if (policy) parent.checkpoint.continuation.projectionGeometry = policy;
    parent.solverSettings.projectionGeometry = policy ? 'fixed' : 'boundary-increment';
    const before = structuredClone(parent), controlled = warmHooks(parent, null);
    const module = await mach.adapter('src/euler/tests/streamtube-coupled-mach-assembly.js', controlled.hooks);
    const result = module.solveCoupledStreamtubeMach({ ...parent.sourceCase, mach: .24 }, { parentResult: parent,
      maxStages: 1, maxWakeRecoveries: 0, onPrepared: packet =>
        assert.equal(packet.parentResult.solverSettings.projectionGeometry, policy ?? 'fixed') });
    assert.deepEqual(controlled.calls, [before.checkpoint]);
    assert.deepEqual(result.checkpoint, before.checkpoint);
    assert.equal(result.solverSettings.projectionGeometry, policy ?? 'fixed'); assert.deepEqual(parent, before);
  }
  const parent = mach.parentFixture(2, false), controlled = warmHooks(parent, 'boundary-increment');
  const module = await mach.adapter('src/euler/tests/streamtube-coupled-mach-assembly.js', controlled.hooks);
  const result = module.solveCoupledStreamtubeMach({ ...parent.sourceCase, mach: .24 }, { parentResult: parent });
  assert.equal(result.solverSettings.projectionGeometry, 'boundary-increment', 'Returned checkpoint, not parent label, controls output.');
});

test('ordinary warm checkpoint preserves archived output with current policy metadata', async () => {
  const parent = mach.parentFixture(2, false), a = warmHooks(parent, null), b = warmHooks(parent, null);
  const current = await mach.adapter('src/euler/tests/streamtube-coupled-mach-assembly.js', a.hooks);
  const previous = await mach.adapter(archive + 'streamtube-coupled-mach-assembly.js.txt', b.hooks);
  const request = { ...parent.sourceCase, mach: .24 }, options = { parentResult: parent };
  const value = current.solveCoupledStreamtubeMach(request, options);
  const expected = previous.solveCoupledStreamtubeMach(request, options);
  expected.solverSettings.shearCoordinate = 'linear'; expected.solverSettings.hkFloorLinearization = 'exact';
  expected.solverSettings.maxMachStep = .05; expected.solverSettings.maxFirstOrderRecoveries = 2;
  assert.deepEqual(value, expected); assert.deepEqual(a.calls, b.calls);
  assert.equal(Object.hasOwn(value.solverSettings, 'projectionGeometry'), false);
});

test('public refinement resumes inherited checkpoint and labels actual returned policy', async () => {
  // Exercise the unmodified public solve/result function with explicitly
  // injected preparation and ISES stubs; no numerical refiner is executed.
  const source = read('src/euler/streamtube-coupled-refinement-assembly.js');
  const functionText = source.slice(source.indexOf('export function solveCoupledStreamtubeRefinement('));
  assert(functionText.startsWith('export function'));
  for (const policy of [undefined, 'boundary-increment']) {
    const parent = { solverSettings: { projectionGeometry: policy ? 'fixed' : 'boundary-increment' }, initialization: { euler: {} } };
    const checkpoint = { continuation: { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible',
      ...(policy ? { projectionGeometry: policy } : {}) }, restart: { options: {}, input: { bodies: [] } } };
    const calls = [], key = `__refinementProjection${serial++}`;
    globalThis[key] = { ...shearPolicy, coupledCheckpointHkPolicy,
      require: (condition, message) => { if (!condition) throw Error(message); },
      prepareCoupledStreamtubeRefinement: () => ({ checkpoint, settings: { normalization: {}, materialTrips: [], elementOrder: [] },
        refinement: {}, plan: { tolerance: 1e-10 } }),
      solveCoupledStreamtubeIses: (_, options) => { calls.push(options); return { checkpoint, history: [{}], families: {},
        mesh: { quality: { valid: true } }, solverInput: checkpoint.restart.input }; },
      streamtubeMeshSnapshot: () => { throw Error('No mesh callback expected'); },
    };
    try {
      const module = await import(encode(`const {coupledResultShearCoordinate,coupledCheckpointHkPolicy,require,prepareCoupledStreamtubeRefinement,solveCoupledStreamtubeIses,streamtubeMeshSnapshot}=globalThis['${key}'];\n` + functionText));
      const before = structuredClone(parent), value = module.solveCoupledStreamtubeRefinement({}, parent, { maxIterations: 0 });
      assert.equal(calls[0].resume, checkpoint); assert.equal(Object.hasOwn(calls[0], 'projectionGeometry'), false);
      assert.equal(value.solverSettings.projectionGeometry, policy ?? 'fixed'); assert.deepEqual(parent, before);
    } finally { delete globalThis[key]; }
  }
});
