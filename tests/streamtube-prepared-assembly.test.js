import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { planCoupledLogarithmicShearRecovery } from '../src/euler/streamtube-coupled-log-shear-recovery.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// Geometry/initial-gas and controlled routing only: no Newton, Jacobian or LU.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { initializeStreamtubeStartup } from '../src/euler/streamtube-startup.js';
import { streamtubeEquationControls } from '../src/euler/streamtube-equation-selection.js';
import { capturePreparedStreamtubeAssembly, restorePreparedStreamtubeAssembly, retargetPreparedStreamtubeAssembly } from '../src/euler/streamtube-prepared-assembly.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { finiteBaseBodyFixture } from './fixtures/finite-base-body.js';
import { streamtubeWakeGap } from '../src/euler/streamtube-wake-geometry.js';

const copy = structuredClone;
function fixture({ elements = 1, ismom, finite = false, perturb = false } = {}) {
  const base = finite ? finiteBaseBodyFixture() : intrinsicBodyFixture({ elements, bodySegments: 4, tubes: 3 });
  const input = { ...base, flowModel: 'compressible', streamwiseMode: 'isentropic',
    normalStencil: 'body-stations', stagnationMotion: 'walls-only', geometryDomain: 'positive-simple',
    ...streamtubeEquationControls(ismom) };
  const caseData = { elements: input.bodies.map((b, i) => ({ name: `body${i}`, points: copy(b.points) })),
    mach: input.mach, alpha: input.alpha ?? 0, referenceChord: 1, gridIntervals: 4, gridTubes: 3,
    ...(ismom === undefined ? {} : { eulerIsmom: ismom }) };
  const system = createStreamtubeBodySystem(input), initial = system.initial.slice();
  if (perturb) {
    for (let i = 0; i < system.layout.densityCount; i++) initial[i] = .002 * Math.sin(i);
    for (let i = system.layout.globalOffset; i < initial.length; i++) initial[i] = .00001 * Math.cos(i);
    for (const { column } of system.layout.positions) initial[column] = .00001 * Math.sin(column);
  }
  const nodes = system.decode(initial).nodes, diagnostics = { gridSmoothing: { converged: true, history: [{ iteration: 0 }] } };
  const mesh = streamtubeMeshSnapshot({ system, initial, nodes, diagnostics });
  assert.equal(mesh.quality.valid, true);
  return { caseData, prepared: { input, system, initial, nodes, initialEuler: { x: initial.slice(), nodes: copy(nodes) },
    mesh, diagnostics, mach: input.mach, referenceChord: 1, momentReference: { x: .25, y: 0 } } };
}
const snapshot = p => ({ data: copy(Object.fromEntries(Object.entries(p).filter(([k]) => k !== 'system'))),
  conditions: copy(p.system.conditions), chart: p.system.geometryChart(), decoded: p.system.decode(p.initial) });

test('fresh Mach geometry reuse preserves a nonuniform two-body seed, masses/globals and caller source', () => {
  const { prepared, caseData } = fixture({ elements: 2, perturb: true }), before = snapshot(prepared), c0 = copy(caseData);
  const packet = capturePreparedStreamtubeAssembly(prepared, caseData);
  const nextCase = { ...copy(caseData), mach: .1 }, next = retargetPreparedStreamtubeAssembly(packet, nextCase).prepared;
  assert.notEqual(next.system, prepared.system);
  const old = before.decoded, value = next.system.decode(next.initial);
  for (const k of ['captured', 'allocation', 'stagnation']) assert.deepEqual(value[k], old[k]);
  assert.deepEqual(next.input, { ...before.data.input, mach: .1 });
  assert.deepEqual(next.initial.slice(0, next.system.layout.densityCount), prepared.initial.slice(0, prepared.system.layout.densityCount));
  assert.deepEqual(next.initial.slice(next.system.layout.globalOffset), prepared.initial.slice(prepared.system.layout.globalOffset));
  assert.ok(next.diagnostics.preparedEulerReuse.maximumNodeDeparture < 2e-14);
  assert.equal(next.diagnostics.preparedEulerReuse.geometryPreparationMach, .2);
  assert.equal(next.mesh.quality.valid, true); assert.equal(next.mesh.flow, undefined);
  const gas = next.system.conditions;
  assert.equal(gas.pInf, 1 / (gas.gamma * .1 * .1));
  assert.equal(gas.h0, 1 / ((gas.gamma - 1) * .1 * .1) + .5);
  assert.notEqual(gas.pInf, prepared.system.conditions.pInf);
  assert.deepEqual(snapshot(prepared), before); assert.deepEqual(caseData, c0);
});

test('retargeted default gas passes the real unchanged startup and retains mesh preparation provenance', () => {
  const { prepared, caseData } = fixture(), packet = capturePreparedStreamtubeAssembly(prepared, caseData);
  const next = retargetPreparedStreamtubeAssembly(packet, { ...caseData, mach: .1 });
  const seed = initializeStreamtubeStartup(next.prepared.system, next.prepared.initial);
  assert.equal(seed.diagnostics.method, 'isentropic');
  assert.ok(seed.flow.sections.flat(2).every(s => s.p > 0 && s.enthalpy > 0 && s.machSquared < 1));
  assert.deepEqual(seed.flow.nodes, next.prepared.system.decode(next.prepared.initial).nodes);
  const later = retargetPreparedStreamtubeAssembly(next, { ...caseData, mach: .05 });
  assert.equal(later.prepared.diagnostics.preparedEulerReuse.geometryPreparationMach, .2);
  assert.equal(later.prepared.diagnostics.preparedEulerReuse.sourceMach, .1);
});

test('each explicit equation selection remains unchanged while its gas Mach is retargeted', () => {
  for (const ismom of [1, 2, 3, 4]) {
    const { prepared, caseData } = fixture({ ismom });
    const result = retargetPreparedStreamtubeAssembly(capturePreparedStreamtubeAssembly(prepared, caseData), { ...caseData, mach: .1 });
    for (const [key, value] of Object.entries(streamtubeEquationControls(ismom))) {
      assert.deepEqual(result.prepared.input[key], value); assert.deepEqual(result.prepared.system.conditions[key], value);
    }
    assert.equal(result.prepared.diagnostics.preparedEulerReuse.equationsChanged, false);
  }
});

test('finite-base inviscid reuse preserves both physical corners and every constant wake gap', () => {
  const { prepared, caseData } = fixture({ finite: true });
  const next = retargetPreparedStreamtubeAssembly(capturePreparedStreamtubeAssembly(prepared, caseData), { ...caseData, mach: .1 }).prepared;
  assert.deepEqual(next.system.displacement, prepared.system.displacement);
  const nodes = next.nodes, te = next.input.bodies[0].trailingIndex;
  for (const [g, j] of [[0, nodes[0][te].length - 1], [1, 0]]) assert.deepEqual(nodes[g][te][j], prepared.nodes[g][te][j]);
  for (let i = te + 1; i <= next.system.layout.nx; i++) {
    const ids = [i - 1, i, Math.min(i + 1, next.system.layout.nx)];
    const gap = streamtubeWakeGap(ids.map(k => nodes[0][k].at(-1)), ids.map(k => nodes[1][k][0])).gap;
    assert.ok(Math.abs(gap - prepared.system.baseGeometry[0].width) < 2e-14);
  }
});

test('packet binding rejects changed controls/data/chart and returned metadata is independently owned', () => {
  const { prepared, caseData } = fixture(), packet = capturePreparedStreamtubeAssembly(prepared, caseData);
  assert.throws(() => restorePreparedStreamtubeAssembly({ ...packet }, caseData), /unrecognized/);
  for (const change of [{ alpha: 1 }, { gridTubes: 4 }, { eulerIsmom: 4 }, { referenceChord: 2 }, { elements: [] }])
    assert.throws(() => retargetPreparedStreamtubeAssembly(packet, { ...caseData, ...change, mach: .1 }), /only Mach/);
  for (const mach of [0, 1, NaN, undefined]) assert.throws(() => retargetPreparedStreamtubeAssembly(packet, { ...caseData, mach }), /only Mach/);
  assert.throws(() => restorePreparedStreamtubeAssembly(packet, { ...caseData, mach: .1 }), /retarget Mach explicitly/);
  const detached = restorePreparedStreamtubeAssembly(packet, caseData);
  detached.mesh.initialization.gasInitialization = { wrong: true }; detached.diagnostics.gridSmoothing.converged = false;
  assert.equal(restorePreparedStreamtubeAssembly(packet, caseData).mesh.initialization.gasInitialization, undefined);
  packet.prepared.initial[0] += .1;
  assert.throws(() => restorePreparedStreamtubeAssembly(packet, caseData), /saved data changed/);
  const fresh = capturePreparedStreamtubeAssembly(prepared, caseData);
  prepared.system.conditions.mach = .1;
  assert.throws(() => retargetPreparedStreamtubeAssembly(fresh, { ...caseData, mach: .1 }), /source system changed/);
});

test('capture rejects mismatched gas Mach, stale geometry and a completed or provisional flow', () => {
  for (const corrupt of [p => { p.input.mach = .1; }, p => { p.nodes[0][1][1].x += .01; },
    p => { p.mesh.flow = {}; }, p => { p.initialEuler.x[0] = .1; }]) {
    const { prepared, caseData } = fixture(); corrupt(prepared);
    assert.throws(() => capturePreparedStreamtubeAssembly(prepared, caseData), /Mach disagree|physical geometry|old flow|stale/);
  }
});

let sequence = 0;
async function routingHarness(failure) {
  const f = fixture(), calls = [], key = `__preparedRoute${++sequence}`;
  const stubs = { streamtubeEquationControls, capturePreparedStreamtubeAssembly, restorePreparedStreamtubeAssembly,
    createInitialStreamtubeTopology() { calls.push('topology'); return f.prepared.input; },
    createPanelStreamtubeGrid() { calls.push('panel'); return f.prepared; },
    prepareStreamtubeMesh(value) { calls.push('mesh'); return value; },
    refineStreamtubeBody() { calls.push('refine'); return f.prepared; }, streamtubeMeshSnapshot,
    initializeStreamtubeStartup() { calls.push('gas'); throw failure; } };
  globalThis[key] = stubs;
  const source = fs.readFileSync(new URL('../src/euler/streamtube-result.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  return { ...f, calls, release: () => delete globalThis[key],
    module: await import('data:text/javascript;base64,' + Buffer.from(source + `\n//# sourceURL=${key}`).toString('base64')) };
}

test('capture precedes typed gas failure; injection repeats no mesh work and preserves diagnostic identity', async () => {
  const failure = Object.assign(new Error('Nonpositive interface pressure.'), { code: 'streamtube-interface-pressure', diagnostics: { i: 2, group: 1, tube: 0 } });
  const h = await routingHarness(failure), meshes = []; let packet;
  try {
    const check = e => e.code === failure.code && e.diagnostics === failure.diagnostics && e.stage === 'gas-initialization' && e.cause === failure;
    assert.throws(() => h.module.solveStreamtubeAssembly(h.caseData, { onEulerPrepared: value => { packet = value; } }), check);
    assert.ok(packet);
    const nextCase = { ...h.caseData, mach: .1 }, next = retargetPreparedStreamtubeAssembly(packet, nextCase);
    assert.throws(() => h.module.solveStreamtubeAssembly(nextCase, { preparedEuler: next, onMesh: (m, stage) => meshes.push(stage) }), check);
    assert.deepEqual(h.calls, ['topology', 'panel', 'mesh', 'refine', 'gas', 'gas']);
    assert.deepEqual(meshes, ['initial']);
    assert.equal(restorePreparedStreamtubeAssembly(packet, h.caseData).mach, .2);
  } finally { h.release(); }
});

test('prepared observer cancellation is unwrapped and capacity wording remains specific to startup', async () => {
  const failure = Object.assign(new Error('capacity'), { code: 'streamtube-sonic-capacity', diagnostics: { capacityRatio: 1.104 } });
  const h = await routingHarness(failure), cancellation = Object.freeze({ cancel: true });
  try {
    assert.throws(() => h.module.solveStreamtubeAssembly(h.caseData, { onEulerPrepared: () => { throw cancellation; } }), e => e === cancellation);
    assert.ok(!h.calls.includes('gas'));
    assert.throws(() => h.module.solveStreamtubeAssembly(h.caseData), e => e.code === failure.code
      && e.stage === 'gas-initialization' && e.diagnostics === failure.diagnostics
      && e.message.includes('starting guess') && !e.message.includes('shocks are not supported'));
  } finally { h.release(); }
});

test('prepared injection rejects grid-construction options before mesh publication or gas work', async () => {
  const h = await routingHarness(new Error('stop after capture')); let packet;
  try {
    assert.throws(() => h.module.solveStreamtubeAssembly(h.caseData, { onEulerPrepared: value => { packet = value; } }), /stop after capture/);
    const before = h.calls.slice();
    for (const grid of [{ wallSubdivisions: 2 }, { refinementInterpolation: 'linear' }, { gridRepair: false }, { interiorInitialization: 'linear' }])
      assert.throws(() => h.module.solveStreamtubeAssembly(h.caseData, { preparedEuler: packet, ...grid,
        onMesh: () => assert.fail('Rejected option must not publish a reused mesh') }), /grid-construction options/);
    assert.deepEqual(h.calls, before);
  } finally { h.release(); }
});

test('coupled adapter forwards prepared options only to its initial Euler call and preserves cancellation', async () => {
  const key = `__preparedCoupled${++sequence}`, packet = {}, observer = () => {}, stop = 'stop before coupled work'; let received;
  globalThis[key] = { ...shearPolicy, planCoupledLogarithmicShearRecovery, streamtubeEquationControls, solveStreamtubeAssembly(_input, options) { received = options; throw stop; } };
  const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-assembly.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  try {
    const module = await import('data:text/javascript;base64,' + Buffer.from(source + `\n//# sourceURL=${key}`).toString('base64'));
    assert.throws(() => module.solveCoupledStreamtubeAssembly(fixture().caseData, { preparedEuler: packet, onEulerPrepared: observer }), e => e === stop);
    assert.equal(received.preparedEuler, packet); assert.equal(received.onEulerPrepared, observer);
  } finally { delete globalThis[key]; }
});

test('prepared Euler handoff retains the optional panel samples and rejects their mutation', () => {
 const {prepared,caseData}=fixture();
 prepared.panelVelocitySeed={method:'uncorrected-panel-velocity',speeds:Array(prepared.system.layout.densityCount).fill(1),maximumDirectionMismatch:0};
 const packet=capturePreparedStreamtubeAssembly(prepared,caseData);
 const restored=restorePreparedStreamtubeAssembly(packet,caseData);
 assert.deepEqual(restored.panelVelocitySeed,prepared.panelVelocitySeed);
 restored.panelVelocitySeed.speeds[0]=2;
 assert.equal(restorePreparedStreamtubeAssembly(packet,caseData).panelVelocitySeed.speeds[0],1);
 packet.prepared.panelVelocitySeed.speeds[0]=2;
 assert.throws(()=>restorePreparedStreamtubeAssembly(packet,caseData),/changed/);
});

test('harmonic shock preparation preserves final-grid boundaries and its reusable unsolved state', async () => {
  const { buildReliabilityCase } = await import('../scripts/validation/solver-reliability-cases.js');
  const { prepareStreamtubeAssembly } = await import('../src/euler/streamtube-result.js');
  const { caseData } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-grid',
    changes: { mach: .74, alpha: 2.68, gridIntervals: 8, gridTubes: 7 } });
  const standard = prepareStreamtubeAssembly(caseData);
  const selected = { ...caseData, eulerStartup: 'harmonic-shock' }, before = copy(selected);
  const harmonic = prepareStreamtubeAssembly(selected);
  assert.equal(harmonic.mesh.quality.valid, true);
  assert.equal(harmonic.mesh.initialization.finalGridSmoothing.converged, true);
  assert.equal(harmonic.mesh.cells.length, standard.mesh.cells.length);
  harmonic.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    if (!i || i === group.length - 1 || !j || j === row.length - 1)
      assert.deepEqual(p, standard.nodes[g][i][j]);
  })));
  const packet = capturePreparedStreamtubeAssembly(harmonic, selected);
  assert.deepEqual(restorePreparedStreamtubeAssembly(packet, selected).nodes, harmonic.nodes);
  assert.deepEqual(selected, before);
  const coupled = prepareStreamtubeAssembly({ ...selected, quadBoundaryLayers: true });
  assert.deepEqual(coupled.nodes, harmonic.nodes);
  assert.equal(coupled.mesh.quality.valid, true);
  assert.throws(() => prepareStreamtubeAssembly({ ...selected, eulerStartup: 'unknown' }), /Unknown Euler startup/);
});
