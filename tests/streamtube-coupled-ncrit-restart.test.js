// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { fileURLToPath } from 'node:url';
import { readFile } from 'node:fs/promises';
import { createServer } from 'node:http';
import path from 'node:path';
import { chromium } from '@playwright/test';
import { prepareCoupledNcritRestart, CoupledNcritRestartError } from '../src/euler/streamtube-coupled-ncrit-restart.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { prepareCoupledMrchduProfiles } from '../src/euler/streamtube-coupled-mrchdu-predictor.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const copy = v => JSON.parse(JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x));
const tolerance = 1e-10;
let cached;
function sourceRoot() {
  if (!cached) {
    // One existing 349-unknown, two-body test fixture solve for the whole file.
    // Every preparation/replay below performs zero global solves.
    const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 }),
      streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
      upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
    const root = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity',
      blThermodynamics: 'historical-common-isentrope', transitionMode: 'automatic', tripFractions: [[1, 1], [1, 1]],
      hkFloorLinearization: 'native', ncrit: 4, maxIterations: 12, tolerance, stepAcceptance: 'admissible' });
    assert.equal(root.converged, true, root.reason); assert.ok(root.x.length < 500);
    assert.equal(root.mesh.quality.valid, true); cached = copy(root.checkpoint);
  }
  return structuredClone(cached);
}
function replay(checkpoint) {
  const r = checkpoint.restart;
  const system = createCoupledStreamtubeBody(r.input, { ...r.options, initialEuler: r.initialEuler, initialBL: r.initialBL });
  return { system, value: system.evaluate(system.initial) };
}
function physicalEqual(a, b) {
  for (const key of ['allocation', 'sections', 'captured', 'stagnation', 'strengths', 'nodes', 'undisplacedNodes', 'residual'])
    assert.deepEqual(a.outer[key], b.outer[key], key);
  a.layers.states.forEach((s, id) => {
    for (const key of ['theta', 'deltaStar', 'ue', 's', 'wakeGap']) assert.equal(s[key], b.layers.states[id][key], `${id}/${key}`);
  });
}

test('Ncrit continuation changes only native auxiliaries and phase on a real qualified two-body source', t => {
  const cp = sourceRoot(), before = structuredClone(cp);
  const prepared = prepareCoupledNcritRestart(5, cp, { requestedNcrit: 9 });
  const { checkpoint, seed, diagnostics: d } = prepared;
  assert.deepEqual(cp, before);
  assert.deepEqual(d.sourcePhase, [1, 1, 1, 1]); assert.deepEqual(d.targetPhase, [2, 2, 2, 2]);
  assert.ok(d.auxiliaryChanges.length > 0); assert.equal(d.targetPhaseInitialization.changed, true);
  assert.deepEqual(checkpoint.restart.input, cp.restart.input);
  assert.deepEqual(checkpoint.restart.initialEuler, cp.restart.initialEuler);
  assert.deepEqual(checkpoint.restart.options, { ...cp.restart.options, ncrit: 5, transitionState: [2, 2, 2, 2] });
  assert.deepEqual(checkpoint.continuation, cp.continuation);
  checkpoint.restart.initialBL.forEach((x, i) => { if (i % 4) assert.equal(x, cp.restart.initialBL[i]); });
  const a = replay(cp), b = replay(copy(checkpoint)); physicalEqual(a.value, b.value);
  assert.deepEqual(b.value.families, checkpoint.families); assert.equal(b.system.conditions.hkFloorLinearization, 'native');
  assert.equal(d.quality.valid, true); assert.equal(d.quality.cellsChecked, 102);
  assert.equal(d.operations.nativeProfileCalls, 0); assert.equal(d.operations.globalJacobians, 0);
  assert.equal(d.operations.globalLinearSolves, 0); assert.equal(d.operations.globalNewtonUpdates, 0);
  assert.equal(seed.converged, false); assert.equal(seed.physicalAcceptance, false); assert.equal(seed.intermediate, true);
  assert.equal(d.sourceConverged, true); assert.equal(d.targetConverged, false); assert.equal(d.provisional, true);
  assert.deepEqual(d.warnings, []); assert.equal(d.exactReplay, true);
  assert.ok(d.targetFamilies.boundaryLayer > tolerance);
  // The unchanged ISES driver can resume the actual serialized complete seed.
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = checkpoint.continuation;
  const zero = solveCoupledStreamtubeIses(undefined, { resume: copy(checkpoint), maxIterations: 0, tolerance,
    iterationGeometry, stepAcceptance, stagnationLimiter });
  assert.equal(zero.linearDiagnostics.solves, 0); assert.equal(zero.converged, false);
  assert.deepEqual(copy(zero.checkpoint), copy(checkpoint)); assert.deepEqual(zero.residual, b.value.residual);
  t.diagnostic(JSON.stringify({ sourceFamilies: cp.families, targetFamilies: d.targetFamilies,
    actualCells: d.quality.cellsChecked, convertedAuxiliaries: d.auxiliaryChanges.length }));
});

test('criterion changes within the same active interval preserve every stored BL auxiliary', () => {
  const cp = sourceRoot(), p = prepareCoupledNcritRestart(4.01, cp, { requestedNcrit: 9 });
  assert.deepEqual(p.diagnostics.targetPhase, p.diagnostics.sourcePhase);
  assert.equal(p.diagnostics.targetPhaseInitialization.changed, false);
  assert.deepEqual(p.diagnostics.auxiliaryChanges, []);
  assert.deepEqual(p.checkpoint.restart.initialBL, cp.restart.initialBL);
  assert.ok(p.diagnostics.targetFamilies.boundaryLayer > tolerance);
});

test('omitted and explicit preserve retain the entire archived pre-MRCHDU result', async () => {
  const archive = await readFile(new URL('../docs/solver-reliability/nlr32-ismom4-smoothed/ncrit-restart-native-draft/preserve.before.js.txt', import.meta.url), 'utf8');
  const previous = await import('data:text/javascript;base64,' + Buffer.from(archive.replace(/from '([^']+)'/g,
    (_, p) => `from '${new URL(p, new URL('../scripts/validation/coupled-ncrit-restart-draft.js', import.meta.url)).href}'`)).toString('base64'));
  const cp = sourceRoot();
  for (const criterion of [4.01, 5, 20]) {
    const options = { requestedNcrit: 20 }, old = previous.prepareCoupledNcritRestart(criterion, cp, options);
    assert.deepEqual(prepareCoupledNcritRestart(criterion, cp, options), old);
    assert.deepEqual(prepareCoupledNcritRestart(criterion, cp, { ...options, blPredictor: 'preserve' }), old);
  }
});

test('explicit native prediction retains warnings, physical masses and density while transferring only boundary increments', () => {
  const cp = sourceRoot(), before = structuredClone(cp), r = cp.restart;
  const prepared = prepareCoupledNcritRestart(5, cp, { requestedNcrit: 9, blPredictor: 'xfoil-mrchdu' });
  assert.deepEqual(cp, before);
  const a = replay(cp), b = replay(copy(prepared.checkpoint));
  const target = createCoupledStreamtubeBody(r.input, { ...r.options, ncrit: 5, initialEuler: r.initialEuler, initialBL: r.initialBL });
  const native = prepareCoupledMrchduProfiles({ bl: target.bl, states: a.value.layers.states, initialBL: r.initialBL, targetMach: r.input.mach });
  const d = prepared.diagnostics;
  assert.deepEqual(d.prediction, native.diagnostics);
  assert.equal(d.operations.nativeProfileCalls, 2); assert.equal(d.operations.globalLinearSolves, 0);
  assert.equal(d.operations.globalJacobians, 0); assert.equal(d.operations.globalNewtonUpdates, 0);
  assert.deepEqual(d.warnings, native.diagnostics.bodies.flatMap(body => body.localConvergenceWarnings));
  assert.equal(d.warnings.length, 2, 'The qualified tiny source exposes two real native local warnings at this target.');
  assert.equal(d.warningFree, false); assert.equal(d.prepared, true);
  assert.equal(d.targetConverged, false); assert.equal(prepared.seed.converged, false); assert.equal(prepared.seed.physicalAcceptance, false);
  assert.equal(d.physicalBLPreserved, false); assert.equal(d.physicalEulerPreserved, false);
  assert.equal(d.physicalDensityPreserved, true); assert.equal(d.physicalMassPreserved, true); assert.equal(d.physicalGlobalsPreserved, true);
  assert.deepEqual(prepared.checkpoint.restart.input, r.input);
  assert.deepEqual(prepared.checkpoint.continuation, cp.continuation);
  assert.deepEqual(b.system.bl.kernel.parameters, { ...a.system.bl.kernel.parameters, ncrit: 5 });
  for (const key of ['allocation', 'captured', 'stagnation', 'strengths']) assert.deepEqual(b.value.outer[key], a.value.outer[key]);
  assert.deepEqual(b.value.outer.sections.map(row => row.map(g => g.map(s => s.rho))),
    a.value.outer.sections.map(row => row.map(g => g.map(s => s.rho))));
  prepared.checkpoint.restart.initialBL.forEach((x, i) => {
    if (i % 4) assert.equal(x, native.initialBL[i], 'Only auxiliary phase conversion follows the native profile prediction.');
  });
  assert.notDeepEqual(b.value.outer.nodes, a.value.outer.nodes);
  b.value.outer.nodes.forEach((group, g) => {
    const masses = a.value.outer.allocation.groups[g].map(t => t.massFlow), total = masses.reduce((x, y) => x + y, 0);
    group.forEach((row, i) => {
      const old = a.value.outer.nodes[g][i], last = row.length - 1;
      let eta = 0;
      for (let j = 1; j < last; j++) {
        eta += masses[j - 1] / total;
        for (const axis of ['x', 'y']) {
          const expected = old[j][axis] + (1 - eta) * (row[0][axis] - old[0][axis]) + eta * (row[last][axis] - old[last][axis]);
          assert.ok(Math.abs(row[j][axis] - expected) <= 8 * Number.EPSILON * Math.max(1, Math.abs(expected)));
        }
      }
    });
  });
  assert.equal(d.quality.valid, true); assert.deepEqual(b.value.families, prepared.checkpoint.families);
  assert.throws(() => prepareCoupledNcritRestart(6, prepared.checkpoint, { blPredictor: 'xfoil-mrchdu' }), /source must already satisfy/);
});

test('terminal fallback uses the native auxiliary handoff and remains a provisional requested-condition guess', () => {
  const cp = sourceRoot(), p = prepareCoupledNcritRestart(20, cp);
  assert.ok(p.diagnostics.targetPhaseInitialization.targets.every(t => t.kind === 'trailing-edge' && t.fraction === 1));
  assert.deepEqual(p.diagnostics.targetPhase, [3, 3, 3, 3]);
  assert.equal(p.seed.selectedNcrit, 20); assert.equal(p.seed.requestedNcrit, 20); assert.equal(p.seed.intermediate, false);
  assert.equal(p.seed.converged, false); assert.equal(p.diagnostics.targetConverged, false);
  physicalEqual(replay(cp).value, replay(p.checkpoint).value);
  assert.throws(() => prepareCoupledNcritRestart(21, p.checkpoint), /source must already satisfy/);
});

test('omitted exact and explicit native floor policies preserve their original checkpoint semantics', () => {
  const native = sourceRoot(), exact = structuredClone(native); delete exact.restart.options.hkFloorLinearization;
  const a = prepareCoupledNcritRestart(5, exact), b = prepareCoupledNcritRestart(5, native);
  assert.equal(Object.hasOwn(a.checkpoint.restart.options, 'hkFloorLinearization'), false);
  assert.equal(b.checkpoint.restart.options.hkFloorLinearization, 'native');
  assert.deepEqual(a.checkpoint.families, b.checkpoint.families);
  assert.deepEqual(a.checkpoint.restart.initialEuler, b.checkpoint.restart.initialEuler);
  assert.deepEqual(a.checkpoint.restart.initialBL, b.checkpoint.restart.initialBL);
});

test('an explicit Armijo checkpoint keeps its update policy through criterion preparation and exact resume', () => {
  const cp = sourceRoot(); cp.continuation.stepAcceptance = 'armijo';
  const before = structuredClone(cp), p = prepareCoupledNcritRestart(5, cp, { requestedNcrit: 9 });
  assert.deepEqual(cp, before); assert.deepEqual(p.checkpoint.continuation, before.continuation);
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = p.checkpoint.continuation;
  const resumed = solveCoupledStreamtubeIses(undefined, { resume: copy(p.checkpoint), maxIterations: 0,
    tolerance, iterationGeometry, stepAcceptance, stagnationLimiter });
  assert.equal(resumed.linearDiagnostics.solves, 0); assert.equal(resumed.stepAcceptance, 'armijo');
  assert.deepEqual(copy(resumed.checkpoint), copy(p.checkpoint));
  assert.deepEqual(resumed.families, p.diagnostics.targetFamilies);
  const unrecognized = structuredClone(cp); unrecognized.continuation.stepAcceptance = 'unqualified-policy';
  assert.throws(() => prepareCoupledNcritRestart(5, unrecognized), /complete ISES maintenance history/);
});

test('invalid direction, incomplete source and forged convergence fail without mutating the checkpoint', () => {
  const cp = sourceRoot();
  for (const [target, options] of [[4, {}], [3, {}], [5, { requestedNcrit: 4.5 }], [Infinity, {}], [5, { tolerance: 0 }], [5, { blPredictor: 'unknown' }]]) {
    const before = structuredClone(cp);
    assert.throws(() => prepareCoupledNcritRestart(target, cp, options), CoupledNcritRestartError);
    assert.deepEqual(cp, before);
  }
  for (const mutate of [
    x => { delete x.continuation; },
    x => { delete x.restart.initialBL; },
    x => { x.restart.options.tripFractions[0][0] = .8; },
    x => { x.restart.options.transitionMode = 'fixed-trip'; },
    x => { x.families.boundaryLayer = .1; },
    x => { x.families.euler = 0; },
    x => { x.restart.initialEuler.x[0] = NaN; },
    x => { x.restart.initialEuler.nodes[0][2][1].x = x.restart.initialEuler.nodes[0][4][1].x; },
    x => { x.restart.initialBL[3] = 100; },
    x => { x.restart.options.unsupportedPhysics = true; },
  ]) {
    const broken = structuredClone(cp); mutate(broken); const before = structuredClone(broken);
    assert.throws(() => prepareCoupledNcritRestart(5, broken), error => error instanceof CoupledNcritRestartError
      && error.diagnostics.targetConverged === false && error.diagnostics.stage === 'source validation');
    assert.deepEqual(broken, before);
  }
});

test('a real browser prepares the complete saved state through HTTP modules and WASM', {
  skip: process.env.MSES_BROWSER_TESTS !== '1',
}, async () => {
  // The existing KLU module selects HTTP fetch in browsers and a dynamic
  // node:fs import only for file URLs. Exercise its real HTTP path, not a
  // bundler that incorrectly resolves the unexecuted Node branch.
  const root = path.resolve(fileURLToPath(new URL('..', import.meta.url)));
  const server = createServer(async (request, response) => {
    try {
      const pathname = new URL(request.url, 'http://localhost').pathname;
      if (pathname === '/') { response.setHeader('Content-Type', 'text/html'); response.end('<!doctype html><title>Ncrit draft test</title>'); return; }
      const file = path.resolve(root, '.' + pathname);
      if (!file.startsWith(root + path.sep)) { response.writeHead(403).end(); return; }
      response.setHeader('Content-Type', file.endsWith('.wasm') ? 'application/wasm' : 'text/javascript');
      response.end(await readFile(file));
    } catch { response.writeHead(404).end(); }
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  let browser;
  try {
    browser = await chromium.launch({ headless: true });
    const page = await browser.newPage(), errors = []; page.on('pageerror', error => errors.push(error.message));
    await page.goto(`http://127.0.0.1:${server.address().port}/`);
    const cp = sourceRoot();
    const result = await page.evaluate(async checkpoint => {
      const { prepareCoupledNcritRestart } = await import('/src/euler/streamtube-coupled-ncrit-restart.js');
      const p = prepareCoupledNcritRestart(5, checkpoint, { requestedNcrit: 9 });
      const native = prepareCoupledNcritRestart(5, checkpoint, { requestedNcrit: 9, blPredictor: 'xfoil-mrchdu' });
      return { checkpoint: p.checkpoint, diagnostics: p.diagnostics, native };
    }, cp);
    assert.deepEqual(errors, []); assert.equal(result.diagnostics.prepared, true);
    assert.equal(result.diagnostics.operations.globalLinearSolves, 0);
    assert.deepEqual(copy(result.checkpoint), copy(prepareCoupledNcritRestart(5, cp, { requestedNcrit: 9 }).checkpoint));
    assert.deepEqual(copy(result.native), copy(prepareCoupledNcritRestart(5, cp, { requestedNcrit: 9, blPredictor: 'xfoil-mrchdu' })));
  } finally {
    await browser?.close(); await new Promise(resolve => server.close(resolve));
  }
});

test('Ncrit preparation accepts verified change convergence and invalidates that evidence for the new target', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2, mach: .2 }),
    streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const source = solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity',
    blThermodynamics: 'historical-common-isentrope', transitionMode: 'automatic', tripFractions: [[1, 1], [1, 1]],
    hkFloorLinearization: 'native', ncrit: 4, maxIterations: 12, tolerance, stepAcceptance: 'admissible', convergence: 'mses' });
  assert.equal(source.converged, true, source.reason);
  assert.equal(source.reason, 'solution changes');
  assert.equal(source.residualConverged, false);
  const before = structuredClone(source.checkpoint);
  const next = prepareCoupledNcritRestart(5, source.checkpoint);
  assert.equal(next.checkpoint.continuation.convergence, 'mses');
  assert.equal(next.checkpoint.convergence, undefined);
  assert.equal(next.seed.converged, false);
  assert.deepEqual(source.checkpoint, before);
  const unproven = structuredClone(source.checkpoint); delete unproven.convergence;
  assert.throws(() => prepareCoupledNcritRestart(5, unproven), /source must already satisfy/);
});
