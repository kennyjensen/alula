// SPDX-License-Identifier: GPL-2.0-or-later
// Real retained BL physics with a controlled target-construction seam.
// This proves the phase contract, not a Mach-only occurrence of stale N.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { pathToFileURL } from 'node:url';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';

const sourcePath = 'docs/rae2822/mses063875-recovered-source-tail/accepted-parent.json';
const checkpoint = JSON.parse(fs.readFileSync(sourcePath)).checkpoint;
const runtime = 'src/euler/tests/streamtube-coupled-flow-restart.js';
const archive = 'docs/coupled-preserve-amplification-reconciliation/before/streamtube-coupled-flow-restart.js.txt';
const encode = text => 'data:text/javascript;base64,' + Buffer.from(text).toString('base64');
const serial = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
let sequence = 0;
async function adapter(file, perturbTarget) {
  const calls = { constructors: 0, evaluations: 0, admissibilityCalls: 0, nativeMarches: 0 }, targets = [];
  const key = `__preserveAmplification${sequence++}`;
  globalThis[key] = {
    create(input, options) {
      calls.constructors++;
      const system = createCoupledStreamtubeBody(input, options), evaluate = system.evaluate;
      system.evaluate = (...args) => { calls.evaluations++; return evaluate(...args); };
      const admissible = system.admissible;
      system.admissible = (...args) => { calls.admissibilityCalls++; return admissible(...args); };
      if (input.mach !== checkpoint.restart.input.mach && perturbTarget) {
        const surface = system.bl.surfaces[0], upstream = surface.ids[surface.transition - 1];
        // Controlled stale target auxiliary only. The accepted source is
        // reconstructed and fully checked before this target is constructed.
        system.initial[system.ne + 4 * upstream] = 0;
        const raw = system.initial.slice(), phase = system.bl.snapshotActive();
        const target = system.bl.activeTargets(raw.subarray(0, system.ne), raw.subarray(system.ne));
        assert(target.every(t => t.from === t.to));
        assert.equal(target[0].amplificationReconciliation, true);
        assert.equal(target[1].amplificationReconciliation, undefined);
        assert.throws(() => system.evaluate(raw), /Transition is outside this active interval/);
        targets.push({ system, raw, phase, target, upstream });
      }
      return system;
    },
    native() { calls.nativeMarches++; throw Error('No native profile march is permitted in this test'); },
  };
  const replacements = {
    './streamtube-coupled.js': encode(`export const createCoupledStreamtubeBody=globalThis['${key}'].create;`),
    './streamtube-coupled-mrchdu-predictor.js': encode(`export const prepareCoupledMrchduProfiles=globalThis['${key}'].native;`),
  };
  const code = fs.readFileSync(file, 'utf8').replace(/from '([^']+)'/g, (_, specifier) => {
    const original = path.resolve(file === runtime ? 'src/euler/tests' : 'src/euler', specifier);
    const resolved = fs.existsSync(original) ? original : path.resolve('src/euler/tests', specifier);
    return `from '${replacements[specifier] ?? pathToFileURL(resolved).href}'`;
  });
  try { return { module: await import(encode(code + `\n//# sourceURL=preserve-reconciliation-${key}.js`)), calls, targets }; }
  finally { delete globalThis[key]; }
}

test('preserved target prepares a flagged same-index N prefix before evaluating its real mixed interval', async t => {
  const cp = structuredClone(checkpoint), before = structuredClone(cp), targetMach = .65;
  const old = await adapter(archive, true);
  assert.throws(() => old.module.initializeCoupledStreamtubeFromFlow(targetMach, cp, { blPredictor: 'preserve' }),
    /Transition is outside this active interval/);
  assert.deepEqual(cp, before);
  const current = await adapter(runtime, true);
  const prepared = current.module.initializeCoupledStreamtubeFromFlow(targetMach, cp, { blPredictor: 'preserve' });
  const { system, raw, phase, target } = current.targets[0], change = prepared.diagnostics.targetPhaseInitialization;
  assert.equal(change.changes.length, 1);
  assert.deepEqual([change.changes[0].body, change.changes[0].side, change.changes[0].kind, change.changes[0].auxiliaryOnly],
    [0, 'upper', 'amplification-reconciliation', true]);
  assert.deepEqual(change.before, phase); assert.deepEqual(change.after, phase);
  assert.deepEqual(prepared.checkpoint.restart.options.transitionState, cp.restart.options.transitionState);
  const allowed = new Set(system.bl.surfaces[0].ids.slice(0, phase[0]).map(id => system.ne + 4 * id));
  let changed = 0;
  prepared.initial.forEach((value, i) => {
    if (!allowed.has(i)) assert.equal(value, raw[i], `Non-target auxiliary or primary field ${i} changed`);
    else if (value !== raw[i]) changed++;
  });
  assert(changed > 0); assert.equal(prepared.diagnostics.physicalGeometryPreserved, true);
  assert.equal(prepared.diagnostics.physicalDensityPreserved, true); assert.equal(prepared.diagnostics.physicalMassPreserved, true);
  assert.equal(prepared.diagnostics.transitionMapPreserved, true); assert.equal(prepared.diagnostics.packedBLPreserved, false);
  assert.deepEqual(prepared.checkpoint.restart.initialEuler, cp.restart.initialEuler);
  assert.deepEqual(prepared.checkpoint.continuation, cp.continuation);
  assert.deepEqual(prepared.system.initial, prepared.initial);
  assert(prepared.value.residual.every(Number.isFinite));
  const next = system.bl.activeTargets(prepared.initial.subarray(0, system.ne), prepared.initial.subarray(system.ne));
  assert(next.every(t => t.from === t.to && !t.amplificationReconciliation));
  const stable = prepared.initial.slice();
  assert.deepEqual(system.bl.updateActive(prepared.initial.subarray(system.ne), prepared.initial.subarray(0, system.ne)), { changed: false, changes: [] });
  assert.deepEqual(prepared.initial, stable);
  const f = prepared.checkpoint.restart, replay = createCoupledStreamtubeBody(f.input, {
    ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  assert.deepEqual(replay.initial, prepared.initial);
  assert.deepEqual(replay.evaluate(replay.initial).residual, prepared.value.residual);
  assert.deepEqual(cp, before); assert.equal(old.calls.nativeMarches + current.calls.nativeMarches, 0);
  t.diagnostic(JSON.stringify({ source: sourcePath, sourceMach: cp.restart.input.mach, targetMach,
    scope: 'Target auxiliary perturbation is controlled; no Mach-only reproduction is claimed.',
    phase, flaggedSurfaces: target.flatMap((v, k) => v.amplificationReconciliation ? [k] : []), changedAuxiliaries: changed,
    constructors: old.calls.constructors + current.calls.constructors + 1,
    explicitFullEvaluations: old.calls.evaluations + current.calls.evaluations + 1,
    admissibilityEvaluations: old.calls.admissibilityCalls + current.calls.admissibilityCalls,
    globalJacobians: 0, linearSolves: 0, newtonUpdates: 0, nativeProfileMarches: 0,
    targetFamilies: prepared.value.families, targetCheckpointReplayExact: true }));
});

test('valid same-Mach and unflagged changed-Mach preserve transfers remain exactly archived', async () => {
  const current = await adapter(runtime, false), old = await adapter(archive, false);
  const extract = p => serial({ state: p.initial, residual: p.value.residual, checkpoint: p.checkpoint, diagnostics: p.diagnostics });
  const cp = structuredClone(checkpoint), before = structuredClone(cp);
  for (const mach of [cp.restart.input.mach, .65]) {
    const expected = old.module.initializeCoupledStreamtubeFromFlow(mach, cp, { blPredictor: 'preserve' });
    const actual = current.module.initializeCoupledStreamtubeFromFlow(mach, cp, { blPredictor: 'preserve' });
    assert.deepEqual(extract(actual), extract(expected));
    assert.equal(actual.diagnostics.targetPhaseInitialization, undefined);
  }
  assert.deepEqual(cp, before); assert.equal(current.calls.nativeMarches + old.calls.nativeMarches, 0);
});
