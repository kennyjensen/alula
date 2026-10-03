// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { selectSurfaceTransition } from '../src/viscous/transition-selection.js';
import { checkAutomaticTransition } from '../src/viscous/transition-interval.js';
import { createCoupledStreamtubeBody, coupledStreamtubeTripEvents } from '../src/euler/streamtube-coupled.js';

test('an unchanged marched transition index does not imply a valid packed-N mixed interval', () => {
  const reference = JSON.parse(fs.readFileSync(new URL('fixtures/fortran/automatic-transition.json', import.meta.url)));
  const profile = reference.profiles.find(p => p.parameters.ncrit === 9 && !p.tripS);
  const kernel = createIntegralKernel({ ...profile.parameters, exactJacobian: true, transitionTolerance: 1e-12 });
  const selected = selectSurfaceTransition(kernel, profile.states), j = selected.index;
  const states = profile.states.map((s, i) => ({ ...s, aux: i < j ? selected.amplification[i] : s.aux }));
  assert.equal(checkAutomaticTransition(kernel, { upstream: states[j - 1], downstream: states[j] }).transition, true);
  states[j - 1].aux = 0;
  assert.ok(states.slice(0, j).every(s => s.aux < kernel.parameters.ncrit));
  assert.equal(selectSurfaceTransition(kernel, states).index, j);
  assert.equal(checkAutomaticTransition(kernel, { upstream: states[j - 1], downstream: states[j] }).transition, false);
});

test('an incompatible undercritical N prefix is prepared atomically without changing the coupled equations', () => {
  const cp = JSON.parse(fs.readFileSync(new URL('../docs/rae2822/mses063875-recovered-source-tail/accepted-parent.json', import.meta.url))).checkpoint;
  const f = cp.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const x = system.initial.slice(), initial = x.slice(), phase = system.bl.snapshotActive();
  const value = system.evaluate(x), events = coupledStreamtubeTripEvents(system, { blUpdate: 'xfoil' });
  assert.deepEqual(value.families, cp.families);
  assert.deepEqual(events.prepare(x, initial), { changed: false });
  assert.deepEqual(x, initial, 'Valid solved amplification remains exactly unchanged.');
  const surface = system.bl.surfaces[0], upstream = surface.ids[surface.transition - 1];
  x[system.ne + 4 * upstream] = 0;
  const before = x.slice();
  const targets = system.bl.activeTargets(x.subarray(0, system.ne), x.subarray(system.ne));
  assert.ok(targets.every(t => t.from === t.to));
  assert.ok(surface.ids.slice(0, surface.transition).every(id => x[system.ne + 4 * id] < f.options.ncrit));
  assert.throws(() => system.evaluate(x), /Transition is outside this active interval/);
  const event = events.prepare(x, before);
  assert.equal(event.changed, true);
  assert.equal(event.changes.length, 1);
  const change = event.changes[0];
  assert.deepEqual([change.body, change.side, change.kind, change.auxiliaryOnly], [0, 'upper', 'amplification-reconciliation', true]);
  assert.deepEqual(system.bl.snapshotActive(), phase);
  const allowed = new Set(surface.ids.slice(0, surface.transition).map(id => system.ne + 4 * id));
  x.forEach((v, i) => { if (!allowed.has(i)) assert.equal(v, before[i], `Unexpected change to unknown ${i}.`); });
  const repaired = system.evaluate(x);
  assert.deepEqual(repaired.outer.nodes, value.outer.nodes);
  assert.ok(repaired.residual.every(Number.isFinite));
  assert.ok(Math.max(...surface.ids.slice(0, surface.transition).map(id => Math.abs(repaired.residual[system.ne + 4 * id]))) < 1e-10);
  assert.equal(system.admissible(x), true);
  const stable = x.slice();
  assert.deepEqual(events.prepare(x, stable), { changed: false });
  assert.deepEqual(x, stable);
});
