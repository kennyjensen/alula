// SPDX-License-Identifier: GPL-2.0-or-later
// Controller-only fixtures: exact geometry and BL are checked separately on
// the real saved RAE/NLR states, not by these injected numerical outcomes.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';

let serial = 0;
async function fixture(outcomes) {
  const key = `__refinedProfile${serial++}`, calls = [];
  globalThis[key] = (input, options) => { calls.push(structuredClone({ input, options }));
    const result = outcomes[calls.length - 1]; if (result instanceof Error) throw result; return result; };
  const source = fs.readFileSync(process.env.MSES_GEOMETRY_REPLAY_RUNTIME
    ? path.join(process.env.MSES_GEOMETRY_REPLAY_RUNTIME, 'src/euler/streamtube-coupled-grid-profile.js')
    : new URL('../src/euler/streamtube-coupled-grid-profile.js', import.meta.url), 'utf8');
  const isolated = `const prepareCoupledGridProfile = globalThis[${JSON.stringify(key)}];\n`
    + source.slice(source.indexOf('export function prepareRefinedCoupledGridProfile('));
  const { prepareRefinedCoupledGridProfile } = await import(`data:text/javascript;base64,${Buffer.from(isolated).toString('base64')}`);
  delete globalThis[key]; return { calls, prepareRefinedCoupledGridProfile };
}
const initial = () => ({ input: { mach: .2 }, options: { ncrit: 5.5 },
  initialEuler: { x: [0, 1], nodes: [[{ x: 1, y: .3 }]], undisplacedNodes: [[{ x: 1, y: .2 }]] },
  initialBL: [.02, .001, .002, 1] });
const error = (failedInvariant = 'source physical nodes', stage = 'source replay') =>
  Object.assign(new Error('Declared numerical replay rejection'), { code: 'COUPLED_GRID_PROFILE_PREPARATION',
    diagnostics: { stage, failedInvariant } });

test('a successful existing refinement preparation is returned unchanged without another attempt', async () => {
  const expected = { sentinel: true }, h = await fixture([expected]);
  const input = initial(), before = structuredClone(input), controls = { wakeWidthIncrement: true };
  assert.equal(h.prepareRefinedCoupledGridProfile(input, controls), expected);
  assert.deepEqual(input, before); assert.deepEqual(h.calls, [{ input, options: controls }]);
});

test('only the typed new-source physical-node invariant enables a chart-preserving retry', async () => {
  const failure = error(), result = { transfer: { profilePreparation: { canonicalReplayExact: true } } };
  const h = await fixture([failure, result]), input = initial(), before = structuredClone(input);
  assert.equal(h.prepareRefinedCoupledGridProfile(input, { reinitializeAmplification: true }), result);
  assert.deepEqual(input, before); assert.equal(h.calls.length, 2);
  assert.deepEqual(h.calls[1].input, { ...input, options: { ...input.options, geometryReplay: 'preserve-undisplaced' } });
  const receipt = result.transfer.profilePreparation.geometryReplayRecovery;
  assert.equal(receipt.accepted, true); assert.equal(receipt.exactSourceGeometryRequired, true);
  assert.equal(receipt.originalGatesChanged, false); assert.deepEqual(receipt.sourceFailure.diagnostics, failure.diagnostics);
  assert.deepEqual(h.calls[1].options, h.calls[0].options);
});

test('packed-state, invalid geometry, native/target replay and unknown exceptions are never retried', async () => {
  for (const failure of [error('source packed state'), error('canonical full residual', 'canonical replay'),
    error(undefined, 'complete target domain'), new Error('source physical nodes'),
    Object.assign(error(), { code: 'observer-cancellation' })]) {
    const h = await fixture([failure]);
    assert.throws(() => h.prepareRefinedCoupledGridProfile(initial()), e => e === failure);
    assert.equal(h.calls.length, 1);
  }
  const failure = error(), h = await fixture([failure]), input = initial();
  input.options.geometryReplay = 'preserve-undisplaced';
  assert.throws(() => h.prepareRefinedCoupledGridProfile(input), e => e === failure);
  assert.equal(h.calls.length, 1);
});

test('an unsuccessful exact retry stays rejected and retains both failure diagnostics', async () => {
  const first = error(), second = error('canonical full residual', 'canonical replay');
  const h = await fixture([first, second]), input = initial(), before = structuredClone(input);
  assert.throws(() => h.prepareRefinedCoupledGridProfile(input), e => {
    assert.equal(e, second); const receipt = e.diagnostics.geometryReplayRecovery;
    assert.equal(receipt.accepted, false); assert.deepEqual(receipt.sourceFailure.diagnostics, first.diagnostics);
    assert.equal(e.diagnostics.failedInvariant, 'canonical full residual'); return true;
  });
  assert.deepEqual(input, before); assert.equal(h.calls.length, 2);
});
