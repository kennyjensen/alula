// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';

// Compile only the producer's marked pure loader block. Its CLI, Worker,
// numerical imports and flow driver are never executed by these tests.
const source = fs.readFileSync('scripts/check-coupled-shock-broadening.js', 'utf8');
const start = source.indexOf('// BEGIN PURE SAVED-INPUT HELPERS');
const end = source.indexOf('// END PURE SAVED-INPUT HELPERS');
assert.ok(start >= 0 && end > start);
const { loadShockBroadeningInput: load, requireAdjacentShockCheckpoints: adjacent } = new Function('assert',
  `${source.slice(start, end)}\nreturn { loadShockBroadeningInput, requireAdjacentShockCheckpoints };`)(assert);
const copy = value => structuredClone(value);
const checkpoint = residual => ({ version: 1, restart: { input: {}, options: {}, initialEuler: { x: [0] }, initialBL: [1] },
  continuation: {}, families: { euler: residual, boundaryLayer: .1, edgeMatching: .01 } });
const entry = (iteration, residual) => ({ iteration, step: iteration ? 1 : 0, residual,
  euler: residual, boundaryLayer: .1, edgeMatching: .01 });
const pair = () => {
  const history = [entry(0, 3), entry(1, 2)];
  return [{ checkpoint: checkpoint(2), details: { history } },
    { checkpoint: checkpoint(1), details: { history: [...copy(history), entry(2, 1)] } }];
};
const controller = () => ({ version: 1, kind: 'mses-shock-broadening', checkpoint: checkpoint(2),
  targetMcrit: .99, effectiveMcrit: .8, previousDensityChange: .161, acceptedUpdates: 4,
  tolerance: 1e-8, nextAction: 'target-check', policy: { formula: 'preserved' } });

test('plain checkpoints load, while own result/state/event wrappers retain complete controller history', () => {
  const cp = checkpoint(2), saved = controller(), before = copy(saved);
  assert.equal(load(cp), cp);
  assert.equal(load({ checkpoint: cp }), cp);
  assert.equal(load({ result: { checkpoint: cp } }), cp);
  for (const input of [saved, { controller: saved, checkpoint: saved.checkpoint },
    { type: 'result', result: { controller: saved, checkpoint: saved.checkpoint } }]) {
    assert.equal(load(input), saved);
    assert.deepEqual(load(input), before);
  }
  assert.deepEqual(saved, before);
});

test('malformed or conflicting saved controllers cannot fall back to a plain checkpoint and reset d', () => {
  const cp = checkpoint(2), saved = controller();
  for (const input of [{ controller: null, checkpoint: cp }, { controller: { checkpoint: cp }, checkpoint: cp },
    { type: 'result', result: { controller: null, checkpoint: cp } }])
    assert.throws(() => load(input), /controller.*complete/);
  assert.throws(() => load({ controller: saved, result: { controller: { ...saved, previousDensityChange: .2 } } }), /Conflicting/);
  assert.throws(() => load({ type: 'timeout' }), /complete coupled checkpoint/);
});

test('adjacent accepted histories require exact prefix, contiguous iterations and checkpoint-family binding', () => {
  const [previous, current] = pair(), before = copy([previous, current]);
  const proof = adjacent(previous, current);
  assert.deepEqual(proof, { previousIteration: 1, currentIteration: 2,
    previousHistoryPath: 'details.history', currentHistoryPath: 'details.history',
    previousHistoryEntries: 2, currentHistoryEntries: 3, historyPrefixExact: true, checkpointFamiliesExact: true });
  const previousH = { checkpoint: previous.checkpoint, h: previous.details.history };
  const currentH = { checkpoint: current.checkpoint, history: current.details.history, h: { iteration: 2 } };
  assert.equal(adjacent(previousH, currentH).historyPrefixExact, true);
  assert.deepEqual([previous, current], before);
});

test('missing, unrelated, skipped or contradictory density histories are rejected before a flow can start', () => {
  const [previous, current] = pair();
  assert.throws(() => adjacent({ checkpoint: previous.checkpoint }, current), /complete iteration histories/);
  assert.throws(() => adjacent({ checkpoint: previous.checkpoint, h: { iteration: 1 } }, current), /complete iteration histories/);
  const unrelated = copy(current); unrelated.details.history[0].step = .5;
  assert.throws(() => adjacent(previous, unrelated), /exact prefix/);
  const skipped = copy(current); skipped.details.history.at(-1).iteration = 3;
  assert.throws(() => adjacent(previous, skipped), /contiguous/);
  const stale = copy(current); stale.checkpoint.families.euler = 1.5;
  assert.throws(() => adjacent(previous, stale), /does not match its checkpoint/);
  assert.throws(() => adjacent(previous, { ...current, h: { iteration: 99 } }), /h iteration disagrees/);
  assert.throws(() => adjacent(previous, { ...current, history: [entry(0, 1)] }), /Conflicting retained/);
  assert.throws(() => adjacent(controller(), current), /already owns its density history/);
});

test('completed four-update control source pair has a verifiable adjacent accepted-history prefix', () => {
  const directory = 'docs/temporary-shock-broadening/fine07146875-step7';
  const report = JSON.parse(fs.readFileSync(`${directory}/report.json`));
  const read = path => JSON.parse(fs.readFileSync(path));
  const proof = adjacent(read(report.densityHistory.previous), read(report.densityHistory.current));
  assert.equal(proof.previousIteration, 5);
  assert.equal(proof.currentIteration, 6);
  assert.equal(proof.previousHistoryEntries, 6);
  assert.equal(proof.currentHistoryEntries, 7);
  assert.equal(proof.historyPrefixExact, true);
  const result = read(`${directory}/result.json`), state = read(`${directory}/state.json`);
  assert.deepEqual(load(result), result.controller);
  assert.deepEqual(load(state), result.controller);
  assert.equal(load(result).previousDensityChange, result.controller.previousDensityChange);
});
