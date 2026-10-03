// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { createBrowserCheckpointJournal } from '../scripts/validation/browser-checkpoint-journal.js';
import { coldBrowserCheckpointObservers } from '../scripts/validation/cold-browser-checkpoint-worker.js';

const fixture = (sequence, iteration) => ({ sequence, type: 'coupled-checkpoint', workerSeconds: sequence / 10,
  checkpoint: { restart: { options: { ncrit: 4 }, initialBL: new Float64Array([.002, .0005]) },
    continuation: { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', history: [{ iteration }] },
    families: { euler: .01, boundaryLayer: .02, edgeMatching: .03 } },
  details: { history: [{ iteration }], startupAttempt: 2, stage: 'coupled' } });

test('an interrupted observation retains a complete last checkpoint, continuation, and periodic archives', t => {
  const output = fs.mkdtempSync(path.join(os.tmpdir(), 'mses-browser-journal-'));
  t.after(() => fs.rmSync(output, { recursive: true, force: true }));
  const journal = createBrowserCheckpointJournal(output), first = fixture(1, 0), before = structuredClone(first);
  journal.record(first); journal.record(fixture(2, 1)); journal.record(fixture(3, 10));
  assert.deepEqual(first, before, 'observing must not mutate the supplied checkpoint');
  const load = name => JSON.parse(fs.readFileSync(path.join(output, name + '.json')));
  assert.deepEqual(load('latest-checkpoint').checkpoint.restart.initialBL, [.002, .0005]);
  assert.equal(load('latest-checkpoint').checkpoint.continuation.history.at(-1).iteration, 10);
  assert.equal(load('initial-checkpoint').details.history.at(-1).iteration, 0);
  assert.deepEqual(fs.readdirSync(path.join(output, 'checkpoints')), ['000001.json', '000003.json']);
  assert.equal(load('journal').lastSequence, 3);
  assert.equal(load('journal').checkpointCount, 3);
  // A malformed/reordered subsequent message cannot replace the durable state
  // or claim an iteration which was not completely written.
  assert.throws(() => journal.record(fixture(5, 6)), /Missing or reordered/);
  assert.equal(load('journal').lastSequence, 3);
  const incomplete = fixture(4, 6); delete incomplete.checkpoint.continuation;
  assert.throws(() => journal.record(incomplete), /incomplete/);
  assert.equal(load('latest-checkpoint').details.history.at(-1).iteration, 10);
});

test('Euler handoff and terminal failure survive without any coupled initialization', t => {
  const output = fs.mkdtempSync(path.join(os.tmpdir(), 'mses-browser-journal-'));
  t.after(() => fs.rmSync(output, { recursive: true, force: true }));
  const journal = createBrowserCheckpointJournal(output);
  journal.record({ sequence: 1, type: 'prepared-euler', packet: { mesh: { nodes: [1, 2] } } });
  journal.record({ sequence: 2, type: 'assembly-checkpoint', observation: { stage: 'euler',
    result: { flow: { checkpoint: { initial: [2, 3] }, bestCheckpoint: { initial: [4, 5] } } } } });
  journal.record({ sequence: 3, type: 'error', message: 'typed precursor failure', inputUnchanged: true });
  const read = name => JSON.parse(fs.readFileSync(path.join(output, name + '.json')));
  assert.deepEqual(read('assembly-000002').observation.result.flow.bestCheckpoint.initial, [4, 5]);
  assert.equal(read('terminal').message, 'typed precursor failure');
  assert.equal(read('journal').checkpointCount, 0);
  assert.equal(read('journal').terminal, 'error');
});

test('prepared Euler observation strips only the process-local system and all hook messages are cloneable', () => {
  const messages = [], observe = coldBrowserCheckpointObservers(message => messages.push(structuredClone(message)));
  const source = { version: 1, caseData: { mach: .2 }, prepared: {
    system: { evaluate: () => {}, geometryChart: () => {}, conditions: { mach: .2 } },
    input: { mach: .2 }, initial: new Float64Array([1, 2]), initialEuler: { x: new Float64Array([1, 2]) },
    nodes: [[[{ x: 0, y: 0 }]]], diagnostics: { gasInitialization: undefined }, mesh: { quality: { valid: true } },
  } };
  assert.throws(() => structuredClone(source), { name: 'DataCloneError' });
  observe.onEulerPrepared(source);
  assert.equal(source.prepared.system.conditions.mach, .2);
  assert.equal(typeof source.prepared.system.evaluate, 'function');
  assert.equal(messages[0].packet.prepared.system, undefined);
  assert.deepEqual(messages[0].packet.prepared.initial, source.prepared.initial);
  assert.deepEqual(messages[0].packet.prepared.nodes, source.prepared.nodes);
  assert.equal(messages[0].processLocalSystemOmitted, true);
  // Public assembly checkpoints contain restart/initialization data, not the
  // initializer's live system; complete iteration checkpoints add continuation.
  const complete = fixture(1, 0);
  observe.onCheckpoint({ stage: 'boundary-layer-initialization', restart: complete.checkpoint.restart,
    initialization: { method: 'MRCHUE surfaces and ISET-style wake guess', history: [{ accepted: true }] } });
  observe.onIterationCheckpoint(complete.checkpoint, complete.details);
  assert.deepEqual(messages[2].checkpoint.continuation, complete.checkpoint.continuation);
  assert.deepEqual(messages[1].observation.restart, complete.checkpoint.restart);
});
