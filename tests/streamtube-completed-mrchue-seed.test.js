import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { initializeCoupledStreamtubeBody } from '../src/euler/streamtube-coupled-initializer.js';

const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/naca-completed-mrchue-seed.json', import.meta.url)));

test('a completed original NACA MRCHUE march retains the exact cold state and full initialization diagnostics', async () => {
  const sourceURL = new URL('../src/euler/streamtube-coupled-initializer.js', import.meta.url);
  const archived = fs.readFileSync(new URL('../docs/solver-reliability/te-wake-seed-improvement/initializer.before.js.txt', import.meta.url), 'utf8')
    .replace(/from '(\.\/[^']+)'/g, (_, path) => `from '${new URL(path.replace('./streamtube-displacement-initializer.js', './streamtube-displacement.js'), sourceURL).href}'`);
  const old = await import('data:text/javascript;base64,' + Buffer.from(archived).toString('base64'));
  const { input, options } = structuredClone(fixture);
  options.initialEuler.x = Float64Array.from(options.initialEuler.x);
  const before = structuredClone({ input, options });
  const expected = old.initializeCoupledStreamtubeBody(input, options);
  const actual = initializeCoupledStreamtubeBody(input, options);
  const rows = actual.initialization.originalBLInitialization.boundaryLayer;
  const surfaceReceipts = rows.filter(row => row.method === 'mrchue' && Array.isArray(row.localConvergenceWarnings));
  assert.equal(surfaceReceipts.length, 2);
  assert.ok(surfaceReceipts.every(row => row.localConvergenceWarnings.length === 0));
  assert.equal(actual.initialization.wakeInitialization, undefined);
  assert.equal(actual.initialization.thicknessFactor, 1);
  assert.equal(actual.initialization.history.length, 1);
  assert.equal(actual.mesh.quality.valid, true);
  assert.deepEqual(actual.initialization, expected.initialization);
  assert.deepEqual(actual.system.initial, expected.system.initial);
  assert.deepEqual(actual.system.conditions, expected.system.conditions);
  assert.deepEqual(actual.system.bl.snapshotActive(), expected.system.bl.snapshotActive());
  assert.deepEqual(actual.mesh, expected.mesh);
  const actualFlow = actual.system.evaluate(actual.system.initial), expectedFlow = expected.system.evaluate(expected.system.initial);
  assert.deepEqual(actualFlow.outer.nodes, expectedFlow.outer.nodes);
  assert.deepEqual(actualFlow.outer.undisplacedNodes, expectedFlow.outer.undisplacedNodes);
  assert.deepEqual(actualFlow.residual, expectedFlow.residual);
  assert.deepEqual(actualFlow.families, expectedFlow.families);
  assert.deepEqual({ input, options }, before, 'Neither initializer may mutate the submitted cold state.');
});
