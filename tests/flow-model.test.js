import test from 'node:test';
import assert from 'node:assert/strict';
import { flowModelRoute } from '../src/ui/flow-model.js';

test('each explicit selector has one canonical Worker solver and panels remain incompressible', () => {
  assert.deepEqual(flowModelRoute('inviscid', .74), { flowModel: 'inviscid', mach: 0 });
  assert.deepEqual(flowModelRoute('coupled', NaN), { flowModel: 'coupled', mach: 0 });
  assert.deepEqual(flowModelRoute('streamtube-grid', .28), { flowModel: 'streamtube-grid', mach: .28, quadBoundaryLayers: false });
  assert.deepEqual(flowModelRoute('streamtube-bl', .35), { flowModel: 'streamtube-grid', mach: .35, quadBoundaryLayers: true });
});

test('removed or invalid UI choices cannot fall through to another solver', () => {
  for (const mode of ['subcritical', '', undefined, 'euler']) assert.throws(() => flowModelRoute(mode), /four supported/);
  for (const mach of [NaN, Infinity, 0, -.1]) assert.throws(() => flowModelRoute('streamtube-bl', mach), /positive/);
});

