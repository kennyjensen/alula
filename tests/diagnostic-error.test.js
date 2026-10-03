// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { serializeDiagnosticError } from '../scripts/validation/diagnostic-error.js';

test('a nested wake initialization failure retains the typed local replay data', () => {
  const upstream = { s: .31, ue: .95, aux: .04, theta: .001, deltaStar: .002 };
  const inner = Object.assign(new Error('BL station inverse initialization failed (wake, s=0.3317364308903284): iteration limit'), {
    code: 'BL_STATION_INVERSE_INITIALIZATION', diagnostics: { regime: 'wake', reynolds: 1e6, mach: .2,
      gamma: 1.4, targetHK: 1.3, upstream, requestedUe: .94, originalSeed: { ...upstream, s: .3317364308903284 },
      directState: { ...upstream }, retainedState: { ...upstream, aux: .05 }, directReason: 'line search failed',
      history: [{ iteration: 30, residual: .1, step: .3 }],
      inverseRecovery: { method: 'xfoil-mrchue-bounded-inverse', history: [{ iteration: 30, residual: .2 }] } },
  });
  const outer = new Error('BL initialization failed at body 2, wake station 492', { cause: inner });
  const result = serializeDiagnosticError(outer);
  assert.equal(result.message, outer.message);
  assert.equal(result.cause.code, inner.code);
  assert.deepEqual(result.cause.diagnostics, inner.diagnostics);
  assert.deepEqual(JSON.parse(JSON.stringify(result)), result);
  upstream.theta = 99;
  assert.equal(result.cause.diagnostics.upstream.theta, .001);
});

test('cyclic and non-JSON diagnostics cannot break error reporting', () => {
  const error = new Error('failure'); error.cause = error;
  error.diagnostics = { big: 12n, bad: Infinity, fn: () => {}, array: new Float64Array([1, NaN]) };
  error.diagnostics.self = error.diagnostics;
  Object.defineProperty(error.diagnostics, 'badGetter', { enumerable: true, get() { throw new Error('getter'); } });
  error.diagnostics.proxy = new Proxy({}, { ownKeys() { throw new Error('enumeration'); } });
  const result = serializeDiagnosticError(error);
  assert.equal(result.cause, '[Circular cause]');
  assert.equal(result.diagnostics.self, '[Circular]');
  assert.equal(result.diagnostics.big, '12');
  assert.equal(result.diagnostics.bad, 'Infinity');
  assert.deepEqual(result.diagnostics.array, [1, 'NaN']);
  assert.equal(result.diagnostics.badGetter, '[Unreadable property]');
  assert.equal(result.diagnostics.proxy, '[Unreadable object]');
  assert.deepEqual(JSON.parse(JSON.stringify(result)), result);
});

test('cause depth and diagnostic arrays are explicitly bounded', () => {
  let error = new Error('innermost');
  for (let i = 0; i < 20; i++) error = new Error(`wrapper ${i}`, { cause: error });
  error.diagnostics = { history: Array.from({ length: 2000 }, (_, iteration) => ({ iteration })) };
  const result = serializeDiagnosticError(error);
  let end = result, depth = 0;
  while (end && typeof end === 'object') { end = end.cause; depth++; }
  assert.equal(depth, 8); assert.equal(end, '[Truncated: cause depth]');
  assert.equal(result.diagnostics.history.length, 1025);
  assert.equal(result.diagnostics.history.at(-1), '[Truncated: 976 array entries]');
  assert.deepEqual(JSON.parse(JSON.stringify(result)), result);
});
