// SPDX-License-Identifier: GPL-2.0-or-later
// Tiny exact-storage checks only: no mesh, equations, Jacobian or LU.
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { sparseMatrixFromRows } from '../src/numerics/sparse-rows.js';

// Bind the oracle to the archived implementation, independent of future
// integration edits to production sparse.js.
const archive = new URL('../docs/performance-assembly/before/sparse.js.txt', import.meta.url);
const original = fs.readFileSync(archive, 'utf8');
const imported = original.replace("from './linear.js'", `from '${new URL('../src/numerics/linear.js', import.meta.url).href}'`);
assert.notEqual(imported, original);
const legacy = await import(`data:text/javascript;base64,${Buffer.from(imported).toString('base64')}`);
const snapshot = rows => rows.map(row => Array.from(row));
const bytes = a => Buffer.from(a.buffer, a.byteOffset, a.byteLength);
function verify(rows) {
  const before = snapshot(rows), expected = legacy.sparseMatrix(rows.map(row => row.keys()));
  rows.forEach((row, i) => { for (const [col, value] of row) legacy.sparseAdd(expected, i, col, value); });
  const actual = sparseMatrixFromRows(rows);
  assert.equal(actual.n, expected.n);
  assert(actual.rowPtr instanceof Int32Array && actual.colIndex instanceof Int32Array && actual.values instanceof Float64Array);
  for (const key of ['rowPtr', 'colIndex', 'values']) assert.deepEqual(bytes(actual[key]), bytes(expected[key]), key);
  assert.deepEqual(snapshot(rows), before, 'Input Maps and insertion order remain unchanged.');
  return actual;
}

test('sorted CSR retains explicit zero entries and forces every missing diagonal', () => {
  const rows = [new Map([[3, -0], [2, 5]]), new Map(), new Map([[1, 0], [3, -2], [2, 7]]), new Map([[0, -4]])];
  const result = verify(rows);
  assert.deepEqual(Array.from(result.rowPtr), [0, 3, 4, 7, 9]);
  assert.deepEqual(Array.from(result.colIndex), [0, 2, 3, 1, 1, 2, 3, 0, 3]);
  assert.deepEqual(Array.from(result.values), [0, 5, 0, 0, 0, 7, -2, -4, 0]);
  result.values[1] = 100; assert.equal(rows[0].get(2), 5);
});

test('signed zero, subnormals and extreme finite FP64 values match archived 0+d semantics bitwise', () => {
  const pool = [-0, 0, Number.MIN_VALUE, -Number.MIN_VALUE, Number.MAX_VALUE, -Number.MAX_VALUE,
    1 + Number.EPSILON, 1 - Number.EPSILON / 2, Math.PI, -Math.E];
  const rows = pool.map((value, i) => new Map([[pool.length - 1 - i, value]]));
  const result = verify(rows);
  for (const value of result.values) if (value === 0) assert(Object.is(value, 0));
  const single = verify([new Map([[-0, -0]])]);
  assert.equal(single.colIndex[0], 0); assert(Object.is(single.values[0], 0));
});

test('preaccumulated cancellation values are serialized without re-summing contributions', () => {
  const row = new Map();
  for (const term of [1e16, 1, -1e16]) row.set(0, (row.get(0) ?? 0) + term);
  assert.equal(row.get(0), 0);
  const result = verify([row]);
  assert.equal(result.values[0], 0, 'Serialization retains the already rounded sum.');
});

test('deterministic small irregular patterns match all archived CSR bytes', () => {
  let state = 0x4a63f219;
  const random = () => { state ^= state << 13; state ^= state >>> 17; state ^= state << 5; return state >>> 0; };
  const pool = [0, -0, Number.MIN_VALUE, -Number.MIN_VALUE, 1e-150, -1e150, Math.PI, 1 + Number.EPSILON, -7];
  for (let example = 0; example < 200; example++) {
    const n = 1 + random() % 17, rows = Array.from({ length: n }, () => new Map());
    for (const row of rows) {
      const updates = random() % (2 * n + 1);
      for (let k = 0; k < updates; k++) row.set(random() % n, pool[random() % pool.length]);
    }
    verify(rows);
  }
});

test('invalid shape, square indices and nonfinite values are rejected without changing Maps', () => {
  for (const rows of [null, {}, new Map(), [], [null], [[]], Array(1), Array(0x80000000)])
    assert.throws(() => sparseMatrixFromRows(rows), /Sparse/);
  for (const col of [-1, 2, .5, '0', NaN, Infinity, 1n, Symbol('column')]) {
    const rows = [new Map([[col, 1]]), new Map()], before = snapshot(rows);
    assert.throws(() => sparseMatrixFromRows(rows), /column/); assert.deepEqual(snapshot(rows), before);
  }
  for (const value of [NaN, Infinity, -Infinity, '1', null, undefined, 1n, {}]) {
    const rows = [new Map([[0, value]])], before = snapshot(rows);
    assert.throws(() => sparseMatrixFromRows(rows), /value/); assert.deepEqual(snapshot(rows), before);
  }
});
