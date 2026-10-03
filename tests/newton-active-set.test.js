import test from 'node:test';
import assert from 'node:assert/strict';
import { solveNewton } from '../src/numerics/newton.js';

test('failed Newton phase transfers roll back metadata and never publish a rejected state', () => {
  for (const failure of ['prepare', 'admissible', 'residual']) {
    let phase = 0, calls = 0; const snapshots = [], initial = Float64Array.of(0);
    const r = solveNewton({ initial, maxIterations: 1, jacobian: () => Float64Array.of(1),
      residual: ([x]) => { if (phase) throw new Error('invalid phase residual'); return [x + 1]; },
      admissible: () => !phase || failure !== 'admissible', onState: s => snapshots.push(s.x),
      activeSet: { snapshot: () => phase, restore: p => { phase = p; }, prepare: candidate => {
        assert.equal(phase, 0); calls++; phase = 1; candidate[0] = 100;
        if (failure === 'prepare') throw new Error('invalid transfer'); return { changed: true };
      } } });
    assert.equal(r.converged, false); assert.equal(r.reason, 'line search failed');
    assert.equal(calls, 21); assert.equal(phase, 0); assert.equal(snapshots.length, 1);
    assert.deepEqual(r.x, initial); assert.deepEqual(initial, Float64Array.of(0));
  }
});

test('a bounded Newton phase event retains its new state without claiming old-phase merit descent', () => {
  let phase = 0;
  const r = solveNewton({ initial: [0], maxIterations: 1, residual: ([x]) => [phase ? x : x + 1],
    jacobian: () => Float64Array.of(1), activeSet: {
      snapshot: () => phase, restore: p => { phase = p; }, prepare: candidate => {
        phase = 1; candidate[0] = 100; return { changed: true, changes: [{ from: 0, to: 1 }] };
      } } });
  assert.equal(phase, 1); assert.equal(r.x[0], 100); assert.equal(r.converged, false);
  assert.equal(r.history[1].activeChange, true); assert.equal(r.history[1].meritComparable, false);
  assert.ok(r.history[1].residual > r.history[0].residual);
});

test('inactive event callbacks preserve the ordinary Newton solution and every accepted step', () => {
  const options = { initial: [3], residual: ([x]) => [x * x - 2], jacobian: ([x]) => Float64Array.of(2 * x) };
  const plain = solveNewton(options), observed = solveNewton({ ...options,
    activeSet: { snapshot: () => 0, restore: () => {}, prepare: () => ({ changed: false }) } });
  assert.deepEqual(observed, plain); assert.equal(observed.converged, true);
  assert.throws(() => solveNewton({ ...options, activeSet: {} }), /active-set/);
});
