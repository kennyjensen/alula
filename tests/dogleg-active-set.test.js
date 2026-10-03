import test from 'node:test';
import assert from 'node:assert/strict';
import { takeDoglegStep } from '../src/numerics/dogleg.js';

const options = () => ({ initial: Float64Array.of(0), matrix: Float64Array.of(1), currentResidual: Float64Array.of(1),
  newtonDirection: Float64Array.of(-1), radius: .1, maxTrials: 3 });

test('all rejected active-set proposals restore phase metadata and leave the accepted state untouched', () => {
  for (const failure of ['prepare', 'admissible', 'residual']) {
    let phase = 0, calls = 0;
    const input = options(), before = input.initial.slice();
    const result = takeDoglegStep({ ...input, activeSet: {
      snapshot: () => phase, restore: value => { phase = value; },
      prepare: candidate => {
        assert.equal(phase, 0, 'each proposal must start with the accepted active set');
        calls++; phase = 1; candidate[0] = 10;
        if (failure === 'prepare') throw new Error('transfer failed');
        return { changed: true, changes: [phase] };
      } }, admissible: () => failure !== 'admissible',
      residual: () => { throw new Error('residual failed'); } });
    assert.equal(result.accepted, false); assert.equal(calls, 3);
    assert.equal(phase, 0); assert.deepEqual(input.initial, before);
  }
});

test('an admissible active event retains its new phase without claiming an old-model merit reduction', () => {
  let phase = 0, calls = 0;
  const input = options(), before = input.initial.slice();
  const result = takeDoglegStep({ ...input, activeSet: {
    snapshot: () => phase, restore: value => { phase = value; },
    prepare: candidate => {
      assert.equal(phase, 0); calls++; phase = 1;
      // The first prepared event fails its independent physical check; the
      // second must start from the original state and active metadata.
      assert.equal(candidate[0], -.1);
      candidate[0] = 100; return { changed: true, changes: [{ from: 0, to: 1 }] };
    } }, admissible: () => calls > 1, residual: () => [100] });
  assert.equal(result.accepted, true); assert.equal(phase, 1); assert.equal(calls, 2);
  assert.equal(result.activeChange, true); assert.equal(result.meritComparable, false);
  for (const key of ['ratio', 'actualReduction', 'predictedReduction', 'linearizedResidualNorm']) assert.equal(result[key], null);
  assert.equal(result.radius, input.radius); assert.equal(result.residual[0], 100);
  assert.deepEqual(input.initial, before);
});

test('no active event preserves ordinary trust-region acceptance exactly', () => {
  const input = { ...options(), residual: ([x]) => [x + 1] };
  const baseline = takeDoglegStep(input);
  const observed = takeDoglegStep({ ...input, activeSet: { snapshot: () => 0, restore: () => {}, prepare: () => ({ changed: false }) } });
  assert.deepEqual(observed, baseline);
});
