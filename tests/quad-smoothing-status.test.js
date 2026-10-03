import test from 'node:test';
import assert from 'node:assert/strict';
import { quadSmoothingLabel } from '../src/ui/quad-smoothing-status.js';

test('partial accepted SLOR stays distinct from completed, in-progress and rejected smoothing', () => {
  const accepted = { attempted: true, converged: false, initialGuessAccepted: true };
  assert.equal(quadSmoothingLabel(accepted), 'SLOR incomplete · admissible initial grid retained');
  assert.equal(quadSmoothingLabel({ ...accepted, retainedOriginal: true }, { summary: true }), 'incomplete · admissible initial grid retained');
  assert.equal(quadSmoothingLabel({ attempted: true }), 'SLOR smoothing in progress');
  assert.equal(quadSmoothingLabel({ attempted: true, retainedOriginal: true }), 'SLOR rejected; original mesh retained');
  assert.equal(quadSmoothingLabel(undefined, { summary: true }), 'not run');
  assert.equal(quadSmoothingLabel({ converged: true, initialGuessAccepted: true }, { summary: true }), 'completed');
  assert.equal(quadSmoothingLabel({ converged: true, harmonicPassages: [1] }), 'SLOR smoothed · harmonic SLOR used');
});
