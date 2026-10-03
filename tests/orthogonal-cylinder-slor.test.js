import test from 'node:test';
import assert from 'node:assert/strict';
import { runOrthogonalCylinderControl } from '../scripts/validation/orthogonal-boundary-slor.js';

test('fixed exponential boundary control refines both harmonic mass and the independent wall-angle derivative', t => {
  // Rates are fixed from the coarsest grid, so refinement resolves the
  // same control-extension function rather than narrowing it with h.
  const first = Math.expm1(.5) / Math.expm1(3), last = 1 - Math.expm1(2.5) / Math.expm1(3);
  const decayRates = { lower: .45 / first, upper: .45 / last }, evidence = [];
  let previousMass = Infinity, previousAngle = Infinity;
  for (const nt of [6, 12, 24]) {
    const result = runOrthogonalCylinderControl({ nt, exponential: true, decayRates });
    assert.equal(result.converged, true, result.reason);
    assert.equal(result.quality.valid, true);
    assert.ok(result.massError < .3 * previousMass, `${previousMass} -> ${result.massError}`);
    assert.ok(result.maximumWallDerivativeShear < .3 * previousAngle, `${previousAngle} -> ${result.maximumWallDerivativeShear}`);
    assert.ok(result.wallShiftOverGap < .004 && result.maximumMovement < .004);
    assert.ok(result.history.at(-1).residual < 1e-10);
    previousMass = result.massError; previousAngle = result.maximumWallDerivativeShear;
    evidence.push({ nt, nx: result.nx, sweeps: result.sweeps, massError: result.massError, derivativeShear: result.maximumWallDerivativeShear,
      firstRowShear: result.wallShiftOverGap, maximumMovement: result.maximumMovement, minimumCornerSine: result.quality.minCornerSine,
      seconds: result.seconds, residual: result.history.at(-1).residual });
  }
  t.diagnostic(JSON.stringify({ decayRates, evidence }));
});
