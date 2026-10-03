import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';
import { streamtubeEdgeVelocity, streamtubeEdgeVelocityTangent } from '../src/euler/streamtube-edge-velocity.js';

test('thin-cell tangential slip and streamline rotation have the predicted quadratic width sensitivity', t => {
  // Exactly straight, parallel banks: the signed bend correction is zero.
  // The lower bank slips tangentially while both banks rotate. For affine
  // physical-node paths, W(t)=(w+t^2*k*tau)/sqrt(1+(t*k)^2), W'(0)=0.
  // Thus a correct first derivative does not bound a finite Newton step.
  const h = .0075, k = -.293, tau = .0019, observations = [];
  for (const w of [7.2e-6, 1.8e-6]) for (const angle of [0, .37]) {
    const transform = p => ({ x: Math.cos(angle) * p.x - Math.sin(angle) * p.y,
      y: Math.sin(angle) * p.x + Math.cos(angle) * p.y });
    const lower = [-h, 0, h].map(x => transform({ x, y: 0 }));
    const upper = [-h, 0, h].map(x => transform({ x, y: w }));
    const dl = [-h, 0, h].map(x => transform({ x: tau, y: k * x }));
    const du = [-h, 0, h].map(x => transform({ x: 0, y: k * x }));
    const input = { lower, upper, densities: [1, 1], massFlow: w, stagnationEnthalpy: 63, gamma: 1.4 };
    const linear = linearizeStreamtubeCell(input), tangent = linear.apply({ lower: dl, upper: du });
    assert.ok(Math.abs(streamtubeEdgeVelocityTangent(linear.value, tangent)) < 1e-9);
    for (const step of [.01, .005, .0025]) {
      const shift = (p, v) => ({ x: p.x + step * v.x, y: p.y + step * v.y });
      const cell = evaluateStreamtubeCell({ ...input, lower: lower.map((p, i) => shift(p, dl[i])), upper: upper.map((p, i) => shift(p, du[i])) });
      const expectedWidth = (w + step * step * k * tau) / Math.sqrt(1 + (step * k) ** 2);
      const expectedSpeed = w / expectedWidth;
      cell.geometry.normalAreas.forEach(width => assert.ok(Math.abs(width / expectedWidth - 1) < 1e-9));
      assert.ok(Math.abs(cell.geometry.pressureCurvature) < 1e-9);
      assert.ok(Math.abs(streamtubeEdgeVelocity(cell).ue - expectedSpeed) < 1e-9);
      assert.ok(expectedSpeed > 1);
      if (angle === 0) observations.push({ w, step, speedError: expectedSpeed - 1, quadraticCoefficient: (expectedSpeed - 1) / step ** 2 });
    }
  }
  const ratio = observations.at(-1).quadraticCoefficient / observations[2].quadraticCoefficient;
  assert.ok(ratio > 3.99 && ratio < 4.02);
  t.diagnostic(JSON.stringify({ observations, normalRefinementAmplification: ratio }));
});
