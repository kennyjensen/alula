import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { streamtubeCornerConstraints } from '../src/euler/streamtube-corner-constraints.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { solveSparseDirect } from '../src/numerics/klu.js';
import { takeDoglegStep } from '../src/numerics/dogleg.js';

test('the retained default GUI state extends Cauchy along a feasible dogleg instead of a tiny projected Newton step', t => {
  const f = JSON.parse(readFileSync(new URL('./fixtures/streamtube-euler-weak-projection.json', import.meta.url)));
  const system = createStreamtubeBodySystem(f.input);
  const x = system.adoptGeometry(Float64Array.from(f.initialEuler.x), f.initialEuler.nodes), before = x.slice();
  const value = system.evaluate(x);
  assert.ok(Math.abs(value.diagnostics.residual - f.expected.residual) < 1e-10);
  const matrix = system.jacobian(x, { sparse: true });
  const linear = solveSparseDirect(matrix, value.residual.map(v => -v));
  assert.ok(linear.relativeResidual < 1e-10);
  const corners = streamtubeCornerConstraints(system, x).map(c => ({ value: c.value, gradient: c.gradient, lower: -.9 * c.value }));
  for (const radius of [8, .125]) {
    const step = takeDoglegStep({ initial: x, currentResidual: value.residual, matrix, newtonDirection: linear.x,
      residual: system.residual, admissible: system.admissible, radius, maxTrials: 20,
      linearizedConstraints: () => corners,
      constraintValues: state => streamtubeCornerConstraints(system, state, { derivatives: false }).map(c => c.value) });
    assert.equal(step.accepted, true, step.reason); assert.equal(step.kind, 'constrained-dogleg');
    assert.equal(step.projection.active, 0, 'the short Cauchy point is already linearly feasible');
    assert.ok(step.actualReduction > 3.8e-5);
    assert.ok(step.actualReduction > 2000 * f.expected.oldWeakReduction);
    assert.ok(step.predictedReduction > step.doglegSegment.cauchyPredictedReduction);
    assert.ok(step.doglegSegment.fraction > 0 && step.doglegSegment.fraction <= step.doglegSegment.maximumFraction);
    assert.ok(step.ratio > .99); assert.ok(step.scaledStepNorm <= radius);
    const next = system.evaluate(step.x);
    const quality = streamtubeMeshSnapshot({ system, nodes: next.nodes }).quality;
    assert.equal(quality.valid, true); assert.ok(next.diagnostics.residual < f.expected.residual);
    assert.ok(next.diagnostics.residual > 1e-10, 'this step does not establish convergence');
    assert.deepEqual(x, before);
    t.diagnostic(JSON.stringify({ radius, reduction: step.actualReduction, ratio: step.ratio,
      residual: next.diagnostics.residual, kind: step.kind, quality }));
  }
});
