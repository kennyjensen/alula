import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createStreamtubeBodySystem, solveStreamtubeBody } from '../src/euler/streamtube-body.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';

test('a retained two-element admissibility stall takes an inward Newton ray and reduces the unchanged Euler residual', t => {
  const fixture = JSON.parse(readFileSync(new URL('./fixtures/streamtube-admissibility-stall.json', import.meta.url)));
  const system = createStreamtubeBodySystem(fixture.input), { initialEuler } = fixture;
  const initial = system.adoptGeometry(Float64Array.from(initialEuler.x), initialEuler.nodes);
  const before = system.evaluate(initial), oldMesh = streamtubeMeshSnapshot({ system, nodes: before.nodes });
  assert.equal(oldMesh.quality.valid, true);
  assert.ok(oldMesh.quality.minCornerSine < 1e-9, 'The fixture must exercise an almost active geometry constraint');
  const result = solveStreamtubeBody(system, { initial, stepMethod: 'dogleg', initialTrustRadius: 1, maxIterations: 1 });
  const mesh = streamtubeMeshSnapshot({ system, nodes: result.nodes });
  assert.equal(result.converged, false, 'One admissible step does not certify a flow root');
  assert.equal(result.history.length, 2);
  assert.equal(result.history[1].stepKind, 'newton-ray');
  assert.ok(result.diagnostics.residual < before.diagnostics.residual);
  assert.ok(result.history[1].actualReduction > 0);
  assert.ok(result.history[1].reductionRatio > 1e-4);
  assert.ok(result.linearDiagnostics.maxRelativeResidual < 1e-10);
  assert.equal(mesh.quality.valid, true);
  assert.ok(mesh.quality.minCornerSine > 1e4 * oldMesh.quality.minCornerSine);
  t.diagnostic(JSON.stringify({ before: before.diagnostics.residual, after: result.diagnostics.residual,
    minimumCornerBefore: oldMesh.quality.minCornerSine, minimumCornerAfter: mesh.quality.minCornerSine,
    step: result.history[1] }));
});

test('the later default Euler stall escapes with a constrained Newton proposal when both ordinary directions fail', t => {
  const f = JSON.parse(readFileSync(new URL('./fixtures/streamtube-euler-projection-stall.json', import.meta.url)));
  const make = () => {
    const system = createStreamtubeBodySystem(f.input);
    return { system, initial: system.adoptGeometry(Float64Array.from(f.initialEuler.x), f.initialEuler.nodes) };
  };
  const a = make(), b = make(), controls = { stepMethod: 'dogleg', initialTrustRadius: 1, maxIterations: 1 };
  assert.ok(Math.abs(a.system.evaluate(a.initial).diagnostics.residual - f.expected.residual) < 1e-10);
  const plain = solveStreamtubeBody(a.system, { ...controls, initial: a.initial, projectedSteps: false });
  const projected = solveStreamtubeBody(b.system, { ...controls, initial: b.initial });
  assert.equal(plain.converged, false); assert.equal(plain.history.length, 1);
  assert.equal(projected.converged, false); assert.equal(projected.history.length, 2);
  const h = projected.history[1]; assert.equal(h.stepKind, 'projected-newton');
  assert.equal(h.projection.converged, true); assert.ok(h.projection.active > 0);
  assert.ok(h.actualReduction > 5e-5); assert.ok(h.reductionRatio > .9);
  assert.ok(projected.diagnostics.residual < .008); assert.ok(projected.linearDiagnostics.maxRelativeResidual < 1e-10);
  const quality = streamtubeMeshSnapshot({ system: b.system, nodes: projected.nodes }).quality;
  assert.equal(quality.valid, true); assert.ok(quality.minCornerSine > 1e-5);
  assert.equal(projected.surfaces.length, 4);
  t.diagnostic(JSON.stringify({ residual: projected.diagnostics.residual, quality, step: h }));
});
