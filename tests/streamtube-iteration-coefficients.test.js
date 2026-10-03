import test from 'node:test';
import assert from 'node:assert/strict';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { streamtubeEquationControls } from '../src/euler/streamtube-equation-selection.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { capturePreparedStreamtubeAssembly } from '../src/euler/streamtube-prepared-assembly.js';
import { solveStreamtubeAssembly } from '../src/euler/streamtube-result.js';

for (const ismom of [undefined, 4]) test(`Euler ${ismom ?? 'isentropic'} reports every accepted step without changing the solve`, () => {
  const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 3, alpha: 2 }),
    flowModel: 'compressible', streamwiseMode: 'isentropic', normalStencil: 'body-stations',
    stagnationMotion: 'walls-only', geometryDomain: 'positive-simple', ...streamtubeEquationControls(ismom) };
  const caseData = { elements: input.bodies.map(b => ({ points: b.points })), mach: input.mach,
    alpha: input.alpha, referenceChord: 1.7, momentReference: { x: .3, y: .04 },
    gridIntervals: 4, gridTubes: 3, ...(ismom === undefined ? {} : { eulerIsmom: ismom }) };
  const system = createStreamtubeBodySystem(input), initial = system.initial, nodes = system.decode(initial).nodes;
  const preparedEuler = capturePreparedStreamtubeAssembly({ input, system, initial, nodes,
    initialEuler: { x: initial, nodes }, mesh: streamtubeMeshSnapshot({ system, initial, nodes }),
    diagnostics: {}, mach: input.mach, referenceChord: caseData.referenceChord, momentReference: caseData.momentReference }, caseData);
  const options = { preparedEuler, maxIterations: 2, adaptiveMcrit: false };
  const baseline = solveStreamtubeAssembly(caseData, options), frames = [], iterations = [];
  const result = solveStreamtubeAssembly(caseData, { ...options,
    onIteration: h => iterations.push(h.iteration),
    onMesh: mesh => { if (mesh.coefficientProgress) frames.push(structuredClone(mesh.coefficientProgress)); } });
  assert.deepEqual(result.flow.x, baseline.flow.x);
  assert.deepEqual(result.flow.nodes, baseline.flow.nodes);
  assert.deepEqual(result.flow.history, baseline.flow.history);
  assert.deepEqual(frames.map(f => f.iteration), iterations);
  assert.ok(frames.length > 1, result.diagnostics.reason);
  for (const frame of frames) {
    assert.equal(frame.kind, 'euler-pressure');
    assert.equal(frame.mach, caseData.mach);
    assert.equal(frame.actualAlpha, caseData.alpha);
    for (const key of ['cl', 'cd', 'cm']) assert.ok(Number.isFinite(frame.coefficients[key]), key);
    assert.equal(frame.coefficients.referenceChord, caseData.referenceChord);
    assert.deepEqual(frame.coefficients.momentReference, caseData.momentReference);
  }
  assert.deepEqual(frames.at(-1).coefficients, result.coefficients);
});
