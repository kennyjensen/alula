import test from 'node:test';
import assert from 'node:assert/strict';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
test('ISES checkpoints preserve the explicit historical hybrid model through zero-update resume', () => {
  const input = { ...intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  const options = { edgeMatching: 'section-velocity', blThermodynamics: 'historical-common-isentrope',
    maxIterations: 0, tolerance: 1e-10, stepAcceptance: 'admissible' };
  const first = solveCoupledStreamtubeIses(input, options);
  assert.ok(first.checkpoint, first.reason);
  assert.equal(first.checkpoint.restart.options.blThermodynamics, options.blThermodynamics);
  assert.deepEqual(first.checkpoint.restart.input.hybrid, input.hybrid);
  const saved = serialize(first.checkpoint), before = structuredClone(saved);
  const resumed = solveCoupledStreamtubeIses(undefined, {
    resume: saved, maxIterations: 0, tolerance: 1e-10, stepAcceptance: 'admissible' });
  assert.equal(resumed.initialRedistribution.resumed, true);
  assert.equal(resumed.linearDiagnostics.solves, 0);
  assert.deepEqual(resumed.families, saved.families);
  assert.deepEqual(serialize(resumed.flow.nodes), saved.restart.initialEuler.nodes);
  assert.equal(resumed.conditions.blThermodynamics, options.blThermodynamics);
  assert.equal(resumed.edgeThermodynamics.model, options.blThermodynamics);
  assert.deepEqual(resumed.checkpoint.restart.options, first.checkpoint.restart.options);
  assert.deepEqual(saved, before);
});

test('ISES driver forwards thermodynamic validation instead of silently dropping an unknown model', () => {
  const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2 }),
    streamwiseMode: 'hybrid', hybrid: { epsilonP: 1e-5 },
    upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } };
  assert.throws(() => solveCoupledStreamtubeIses(input, { blThermodynamics: 'unsupported', maxIterations: 0 }), /Unknown coupled BL thermodynamic model/);
});
