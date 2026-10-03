import test from 'node:test';
import assert from 'node:assert/strict';
import { msesTemporaryMcrit } from '../src/euler/streamtube-shock-audit.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

test('MSES preferred dissipation formula, limiting values and monotonicity', () => {
  let previous = .99;
  for (const d of [0, 1e-8, 1e-4, .001, .01, .05, .15, .5, 1, 10, 100, Number.MAX_VALUE]) {
    const m = msesTemporaryMcrit(.99, d);
    assert(m >= .75 && m <= previous); previous = m;
    if (d <= 1) {
      const r = d ** 3 / (.15 * (d ** 2 + (.15 / 4) ** 2));
      assert(Math.abs(m - (.75 + .24 * Math.exp(-r * r))) < 3e-16);
    }
  }
  assert.equal(msesTemporaryMcrit(.99, 0), .99);
  for (const x of [NaN, Infinity, -1]) assert.throws(() => msesTemporaryMcrit(.99, x));
  for (const x of [NaN, .74, 1.01]) assert.throws(() => msesTemporaryMcrit(x, .1));
});

test('temporary equations preserve restart state and final acceptance uses prescribed Mcrit', () => {
  const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2 }), streamwiseMode: 'hybrid',
    hybrid: { ismom: 4, epsilonP: 1e-5 }, upwind: { mcrit: .99, mucon: 1, boundary: { kind: 'unfiltered-first-two' } } };
  const options = { edgeMatching: 'section-velocity', blThermodynamics: 'historical-common-isentrope',
    iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', maxIterations: 3, tolerance: 1e-10,
    dissipationEnhancement: true };
  const before = structuredClone(input), r = solveCoupledStreamtubeIses(input, options);
  assert.deepEqual(input, before);
  assert.equal(r.dissipationEnhancement.finalEquations, true);
  assert.equal(r.checkpoint.restart.input.upwind.mcrit, .99);
  assert(r.history.slice(1).some(h => h.dissipation.mcrit < .99), 'schedule must actually be exercised');
  const c = r.checkpoint, f = c.restart;
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  assert.deepEqual(system.evaluate(system.initial).families, r.families);
  // A parameter-only reconstruction must not move the grid or relabel transition.
  const temporary = createCoupledStreamtubeBody({ ...f.input, upwind: { ...f.input.upwind, mcrit: .8 } },
    { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const originalValue = system.evaluate(system.initial), temporaryValue = temporary.evaluate(temporary.initial);
  assert.deepEqual(temporary.initial, system.initial);
  assert.deepEqual(temporaryValue.outer.nodes, originalValue.outer.nodes);
  assert.deepEqual(temporaryValue.outer.undisplacedNodes, originalValue.outer.undisplacedNodes);
  assert.deepEqual(temporary.bl.snapshotActive(), system.bl.snapshotActive());
  assert.equal(r.converged, Math.max(...Object.values(r.families)) <= options.tolerance && r.mesh.quality.valid);
  const replay = solveCoupledStreamtubeIses(undefined, { resume: c, ...options, maxIterations: 0 });
  assert.deepEqual(replay.families, r.families);
  assert.deepEqual(replay.checkpoint.continuation.dissipationEnhancement, c.continuation.dissipationEnhancement);
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { resume: c, ...options, dissipationEnhancement: false }), /conflicting/);
  const broken = structuredClone(c); broken.continuation.dissipationEnhancement.previousDensityChange = NaN;
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { resume: broken, ...options }), /dissipation/);
});

test('periodic prescribed residual is an independent evaluation of the unchanged physical checkpoint', () => {
  const input = { ...intrinsicBodyFixture({ bodySegments: 4, tubes: 2 }), streamwiseMode: 'hybrid',
    hybrid: { ismom: 4, epsilonP: 1e-5 }, upwind: { mcrit: .99, mucon: 1, boundary: { kind: 'unfiltered-first-two' } } };
  let checked = 0;
  solveCoupledStreamtubeIses(input, { edgeMatching: 'section-velocity', blThermodynamics: 'historical-common-isentrope',
    iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', iterationRecovery: true,
    dissipationEnhancement: true, maxIterations: 6, tolerance: 1e-20,
    onCheckpoint: (cp, { history }) => {
      const h = history.at(-1);
      if (h.iteration !== 5) return;
      assert(h.prescribedResidual && !h.prescribedResidual.unavailable);
      const f = cp.restart, system = createCoupledStreamtubeBody({ ...f.input, upwind: { ...f.input.upwind, mcrit: .99 } },
        { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
      const value = system.evaluate(system.initial);
      for (const key of ['euler', 'boundaryLayer', 'edgeMatching']) assert.equal(h.prescribedResidual[key], value.families[key]);
      assert.deepEqual(system.initial, Float64Array.from([...f.initialEuler.x, ...f.initialBL]));
      checked++;
    } });
  assert.equal(checked, 1);
});
