import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { solveCoupledStreamtubeIses } from '../src/euler/streamtube-coupled-ises.js';

const plain = value => JSON.parse(JSON.stringify(value, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x));
const input = { ...intrinsicBodyFixture({ elements: 1, bodySegments: 4, tubes: 2, contourPanels: 40 }),
  wakeGeometry: 'independent-banks', wakeOutlet: 'banks', wakeDisplacementMotion: 'te-center' };
const options = { transitionMode: 'automatic', tripFractions: [[1, 1]], ncrit: 100, edgeMatching: 'section-velocity',
  blUpdate: 'xfoil', projectionGeometry: 'boundary-increment', iterationGeometry: 'convex', stepAcceptance: 'admissible',
  stagnationLimiter: 'listing', maxIterations: 1, tolerance: 1e-10 };
let ordinary;

test('private opt-in with no raw event retains exact real state, residual, mesh and accepted step', () => {
  ordinary = solveCoupledStreamtubeIses(input, options);
  const r = solveCoupledStreamtubeIses(input, { ...options, eventProfile: 'xfoil-mrchdu' });
  assert.equal(ordinary.x.length, 169); assert.equal(r.history.length, 2, r.reason);
  assert.equal(r.history[1].eventProfile.active, false);
  for (const key of ['x', 'residual', 'families', 'flow', 'boundaryLayer', 'solverInput', 'coupledOptions'])
    assert.deepEqual(plain(r[key]), plain(ordinary[key]), key);
  assert.deepEqual(r.mesh.nodes, ordinary.mesh.nodes); assert.deepEqual(r.mesh.quality, ordinary.mesh.quality);
  assert.equal(r.history[1].step, ordinary.history[1].step); assert.equal(r.history[1].backtracks, ordinary.history[1].backtracks);
  assert.equal(Object.hasOwn(ordinary, 'eventProfile'), false);
  assert.equal(Object.hasOwn(ordinary.checkpoint.continuation, 'eventProfile'), false);
  assert.equal(r.checkpoint.continuation.eventProfile, 'xfoil-mrchdu');
});

test('explicitly upgrading a checkpoint policy preserves zero-step physical replay and subsequent resume policy', () => {
  const old = plain(ordinary.checkpoint), upgraded = structuredClone(old);
  upgraded.continuation.eventProfile = 'xfoil-mrchdu';
  const { iterationGeometry, stepAcceptance, stagnationLimiter } = old.continuation;
  const r = solveCoupledStreamtubeIses(undefined, { resume: upgraded, maxIterations: 0,
    iterationGeometry, stepAcceptance, stagnationLimiter });
  for (const key of ['x', 'residual', 'families', 'flow', 'boundaryLayer']) assert.deepEqual(plain(r[key]), plain(ordinary[key]), key);
  assert.deepEqual(r.checkpoint.continuation, upgraded.continuation);
  assert.deepEqual(old, plain(ordinary.checkpoint));
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { resume: old, eventProfile: 'xfoil-mrchdu', maxIterations: 0,
    iterationGeometry, stepAcceptance, stagnationLimiter }), /event profile controls/);
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { resume: upgraded, eventProfile: 'none', maxIterations: 0,
    iterationGeometry, stepAcceptance, stagnationLimiter }), /event profile controls/);
  const malformed = structuredClone(old); malformed.continuation.eventProfile = null;
  assert.throws(() => solveCoupledStreamtubeIses(undefined, { resume: malformed, maxIterations: 0,
    iterationGeometry, stepAcceptance, stagnationLimiter }), /iteration controls/);
});

test('malformed and unsupported private policies reject before any checkpoint is published', () => {
  let published = 0;
  for (const eventProfile of [null, false, 'mset']) assert.throws(() => solveCoupledStreamtubeIses(input,
    { ...options, maxIterations: 0, eventProfile, onCheckpoint: () => published++ }), /iteration controls/);
  for (const patch of [{ transitionMode: 'fixed-trip' }, { blUpdate: 'giles' }, { projectionGeometry: 'fixed' }])
    assert.throws(() => solveCoupledStreamtubeIses(input, { ...options, ...patch, maxIterations: 0,
      eventProfile: 'xfoil-mrchdu', onCheckpoint: () => published++ }), /iteration controls/);
  assert.equal(published, 0);
});

test('a later rejected grid trial retains completed native profile warnings and changes', async () => {
  const base = new URL('../src/euler/streamtube-coupled-ises.js', import.meta.url);
  let gridChecks = 0;
  const diagnostics = { active: true, initialGuessOnly: true, prediction: { bodies: [{ localConvergenceWarnings: ['native local warning'] }] } };
  globalThis.__eventProfileRejectionTest = {
    prepareCoupledTransitionProfileTrial: (system, x, { decoded, baseNodes }) => ({ x, nodes: baseNodes, decoded, diagnostics }),
    prepareConvexWakeGridUpdate: (_, nodes) => {
      if (!gridChecks++) throw new Error('deliberate later grid rejection');
      return { nodes };
    },
  };
  const source = fs.readFileSync(base, 'utf8').replace(/^import \{ (prepareCoupledTransitionProfileTrial|prepareConvexWakeGridUpdate) \} from '[^']+';/gm,
    (_, name) => `const { ${name} } = globalThis.__eventProfileRejectionTest;`)
    .replace(/from\s+(['"])(\.[^'"]+)\1/g, (_, quote, path) => `from ${quote}${new URL(path, base).href}${quote}`);
  const module = await import('data:text/javascript;base64,' + Buffer.from(source + '\n//# sourceURL=event-profile-rejection-test.js').toString('base64'));
  delete globalThis.__eventProfileRejectionTest;
  const r = module.solveCoupledStreamtubeIses(input, { ...options, eventProfile: 'xfoil-mrchdu' });
  assert.equal(r.history.length, 2, r.reason); assert.equal(r.history[1].backtracks, 1);
  assert.equal(r.history[1].rejections[0].stage, 'Newton grid step');
  assert.deepEqual(r.history[1].rejections[0].eventProfile, diagnostics);
  assert.equal(r.mesh.quality.valid, true);
});
