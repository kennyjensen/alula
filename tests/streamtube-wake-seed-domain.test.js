import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createHash } from 'node:crypto';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { initializeStreamtubeBLStation } from '../src/euler/streamtube-boundary-layers.js';

const mach = .2, reynolds = 1e6;
const upstream = { s: 1, aux: .03, theta: .0005, deltaStar: .0005065, ue: .8 };
function setup(gap = false) {
  const kernel = createIntegralKernel({ reynolds, mach, exactJacobian: true });
  const previous = gap ? { ...upstream, deltaStar: upstream.deltaStar + .0001, wakeGap: .0001 } : { ...upstream };
  const input = { upstream: previous, s: 1.01, ue: .98, regime: 'wake', reynolds, mach,
    ...(gap ? { wakeGap: .00008 } : {}), interval: kernel.interval, properties: s => kernel.station(s, 'wake') };
  return { kernel, input };
}

for (const gap of [false, true]) test(`accelerating ${gap ? 'finite-gap ' : ''}wake seed reaches the unchanged compressible BL equations`, () => {
  const { kernel, input } = setup(gap), before = structuredClone({ ...input, interval: undefined, properties: undefined });
  assert.ok(kernel.station(input.upstream, 'wake').rawHk > 1);
  const copied = { ...input.upstream, s: input.s, ue: input.ue,
    deltaStar: input.upstream.deltaStar - (input.upstream.wakeGap ?? 0) + (input.wakeGap ?? 0),
    ...(gap ? { wakeGap: input.wakeGap } : {}) };
  assert.ok(copied.deltaStar - (copied.wakeGap ?? 0) > copied.theta);
  assert.throws(() => kernel.interval({ upstream: input.upstream, downstream: copied, regime: 'wake' }), /raw kinematic shape/);
  const result = initializeStreamtubeBLStation(input);
  assert.equal(result.converged, true);
  if (result.mode === 'direct') assert.equal(result.state.ue, input.ue);
  else assert.ok(Math.abs(kernel.station(result.state, 'wake').rawHk - result.targetHK) < 1e-10);
  assert.ok(result.seedProjection.rawHk < 1);
  assert.ok(Math.abs(result.seedProjection.projectedRawHk - 1.00005) < 1e-14);
  assert.equal(result.seedProjection.wakeGap, input.wakeGap ?? 0);
  assert.equal(result.state.wakeGap ?? 0, input.wakeGap ?? 0);
  const block = kernel.interval({ upstream: input.upstream, downstream: result.state, regime: 'wake' });
  assert.ok(Math.max(...block.residual.map(Math.abs)) < 1e-10);
  assert.ok(kernel.station(result.state, 'wake').rawHk > 1);
  assert.deepEqual({ ...input, interval: undefined, properties: undefined }, before);
});

test('already-admissible copied wake seeds retain their exact existing result', () => {
  const { input } = setup(); input.ue = upstream.ue;
  const explicit = initializeStreamtubeBLStation(input), legacy = initializeStreamtubeBLStation({ ...input, mach: undefined });
  assert.equal(explicit.seedProjection, undefined); assert.deepEqual(explicit, legacy);
});

test('wake initialization does not repair negative thickness or hide equation errors', () => {
  const { input } = setup();
  assert.throws(() => initializeStreamtubeBLStation({ ...input, initialState: { theta: -.0005, deltaStar: -.0006 } }), /positive|Inadmissible/i);
  const error = new Error('Unexpected equation implementation failure');
  assert.throws(() => initializeStreamtubeBLStation({ ...input, interval: () => { throw error; } }), e => e === error);
});

test('native bounded inverse recovery solves the retained decelerating wake and matches original Fortran', () => {
  const fixture = JSON.parse(fs.readFileSync(new URL('./fixtures/fortran/wake-inverse-initialization.json', import.meta.url)));
  for (const [file, hash] of Object.entries(fixture.provenance.sourceHashes))
    assert.equal(createHash('sha256').update(fs.readFileSync(new URL(`../${file}`, import.meta.url))).digest('hex'), hash);
  const { parameters, ...input } = fixture.input, kernel = createIntegralKernel({ ...parameters, exactJacobian: true });
  for (const c of fixture.cases) {
    const block = kernel.interval({ upstream: input.upstream, downstream: c.state, regime: 'wake' });
    block.residual.forEach((v, i) => assert.ok(Math.abs(v - c.nativeResidual[i]) < 1e-12));
  }
  const result = initializeStreamtubeBLStation({ ...input, ...parameters,
    interval: kernel.interval, properties: state => kernel.station(state, 'wake') });
  assert.equal(result.converged, true); assert.equal(result.mode, 'inverse');
  assert.equal(result.inverseRecovery.method, 'xfoil-mrchue-bounded-inverse');
  assert.equal(result.inverseRecovery.reason, 'iteration limit');
  assert.ok(result.inverseRecovery.history.at(-1).residual > .004);
  assert.ok(result.history.some(h => h.residual > result.history[0].residual), 'Native initialization can cross the local merit barrier.');
  assert.ok(result.history.every(h => h.iteration === 0 || h.step * h.maximumRelativeIncrement <= .3 + 1e-15));
  const expected = fixture.cases.at(-1).state;
  for (const key of ['aux', 'theta', 'deltaStar', 'ue']) assert.ok(Math.abs(result.state[key] - expected[key]) < 1e-10);
  const block = kernel.interval({ upstream: input.upstream, downstream: result.state, regime: 'wake' });
  assert.ok(Math.max(...block.residual.map(Math.abs)) < 1e-10);
  assert.ok(Math.abs(kernel.station(result.state, 'wake').rawHk - fixture.targetHK) < 1e-10);
});
