import test from 'node:test';
import assert from 'node:assert/strict';
import { proposeLogarithmicShearUpdate as propose, proposeCoupledXfoilBLUpdate as nativeProposal } from '../src/euler/streamtube-coupled-xfoil-update.js';

// Real proposal and density limiter; controlled geometry adapter only.
// No global system, mesh, Euler solve or saved checkpoint is constructed.
function stub({ regime = 'wake', finite = true, mach = .185 } = {}) {
  const system = { ne: 1, n: 9, bl: { scale: .01, hasFiniteBase: finite,
    stations: [{ id: 0, kind: 'surface', regime: 'laminar' },
      { id: 1, kind: regime === 'wake' || regime === 'trailing-edge' ? 'wake' : 'surface', regime }],
    kernel: { parameters: { mach, gamma: 1.4 } }, thicknesses: x => Array.from(x),
    geometry: () => ({ coordinates: [{}, { wakeGap: regime === 'wake' || regime === 'trailing-edge' ? .002 : 0 }] }) },
  euler: { layout: { n: 1, densityCount: 1, independentWakeBanks: true,
    globals: { stagnation: [] }, positions: [] }, curves: [], decode: () => ({ nodes: [] }),
    setDisplacement: x => { system.lastDisplacement = x; } } };
  return system;
}
const state = (auxiliary = .03) => Float64Array.from([0, 2, 1, 2, 1, auxiliary, 1, 2, 1]);
const delta = () => new Float64Array(9);
const domainCode = 'COUPLED_XFOIL_CANDIDATE_DOMAIN';
function unchangedInput(x, d, beforeX, beforeD, system) {
  assert.deepEqual(x, beforeX); assert.deepEqual(d, beforeD);
  assert.deepEqual(system.lastDisplacement, Array.from(x.subarray(1)));
}

test('above-cap starting shear reaches the existing native cap even for a decreasing direction', () => {
  for (const regime of ['leading-transition', 'transition', 'turbulent', 'trailing-edge', 'wake'])
    for (const maximumStep of [1, .25]) {
      const s = stub({ regime }), x = state(.3), d = delta(); d[5] = -.01;
      const beforeX = x.slice(), beforeD = d.slice(), p = propose(s, x, d, { maximumStep });
      const unprojected = .3 * Math.exp(maximumStep * (-.01 / .3));
      // Reducing this step cannot get an above-cap starting state under .25.
      // The established native projection must remain available.
      assert(unprojected > .25); assert.equal(p.step, maximumStep); assert.equal(p.x[5], .25);
      assert.equal(p.projection.active, true); assert.equal(p.meritComparable, false);
      assert.equal(p.projection.equationsChanged, false);
      assert.deepEqual(p.projection.auxiliaryChanges, [{ ...s.bl.stations[1], before: unprojected,
        after: .25, correction: .25 - unprojected }]);
      assert.equal(p.projection.displacementChanges.length, 0);
      for (let k = 0; k < x.length; k++) if (k !== 5) assert.equal(p.x[k], x[k]);
      unchangedInput(x, d, beforeX, beforeD, s);
    }
});

test('below the cap only the shear trial map differs from the original native proposal', () => {
  for (const finite of [false, true]) for (const auxiliary of [.001, .03, .24]) for (const maximumStep of [1, .25]) {
    const x = state(auxiliary), d = Float64Array.from([.2, .3, -.01, .02, .01, auxiliary * .01, -.01, .02, .01]);
    const p = propose(stub({ finite }), x, d, { maximumStep });
    const native = nativeProposal(stub({ finite }), x, d, { maximumStep });
    assert.equal(p.step, native.step); assert.equal(p.viscousStep, native.viscousStep);
    assert.deepEqual(p.limiter, native.limiter); assert.deepEqual(p.undampedUpdate, native.undampedUpdate);
    assert.equal(p.x[5], auxiliary * Math.exp(p.step * (d[5] / auxiliary)));
    for (let k = 0; k < x.length; k++) if (k !== 5) assert.equal(p.x[k], native.x[k]);
    assert.equal(p.projection.active, false); assert.equal(p.meritComparable, undefined);
    assert.deepEqual(p.projection.auxiliaryChanges, []); assert.deepEqual(p.projection.displacementChanges, []);
  }
});

test('a synthetic newly mixed above-cap station is projected by the next proposal', () => {
  // Simulate only the output of an event; actual event selection is separate.
  const s = stub({ regime: 'laminar', finite: false }), x = state(9), d = delta();
  const laminar = propose(s, x, d); assert.equal(laminar.x[5], 9); assert.equal(laminar.projection.active, false);
  s.bl.stations[1].regime = 'transition'; x[5] = .29;
  const p = propose(s, x, d); assert.equal(p.x[5], .25); assert.equal(p.projection.active, true);
  assert.equal(p.projection.auxiliaryChanges[0].regime, 'transition');
  const next = propose(s, p.x, d); assert.deepEqual(next.x, p.x); assert.equal(next.projection.active, false);
});

test('laminar and similarity amplification retains its additive normalized limit', () => {
  for (const regime of ['laminar', 'similarity']) {
    const s = stub({ regime, finite: false }), x = state(9), d = delta(); d[5] = 40;
    const p = propose(s, x, d); assert.equal(p.step, .375); assert.equal(p.x[5], 24);
    assert.equal(p.viscousLimiter.variable, 'auxiliary'); assert.deepEqual(p.logarithmicShear, []);
    assert.deepEqual(p.projection.auxiliaryChanges, []); assert.equal(p.projection.active, false);
  }
});

test('underflow and overflow remain typed candidate failures without changing accepted inputs', () => {
  for (const direction of [-30, 30]) {
    const s = stub(), x = state(.03), d = delta(); d[5] = direction;
    const beforeX = x.slice(), beforeD = d.slice();
    assert.throws(() => propose(s, x, d), error => error.code === domainCode && error.recoverable === true
      && error.step === 1 && error.diagnostics.reason === 'logarithmic shear coordinate');
    unchangedInput(x, d, beforeX, beforeD, s);
  }
});

test('density, thickness and speed bounds plus DSLIM remain active with the shear cap', () => {
  const s = stub(), x = state(.3), d = delta();
  x[3] = 1.08; d[3] = -.09; d[5] = -.01; d[0] = .1;
  const p = propose(s, x, d);
  assert.equal(p.step, 1); assert.equal(p.x[0], Math.log1p(.1)); assert.equal(p.x[5], .25);
  assert.equal(p.projection.auxiliaryChanges.length, 1); assert.equal(p.projection.displacementChanges.length, 1);
  assert(p.projection.displacementChanges[0].projectedRawHk >= 1.02 - 1e-15);
  for (const [column, value, expectedStep, variable] of [[2, -2, .25, 'theta'], [3, -4, .25, 'delta-star'], [4, 1.5, .25, 'edge-speed']]) {
    const direction = delta(); direction[column] = value;
    const trial = propose(stub(), state(.3), direction);
    assert.equal(trial.step, expectedStep); assert.equal(trial.viscousLimiter.variable, variable); assert.equal(trial.x[5], .25);
  }
  const direction = delta(); direction[0] = -2;
  const density = propose(stub(), state(.3), direction);
  assert.equal(density.step, .25); assert.equal(density.limiter.kind, 'density-decrease');
});

test('physical proposal rejection and invalid accepted-state guards remain in force', () => {
  for (const [kind, expectedReason] of [['speed', 'nonpositive edge speed'], ['thermal', 'nonpositive thermal state'], ['geometry', 'candidate geometry']]) {
    const s = stub(), x = state(.3), d = delta(); d[5] = -.01;
    if (kind === 'speed') { x[4] = .1; d[4] = -.2; }
    if (kind === 'thermal') { s.bl.kernel.parameters.mach = .6; x[4] = 3.8; d[4] = .2; }
    if (kind === 'geometry') s.bl.geometry = () => { throw new Error('Collapsed BL wake interval.'); };
    const beforeX = x.slice(), beforeD = d.slice();
    assert.throws(() => propose(s, x, d), error => error.code === domainCode && error.recoverable === true
      && error.diagnostics.reason === expectedReason);
    unchangedInput(x, d, beforeX, beforeD, s);
  }
  const s = stub();
  assert.throws(() => propose(s, state(0), delta()), error => !error.recoverable && /accepted/.test(error.message));
  s.euler.layout.independentWakeBanks = false;
  assert.throws(() => propose(s, state(), delta()), /independently solved wake banks/);
});
