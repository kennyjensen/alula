import test from 'node:test';
import assert from 'node:assert/strict';
import { xfoilBLStepLimit, proposeCoupledXfoilBLUpdate, projectXfoilDisplacement, XFOIL_CANDIDATE_DOMAIN } from '../src/euler/streamtube-coupled-xfoil-update.js';

import { proposeCoupledXfoilBLUpdate as researchProposal } from '../scripts/validation/xfoil-coupled-update.js';
import { dslim } from '../src/viscous/xfoil/xbl.js';
import { hkin } from '../src/viscous/xfoil/xblsys.js';

function stub() {
  const system = { ne: 1, n: 9, bl: { scale: .01, hasFiniteBase: true,
    stations: [{ id: 0, kind: 'surface', regime: 'laminar' }, { id: 1, kind: 'wake', regime: 'wake' }],
    kernel: { parameters: { mach: .185, gamma: 1.4 } }, thicknesses: x => Array.from(x),
    geometry: () => ({ coordinates: [{}, { wakeGap: .002 }] }) }, euler: {
    layout: { n: 1, densityCount: 1, independentWakeBanks: true, globals: { stagnation: [] }, positions: [] }, curves: [],
    decode: () => ({ nodes: [] }), setDisplacement: x => { system.lastDisplacement = x; } } };
  return system;
}
const state = () => Float64Array.from([0, 2, 1, 2, 1, .2, 1, 2, 1]);

test('native UPDATE normalization and DHI/DLO select one scalar for each BL variable', t => {
  const s = stub(), x = state(); let cases = 0;
  for (const [column, change, expected, variable] of [
    [1, 40, .375, 'auxiliary'], [1, -20, .25, 'auxiliary'],
    [5, .8, .375, 'auxiliary'], [5, -.4, .25, 'auxiliary'],
    [2, 4, .375, 'theta'], [2, -2, .25, 'theta'],
    [3, 8, .375, 'delta-star'], [3, -4, .25, 'delta-star'],
    [4, 1.5, .25, 'edge-speed'], [4, -1.5, .25, 'edge-speed'],
  ]) {
    const d = new Float64Array(s.n); d[column] = change; const r = xfoilBLStepLimit(s, x, d);
    assert.equal(r.step, expected); assert.equal(r.limiter.variable, variable); cases++;
  }
  const d = new Float64Array(s.n); d[1] = 40; d[4] = 3;
  assert.equal(xfoilBLStepLimit(s, x, d).step, .125);
  t.diagnostic(JSON.stringify({ independentCases: cases, combinedGlobalCase: true }));
});

test('Euler density bounds and caller backtracking further reduce the single native BL step', () => {
  const s = stub(), x = state(), d = new Float64Array(s.n); d[0] = -2; d[1] = 40;
  const p = proposeCoupledXfoilBLUpdate(s, x, d);
  assert.equal(p.viscousStep, .375); assert.equal(p.step, .25); assert.equal(p.x[0], Math.log(.5));
  assert.equal(p.limiter.kind, 'density-decrease');
  const backtracked = proposeCoupledXfoilBLUpdate(s, x, d, { maximumStep: .125 });
  assert.equal(backtracked.step, .125); assert.equal(backtracked.x[1], x[1] + .125 * d[1]);
});

test('native Ctau cap and DSLIM preserve the other global-step fields and physical wake gap', () => {
  const s = stub(), x = state(), d = new Float64Array(s.n);
  x[3] = 1.08; d[3] = -.09; d[5] = .1; d[0] = .1;
  const before = Array.from(x), direction = Array.from(d), p = proposeCoupledXfoilBLUpdate(s, x, d);
  // Delta crosses compressible H>1 but native relative bounds allow a full
  // step. DSLIM handles it afterward; the Giles shape cap is absent.
  assert.equal(p.step, 1); assert.equal(p.viscousStep, 1);
  assert.equal(p.projection.auxiliaryChanges.length, 1); assert.equal(p.x[5], .25);
  assert.equal(p.projection.displacementChanges.length, 1); assert.ok(p.projection.displacementChanges[0].projectedRawHk >= 1.02 - 1e-15);
  for (const k of [1, 2, 4, 6, 8]) assert.equal(p.x[k], x[k] + d[k]);
  assert.equal(p.x[0], Math.log1p(.1)); assert.equal(p.meritComparable, false);
  assert.deepEqual(Array.from(x), before); assert.deepEqual(Array.from(d), direction);
  assert.deepEqual(s.lastDisplacement, before.slice(1));
});

test('inactive updates are identity and invalid accepted phases reject before use', () => {
  const s = stub(), x = state(), d = new Float64Array(s.n), p = proposeCoupledXfoilBLUpdate(s, x, d);
  assert.deepEqual(p.x, x); assert.equal(p.projection.active, false); assert.equal(p.meritComparable, undefined);
  x[5] = 0; assert.throws(() => proposeCoupledXfoilBLUpdate(s, x, d), /accepted XFOIL/);
});

test('DSLIM independently reproduces the native scalar on surface, wake and finite gap states', t => {
  let count = 0, maximumRelativeError = 0;
  for (const wake of [false, true]) for (const mach of [0, .185, .6]) for (const targetHk of [.99, 1.00001, 1.3]) {
    const theta = .0013, ue = 1.3, gamma = 1.4, wakeGap = wake ? .0024 : 0;
    const factor = mach ** 2 / (1 + .5 * (gamma - 1) * mach ** 2);
    const m2 = factor * ue ** 2 / (1 - .5 * (gamma - 1) * factor * ue ** 2);
    const deltaStar = wakeGap + theta * ((1 + .113 * m2) * targetHk + .29 * m2);
    const input = { theta, deltaStar, ue, wakeGap, wake, mach, gamma }, result = projectXfoilDisplacement(input);
    const native = dslim({ hkin }, deltaStar - wakeGap, theta, ue, m2, wake ? 1.00005 : 1.02) + wakeGap;
    const error = Math.abs(native - result.deltaStar) / native; maximumRelativeError = Math.max(maximumRelativeError, error);
    assert.ok(error < 2e-15); assert.ok(result.projectedRawHk >= result.minimumHk - 4e-15);
    if (targetHk === 1.3) { assert.equal(result.active, false); assert.equal(result.deltaStar, deltaStar); }
    assert.equal(result.massDefectChange, ue * (result.deltaStar - deltaStar));
    const scaled = projectXfoilDisplacement({ ...input, theta: theta * 1000, deltaStar: deltaStar * 1000, wakeGap: wakeGap * 1000 });
    assert.ok(Math.abs(scaled.deltaStar / 1000 - result.deltaStar) < 4e-18);
    const again = projectXfoilDisplacement({ ...input, deltaStar: result.deltaStar });
    assert.ok(Math.abs(again.deltaStar - result.deltaStar) < 4e-18); count++;
  }
  t.diagnostic(JSON.stringify({ cases: count, maximumRelativeError }));
});

test('valid runtime proposals exactly preserve the tested research operation and moving-gap restoration', () => {
  for (const maximumStep of [1, .25]) for (const finite of [false, true]) {
    const s = stub(), old = stub(), x = state(), d = new Float64Array(s.n);
    s.bl.hasFiniteBase = old.bl.hasFiniteBase = finite;
    x[3] = 1.08; x[7] = 1.2; d[3] = -.09; d[7] = -.2; d[5] = .1; d[0] = .1;
    // The proposed bank coordinates set the physical gap. It is not the
    // accepted-coordinate gap and is not scaled by the delta projection.
    s.bl.geometry = old.bl.geometry = z => ({ coordinates: [{}, { wakeGap: .002 + .001 * z[0] }] });
    const p = proposeCoupledXfoilBLUpdate(s, x, d, { maximumStep });
    assert.deepEqual(p, researchProposal(old, x, d, { maximumStep }));
    assert.deepEqual(s.lastDisplacement, Array.from(x.slice(s.ne)));
    if (finite && maximumStep === 1) {
      const wake = p.projection.displacementChanges.find(c => c.kind === 'wake');
      assert.equal(wake.wakeGap, (.002 + .001 * p.x[0]) / s.bl.scale);
      assert.equal(wake.massDefectChange, p.x[8] * (p.x[7] - (x[7] + d[7])));
    }
  }
});

test('only recoverable candidate speed/thermal/geometry domains carry a retry step; state and context survive', () => {
  const run = (kind, maximumStep = 1) => {
    const s = stub(), x = state(), d = new Float64Array(s.n);
    if (kind === 'speed') { x[4] = .1; d[4] = -.2; }
    if (kind === 'thermal') { s.bl.kernel.parameters.mach = .6; x[4] = 3.8; d[4] = .2; }
    if (kind === 'geometry') s.bl.geometry = () => { throw new Error('Collapsed BL wake interval.'); };
    const before = x.slice(), beforeDirection = d.slice(); let error;
    try { proposeCoupledXfoilBLUpdate(s, x, d, { maximumStep }); } catch (e) { error = e; }
    assert.deepEqual(x, before); assert.deepEqual(d, beforeDirection);
    assert.deepEqual(s.lastDisplacement, Array.from(x.slice(s.ne)));
    return error;
  };
  for (const [kind, reason] of [['speed', 'nonpositive edge speed'], ['thermal', 'nonpositive thermal state'], ['geometry', 'candidate geometry']]) {
    const error = run(kind); assert.equal(error.code, XFOIL_CANDIDATE_DOMAIN);
    assert.equal(error.recoverable, true); assert.equal(error.step, 1); assert.equal(error.diagnostics.reason, reason);
    if (kind === 'geometry') assert.equal(error.cause.message, 'Collapsed BL wake interval.');
    else assert.equal(run(kind, .25), undefined);
  }
  const s = stub(), x = state(), d = new Float64Array(s.n);
  const programming = new TypeError('missing geometry mapping'); s.bl.geometry = () => { throw programming; };
  assert.throws(() => proposeCoupledXfoilBLUpdate(s, x, d), e => e === programming && !e.recoverable);
  assert.deepEqual(s.lastDisplacement, Array.from(x.slice(s.ne)));
  x[4] = -1;
  assert.throws(() => proposeCoupledXfoilBLUpdate(s, x, d), e => !e.recoverable && /accepted/.test(e.message));
  x[4] = 100;
  assert.throws(() => proposeCoupledXfoilBLUpdate(s, x, d), e => !e.recoverable && /accepted.*thermal/.test(e.message));
  x[4] = 1;
  assert.throws(() => proposeCoupledXfoilBLUpdate(s, x, d, { maximumStep: 0 }), e => !e.recoverable && /maximum/.test(e.message));
});
