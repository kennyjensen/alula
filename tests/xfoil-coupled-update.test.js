import test from 'node:test';
import assert from 'node:assert/strict';
import { xfoilBLStepLimit, proposeCoupledXfoilBLUpdate } from '../scripts/validation/xfoil-coupled-update.js';

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
