import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubePhysicalSectionGas, streamtubeBLEdgeThermodynamics } from '../src/euler/streamtube-bl-edge-thermodynamics.js';

const gamma = 1.4, h0 = 3.5, rhoTotal = 1, pInf = 1;
const conditions = { gamma, h0, rhoTotal, pInf, flowModel: 'compressible' };
const gas = mach => {
  const enthalpy = h0 / (1 + .5 * (gamma - 1) * mach ** 2);
  const q = mach * Math.sqrt((gamma - 1) * enthalpy), rho = (enthalpy / h0) ** (1 / (gamma - 1));
  return { rho, q, p: (gamma - 1) / gamma * rho * enthalpy, enthalpy };
};
const near = (a, b, eps = 2e-14) => assert.ok(Math.abs(a - b) < eps, `${a} != ${b}`);

test('the common-isentrope diagnostic remains continuous through local sonic speed', () => {
  for (const mach of [.2, .999999, 1, 1.000001, 2]) {
    const r = streamtubePhysicalSectionGas(gas(mach), conditions);
    near(r.machSquared, mach * mach); near(r.entropyOverR, 0);
    near(r.stagnationPressureRatio, 1); near(r.commonIsentropeRelativeDensityError, 0);
  }
});

test('normal-shock entropy exposes the historical BL density approximation without changing physical gas', () => {
  const mach = 2, a = gas(mach);
  // Independent normal-shock relations, retaining upstream h0.
  const densityRatio = (gamma + 1) * mach ** 2 / ((gamma - 1) * mach ** 2 + 2);
  const pressureRatio = 1 + 2 * gamma / (gamma + 1) * (mach ** 2 - 1);
  const b = { rho: a.rho * densityRatio, q: a.q / densityRatio, p: a.p * pressureRatio };
  b.enthalpy = gamma / (gamma - 1) * b.p / b.rho;
  near(a.rho * a.q, b.rho * b.q);
  near(a.p + a.rho * a.q ** 2, b.p + b.rho * b.q ** 2);
  near(b.enthalpy + .5 * b.q ** 2, h0);
  const before = structuredClone(b), r = streamtubePhysicalSectionGas(b, conditions);
  const sigma = (Math.log(pressureRatio) - gamma * Math.log(densityRatio)) / (gamma - 1);
  near(r.entropyOverR, sigma); near(r.stagnationPressureRatio, Math.exp(-sigma));
  near(r.commonIsentropeRelativeDensityError, Math.exp(sigma) - 1);
  assert.ok(r.commonIsentropeRelativeDensityError > .38);
  assert.deepEqual(b, before);
});

test('wake diagnostics preserve unequal bank entropy and actual pressure independently of BL and transport speeds', () => {
  const upper = gas(1.2), lower = { ...gas(.8) };
  // A lower-total-pressure isentrope, with the same h0 but different bank gas.
  lower.rho *= .9; lower.p *= .9;
  const make = s => ({ states: [{ ...s }, { ...s }], interfacePressure: { lower: s.p, upper: s.p } });
  const flow = { cells: [[[make(lower)], [make(upper)]]], transportSpeeds: [[99, 88]] };
  const euler = { conditions, layout: { nx: 2, tubes: [1, 1] } };
  const bl = { stations: [{ id: 0, kind: 'wake', body: 0, i: 2 }] };
  const states = [{ ue: .7 }], before = structuredClone({ flow, states });
  const d = streamtubeBLEdgeThermodynamics({ flow, euler, bl, states }), s = d.stations[0];
  assert.equal(s.banks.length, 2); assert.deepEqual(s.banks.map(b => b.side), ['upper', 'lower']);
  near(s.banks[0].interfacePressure, upper.p); near(s.banks[1].interfacePressure, lower.p);
  near(d.maxWakeBankEntropyDifference, -Math.log(.9)); near(d.maxDensityMismatch, 1 / .9 - 1);
  assert.equal(s.boundaryLayer.ue, .7);
  assert.deepEqual(s.banks[0].sections.map(g => g.q), [upper.q, upper.q]);
  assert.deepEqual({ flow, states }, before);
});

test('invalid thermal and nonfinite diagnostic states are rejected', () => {
  assert.throws(() => streamtubePhysicalSectionGas({ ...gas(1), enthalpy: 0 }, conditions), /positive/);
  assert.throws(() => streamtubePhysicalSectionGas({ ...gas(1), rho: NaN }, conditions), /positive/);
});
