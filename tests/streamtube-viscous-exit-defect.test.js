// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeViscousExitDefect } from '../src/euler/streamtube-viscous-exit-defect.js';

const close = (a, b, tolerance = 3e-13) => assert.ok(Number.isFinite(a) && Number.isFinite(b)
  && Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);
const fixture = () => ({ gamma: 1.4, referenceChord: .8, gasSource: 'measured',
  freestream: { density: 1.2, speed: 90, pressure: 101000 },
  wakes: [{ density: 1.17, speed: 91, pressure: 101000, theta: .0012, deltaStar: .0018 }] });

test('equal exit and freestream pressure gives the exact measured wake momentum defect without mutation', () => {
  const input = fixture(), before = structuredClone(input), out = streamtubeViscousExitDefect(input);
  const w = input.wakes[0], expected = w.density * w.speed ** 2 * w.theta;
  assert.equal(out.wakes[0].correctionFactor, 1); assert.equal(out.pressureCorrection, 0);
  assert.equal(out.momentumDeficit, out.exitMomentumDeficit);
  close(out.momentumDeficit, expected);
  close(out.viscousDragCoefficient, 2 * expected / (1.2 * 90 ** 2 * .8));
  assert.deepEqual(input, before);
  out.wakes[0].exit.theta = 200; out.freestream.speed = 1;
  assert.deepEqual(input, before);
  assert.equal(out.includesEulerWaveDefect, false); assert.equal(out.includesSolidPressureIntegral, false);
  assert.equal(Object.hasOwn(out, 'cd'), false); assert.equal(out.physicalAcceptance, false);
});

test('prescribed pressure power agrees with an independent analytic case in both pressure directions', () => {
  // gamma=2, M_inf=1, H_exit=1: H_inf=2, H_avg=3/2,
  // exponent=3/4. Pressure recovery by 16 gives the exact factor8.
  const input = { gamma: 2, referenceChord: 2, gasSource: 'physical-euler',
    freestream: { density: 1, speed: Math.sqrt(2), pressure: 1 },
    wakes: [{ density: 3, speed: 2, pressure: 1 / 16, theta: .01, deltaStar: .01 }] };
  let out = streamtubeViscousExitDefect(input);
  close(out.freestream.machSquared, 1); close(out.wakes[0].averageShape, 1.5);
  close(out.wakes[0].pressureExponent, .75); close(out.wakes[0].correctionFactor, 8);
  close(out.momentumDeficit, .96); close(out.pressureCorrection, .84);
  close(out.viscousDragCoefficient, .48);
  input.wakes[0].pressure = 16;
  out = streamtubeViscousExitDefect(input);
  close(out.wakes[0].correctionFactor, 1 / 8); close(out.momentumDeficit, .015);
  close(out.pressureCorrection, -.105); assert.ok(out.pressureCorrection < 0);
});

test('independent element wakes add, and subdividing one thickness profile preserves the total', () => {
  const input = fixture();
  const second = { density: 1.13, speed: 88, pressure: 100700, theta: .0021, deltaStar: .0033 };
  const first = streamtubeViscousExitDefect(input);
  const other = streamtubeViscousExitDefect({ ...input, wakes: [second] });
  const sum = streamtubeViscousExitDefect({ ...input, wakes: [...input.wakes, second] });
  close(sum.momentumDeficit, first.momentumDeficit + other.momentumDeficit);
  close(sum.viscousDragCoefficient, first.viscousDragCoefficient + other.viscousDragCoefficient);
  const parts = [.125, .375, .5].map(f => ({ ...second, theta: f * second.theta, deltaStar: f * second.deltaStar }));
  const divided = streamtubeViscousExitDefect({ ...input, wakes: parts });
  close(divided.momentumDeficit, other.momentumDeficit);
  close(divided.viscousDragCoefficient, other.viscousDragCoefficient);
});

test('dimensional density, speed and length changes preserve coefficients and reference chord scales inversely', () => {
  const input = fixture(); input.wakes[0].pressure *= .999;
  const baseline = streamtubeViscousExitDefect(input), R = 7, U = 3, L = 5;
  const scale = s => ({ density: R * s.density, speed: U * s.speed, pressure: R * U ** 2 * s.pressure,
    ...(s.theta === undefined ? {} : { theta: L * s.theta, deltaStar: L * s.deltaStar }) });
  const out = streamtubeViscousExitDefect({ ...input, referenceChord: L * input.referenceChord,
    freestream: scale(input.freestream), wakes: input.wakes.map(scale) });
  close(out.viscousDragCoefficient, baseline.viscousDragCoefficient);
  close(out.momentumDeficit, baseline.momentumDeficit * R * U ** 2 * L);
  close(out.freestream.machSquared, baseline.freestream.machSquared);
  close(out.wakes[0].pressureExponent, baseline.wakes[0].pressureExponent);
  const chord = streamtubeViscousExitDefect({ ...input, referenceChord: 4 * input.referenceChord });
  assert.equal(chord.momentumDeficit, baseline.momentumDeficit);
  close(chord.viscousDragCoefficient, baseline.viscousDragCoefficient / 4);
});

test('small signed pressure changes retain their logarithmic correction', () => {
  for (const sign of [-1, 1]) {
    const input = fixture(); input.freestream.pressure = 1;
    input.wakes[0].pressure = 1 + sign * 2 ** -42;
    const out = streamtubeViscousExitDefect(input), wake = out.wakes[0];
    const k = wake.averageShape / (input.gamma * out.freestream.machSquared);
    const expectedLog = -Math.log1p(sign * 2 ** -42);
    close(wake.logPressureRecovery, expectedLog, 1e-27);
    const expectedCorrection = wake.exitMomentumDeficit * k * expectedLog;
    assert.ok(Math.abs(wake.pressureCorrection / expectedCorrection - 1) < 1e-12);
    assert.equal(Math.sign(wake.pressureCorrection), -sign);
  }
});

test('constant-H-average and constant-Mach wake momentum ODE independently integrates to the reported defect', () => {
  const input = fixture(); input.wakes[0].pressure = 100500;
  const out = streamtubeViscousExitDefect(input), wake = input.wakes[0];
  const m2 = input.freestream.density * input.freestream.speed ** 2 / (input.gamma * input.freestream.pressure);
  const H = (wake.deltaStar / wake.theta + 1 + (input.gamma - 1) * m2) / 2;
  // RK4 in physical pressure, integrating dD/dp=H*D/(gamma*M_inf²*p).
  const steps = 128, dp = (input.freestream.pressure - wake.pressure) / steps;
  let D = wake.density * wake.speed ** 2 * wake.theta, p = wake.pressure;
  const rhs = (pressure, defect) => H * defect / (input.gamma * m2 * pressure);
  for (let i = 0; i < steps; i++) {
    const a = rhs(p, D), b = rhs(p + dp / 2, D + dp * a / 2);
    const c = rhs(p + dp / 2, D + dp * b / 2), d = rhs(p + dp, D + dp * c);
    D += dp * (a + 2 * b + 2 * c + d) / 6; p += dp;
  }
  close(out.momentumDeficit, D);
});

test('explicit gas choice is retained without replacing supplied density or pressure', () => {
  const input = fixture(), baseline = streamtubeViscousExitDefect(input);
  for (const gasSource of ['historical-common-isentrope', 'physical-euler', 'measured']) {
    const out = streamtubeViscousExitDefect({ ...input, gasSource });
    assert.equal(out.gasSource, gasSource); assert.equal(out.wakes[0].gasSource, gasSource);
    assert.equal(out.wakes[0].exit.density, input.wakes[0].density);
    assert.equal(out.wakes[0].exit.pressure, input.wakes[0].pressure);
    assert.equal(out.momentumDeficit, baseline.momentumDeficit);
  }
});

test('invalid states and unrepresentable derived corrections fail explicitly without clipping', () => {
  const input = fixture();
  for (const key of ['density', 'speed', 'pressure', 'theta', 'deltaStar']) for (const value of [0, -1, Infinity, NaN, undefined]) {
    const bad = structuredClone(input); bad.wakes[0][key] = value;
    assert.throws(() => streamtubeViscousExitDefect(bad), /positive|finite/);
  }
  for (const key of ['density', 'speed', 'pressure']) for (const value of [0, -1, Infinity, NaN, undefined]) {
    const bad = structuredClone(input); bad.freestream[key] = value;
    assert.throws(() => streamtubeViscousExitDefect(bad), /positive|finite/);
  }
  for (const extra of [{ gamma: 1 }, { gamma: NaN }, { referenceChord: 0 }, { referenceChord: Infinity },
    { wakes: [] }, { wakes: null }, { freestream: null }, { gasSource: undefined }, { gasSource: 'mixed' }])
    assert.throws(() => streamtubeViscousExitDefect({ ...input, ...extra }));
  for (const pressure of [1e-300, 1e300]) {
    const bad = structuredClone(input); bad.wakes[0].pressure = pressure;
    assert.throws(() => streamtubeViscousExitDefect(bad), /positive|finite/);
  }
  const overflow = structuredClone(input); overflow.freestream.speed = 1e300;
  assert.throws(() => streamtubeViscousExitDefect(overflow), /positive|finite/);
});
