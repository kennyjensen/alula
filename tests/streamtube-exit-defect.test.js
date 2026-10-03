// SPDX-License-Identifier: GPL-2.0-or-later
import test from 'node:test';
import assert from 'node:assert/strict';
import { streamtubeExitDefect } from '../src/euler/streamtube-exit-defect.js';

const close = (actual, expected, tolerance = 2e-12) => assert.ok(Number.isFinite(actual) && Number.isFinite(expected)
  && Math.abs(actual - expected) <= tolerance * Math.max(1, Math.abs(actual), Math.abs(expected)), `${actual} != ${expected}`);

function reference(mach = .65, gamma = 1.4) {
  const pressure = 1, density = 1, speed = mach * Math.sqrt(gamma * pressure / density);
  return { pressure, density, speed, stagnationEnthalpy: gamma / (gamma - 1) * pressure / density + .5 * speed ** 2 };
}

// Independent ideal-gas isentrope: prescribe a Mach number and recover
// static temperature, then use p/rho^gamma = constant. No solver imports.
function onIsentrope(source, mach, massFlow, gamma = 1.4) {
  const hSource = gamma / (gamma - 1) * source.pressure / source.density;
  const h = source.stagnationEnthalpy / (1 + .5 * (gamma - 1) * mach ** 2);
  const density = source.density * (h / hSource) ** (1 / (gamma - 1));
  return { pressure: source.pressure * (density / source.density) ** gamma, density,
    speed: mach * Math.sqrt((gamma - 1) * h), stagnationEnthalpy: source.stagnationEnthalpy, massFlow };
}

function shockFixture(massFlow = 2.1) {
  const freestream = reference(2);
  // Exact gamma=7/5, M1=2 Rankine-Hugoniot ratios: rho2/rho1=8/3,
  // p2/p1=9/2. Upstream mass, total enthalpy and normal momentum agree.
  const section = { pressure: 4.5, density: 8 / 3, speed: freestream.speed * 3 / 8,
    stagnationEnthalpy: freestream.stagnationEnthalpy, massFlow };
  return { freestream, sections: [section], referenceChord: .8 };
}

test('an isentropic exit gives zero defect at varied pressure and local subsonic/sonic/supersonic speeds', () => {
  const freestream = reference(), sections = [.25, .8, 1, 1.6, 2].map((mach, i) => onIsentrope(freestream, mach, .3 + i));
  const input = { freestream, sections, referenceChord: 1.7 }, before = structuredClone(input);
  const result = streamtubeExitDefect(input);
  close(result.momentumDeficit, 0); close(result.eulerWaveDragCoefficient, 0);
  result.sections.forEach((s, i) => {
    close(s.recovered.speed, freestream.speed); close(s.recovered.pressure, freestream.pressure);
    close(s.recovered.density, freestream.density); close(s.entropyOverR, 0);
    close(s.stagnationPressureRatio, 1); close(s.momentumDeficit, 0);
    close(s.recovered.area * s.recovered.density * s.recovered.speed, sections[i].massFlow);
    close(s.recovered.enthalpy + .5 * s.recovered.speed ** 2, freestream.stagnationEnthalpy);
  });
  assert.deepEqual(input, before);
  assert.equal(result.dragKind, 'euler-exit-entropy');
  assert.equal(Object.hasOwn(result, 'cd'), false);
  assert.equal(result.includesViscousWakeDefect, false);
  assert.equal(result.includesFinitePlanePressureThrust, false);
  assert.equal(result.physicalAcceptance, false);
});

test('normal-shock total-pressure loss gives the independently calculated positive entropy drag', () => {
  const input = shockFixture(), { freestream: free, sections: [exit], referenceChord: chord } = input;
  const result = streamtubeExitDefect(input), section = result.sections[0];
  // Closed-form normal-shock total-pressure ratio, independent of the
  // exit extrapolator. The unshocked reference has M_infinity=2 here.
  const totalPressureRatio = (8 / 3) ** 3.5 / 4.5 ** 2.5;
  const pTotalUpstream = 1.8 ** 3.5, pTotalDownstream = pTotalUpstream * totalPressureRatio;
  const qRecovered = Math.sqrt(2 * free.stagnationEnthalpy * (1 - (1 / pTotalDownstream) ** (2 / 7)));
  const expectedDeficit = exit.massFlow * (free.speed - qRecovered);
  close(totalPressureRatio, .7208738614847455);
  close(section.stagnationPressureRatio, totalPressureRatio);
  close(section.entropyOverR, -Math.log(totalPressureRatio));
  close(section.recovered.speed, qRecovered);
  close(result.momentumDeficit, expectedDeficit);
  close(result.eulerWaveDragCoefficient, expectedDeficit / (.5 * free.density * free.speed ** 2 * chord));
  assert.ok(result.momentumDeficit > 0);
  const gamma = 1.4;
  close(section.recovered.pressure / section.recovered.density ** gamma, exit.pressure / exit.density ** gamma);
  close(section.recovered.density, free.density * totalPressureRatio ** ((gamma - 1) / gamma));
  // The exit's actual stagnation density is retained, not reset to that
  // of the upstream reference isentrope.
  close(section.stagnationDensity / result.freestream.stagnationDensity, totalPressureRatio);
});

test('moving the exit along its own post-shock isentrope preserves the extrapolated result', () => {
  const input = shockFixture(), baseline = streamtubeExitDefect(input);
  for (const mach of [.3, .9, 1.3]) {
    const moved = onIsentrope(input.sections[0], mach, input.sections[0].massFlow);
    const result = streamtubeExitDefect({ ...input, sections: [moved] });
    close(result.momentumDeficit, baseline.momentumDeficit);
    close(result.sections[0].entropyOverR, baseline.sections[0].entropyOverR);
    close(result.sections[0].recovered.speed, baseline.sections[0].recovered.speed);
    close(result.sections[0].recovered.density, baseline.sections[0].recovered.density);
  }
});

test('physical mass weighting is invariant to partition and combines all inviscid tubes', () => {
  const input = shockFixture(), shocked = input.sections[0];
  const unshocked = onIsentrope(input.freestream, .8, .7);
  const baseline = streamtubeExitDefect({ ...input, sections: [shocked, unshocked] });
  const fractions = [.07, .18, .25, .5];
  const partitioned = fractions.flatMap(f => [shocked, unshocked].map(s => ({ ...s, massFlow: f * s.massFlow })));
  const result = streamtubeExitDefect({ ...input, sections: partitioned.reverse() });
  close(result.momentumDeficit, baseline.momentumDeficit);
  close(result.eulerWaveDragCoefficient, baseline.eulerWaveDragCoefficient);
  close(result.massFlow, shocked.massFlow + unshocked.massFlow);
  assert.equal(result.sections.length, 8);
});

test('negative entropy errors remain signed instead of being filtered or clipped', () => {
  const freestream = reference(), density = 1.01;
  const section = { pressure: freestream.pressure, density, stagnationEnthalpy: freestream.stagnationEnthalpy,
    speed: Math.sqrt(2 * (freestream.stagnationEnthalpy - 3.5 / density)), massFlow: .8 };
  const result = streamtubeExitDefect({ freestream, sections: [section] });
  assert.ok(result.sections[0].entropyOverR < 0);
  assert.ok(result.momentumDeficit < 0); assert.ok(result.eulerWaveDragCoefficient < 0);
  close(result.momentumDeficit, section.massFlow * (freestream.speed - section.speed));
  assert.equal(result.sections[0].recovered.speed, section.speed);
});

test('measured exit pressure and angle are diagnostics and do not insert a finite-plane force correction', () => {
  const input = shockFixture(), baseline = streamtubeExitDefect(input), angle = .23;
  const result = streamtubeExitDefect({ ...input, freestream: { ...input.freestream, direction: { x: 0, y: 2 } },
    sections: [{ ...input.sections[0], direction: { x: -3 * Math.sin(angle), y: 3 * Math.cos(angle) } }] });
  close(result.momentumDeficit, baseline.momentumDeficit);
  close(result.sections[0].exitAngleRadians, angle);
  close(result.maxAbsoluteExitAngleRadians, angle);
  close(result.maxRelativeExitPressureDeparture, 3.5);
  assert.equal(result.measuredExitDirections, 1);
  assert.equal(baseline.maxAbsoluteExitAngleRadians, null);
  assert.equal(baseline.measuredExitDirections, 0);
});

test('dimensional density/velocity/length scaling and explicit reference chord give consistent coefficients', () => {
  const input = shockFixture(), baseline = streamtubeExitDefect(input);
  const densityScale = 7, speedScale = 3, lengthScale = 5;
  const scale = s => ({ ...s, pressure: s.pressure * densityScale * speedScale ** 2,
    density: s.density * densityScale, speed: s.speed * speedScale,
    stagnationEnthalpy: s.stagnationEnthalpy * speedScale ** 2,
    ...(s.massFlow === undefined ? {} : { massFlow: s.massFlow * densityScale * speedScale * lengthScale }) });
  const result = streamtubeExitDefect({ freestream: scale(input.freestream), sections: input.sections.map(scale),
    referenceChord: input.referenceChord * lengthScale });
  close(result.eulerWaveDragCoefficient, baseline.eulerWaveDragCoefficient);
  close(result.momentumDeficit, baseline.momentumDeficit * densityScale * speedScale ** 2 * lengthScale);
  close(result.sections[0].recovered.area, baseline.sections[0].recovered.area * lengthScale);
  close(result.sections[0].entropyOverR, baseline.sections[0].entropyOverR);
  const changedChord = streamtubeExitDefect({ ...input, referenceChord: input.referenceChord * 4 });
  assert.equal(changedChord.momentumDeficit, baseline.momentumDeficit);
  close(changedChord.eulerWaveDragCoefficient, baseline.eulerWaveDragCoefficient / 4);
});

test('impossible thermal recovery, zero recovered mass capacity and incompatible gas are rejected', () => {
  const input = shockFixture(), freestream = reference();
  // Same h0 and local q, but p0_exit is too low to reach p_infinity with
  // positive kinetic energy; no inverse-isentrope branch can repair it.
  const impossible = { ...freestream, pressure: .01, density: .01, massFlow: 1 };
  assert.throws(() => streamtubeExitDefect({ freestream, sections: [impossible] }), /recover|extrapolat|capacity/i);
  // gamma=2 gives an exactly zero recovered speed: p0_exit=p_infinity.
  // Nonzero mass would require infinite recovered streamtube area.
  assert.throws(() => streamtubeExitDefect({ gamma: 2,
    freestream: { pressure: 1, density: 1, speed: 2, stagnationEnthalpy: 4 },
    sections: [{ pressure: .25, density: .25, speed: 2, stagnationEnthalpy: 4, massFlow: 1 }] }), /recover|capacity/i);
  for (const key of ['pressure', 'density', 'speed', 'stagnationEnthalpy', 'massFlow']) {
    for (const value of [0, -1, Infinity, NaN, undefined]) {
      const bad = structuredClone(input); bad.sections[0][key] = value;
      assert.throws(() => streamtubeExitDefect(bad), /positive|finite|enthalpy/i);
    }
  }
  for (const extra of [{ gamma: 1 }, { gamma: NaN }, { referenceChord: 0 }, { referenceChord: Infinity },
    { sections: [] }, { sections: null }, { freestream: null }])
    assert.throws(() => streamtubeExitDefect({ ...input, ...extra }));
  const wrongEnergy = structuredClone(input); wrongEnergy.sections[0].speed *= .9;
  assert.throws(() => streamtubeExitDefect(wrongEnergy), /consistent|enthalpy|energy/i);
  const differentH0 = structuredClone(input);
  differentH0.sections[0].stagnationEnthalpy += 1;
  differentH0.sections[0].speed = Math.sqrt(differentH0.sections[0].speed ** 2 + 2);
  assert.throws(() => streamtubeExitDefect(differentH0), /common|adiabatic|enthalpy/i);
  assert.throws(() => streamtubeExitDefect({ ...input, sections: [{ ...input.sections[0], direction: { x: 0, y: 1 } }] }), /direction/i);
  assert.throws(() => streamtubeExitDefect({ ...input, freestream: { ...input.freestream, direction: { x: 0, y: 0 } } }), /direction/i);
});
