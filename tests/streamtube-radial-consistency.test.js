import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';

const close = (actual, expected, tolerance = 3e-13) => assert.ok(
  Math.abs(actual - expected) <= tolerance * Math.max(1, Math.abs(actual), Math.abs(expected)),
  `${actual} != ${expected}`);

// Independent smooth 2-D Euler solution, without production gas helpers.
// At r=1, rho=1 and sound speed=1, so q=Mach and p=1/gamma.
// rho(q)*q*r=K increases strictly on 0<q<q_sonic; bisection selects
// its unique subsonic root. Reversing the radial direction gives a sink.
function radialSolution(gamma, referenceMach) {
  const hReference = 1 / (gamma - 1), h0 = hReference + referenceMach ** 2 / 2;
  const K = referenceMach, pressureReference = 1 / gamma;
  const sonicSpeed = Math.sqrt(2 * (gamma - 1) * h0 / (gamma + 1));
  const density = q => ((h0 - q * q / 2) / hReference) ** (1 / (gamma - 1));
  const at = r => {
    assert.ok(r > 0 && density(sonicSpeed) * sonicSpeed * r > K, 'Subsonic radial branch must exist.');
    let left = 0, right = sonicSpeed;
    for (let iteration = 0; iteration < 64; iteration++) {
      const middle = (left + right) / 2;
      if (density(middle) * middle * r < K) left = middle;
      else right = middle;
    }
    const q = (left + right) / 2, rho = density(q), p = pressureReference * rho ** gamma;
    return { q, rho, p, machSquared: q * q * rho / (gamma * p) };
  };
  return { at, K, h0, pressureReference };
}

function sampleRadialCell({ gamma, referenceMach, direction, h }) {
  const exact = radialSolution(gamma, referenceMach), halfAngle = h / 2;
  // Both the radial interval and angular width halve at each refinement.
  const radii = [1 - direction * h, 1, 1 + direction * h];
  const point = (r, angle) => ({ x: r * Math.cos(angle), y: r * Math.sin(angle) });
  // For inward flow, exchange the banks as well as reversing radial order
  // so all signed quadrilateral areas stay positive.
  const lower = radii.map(r => point(r, -direction * halfAngle));
  const upper = radii.map(r => point(r, direction * halfAngle));
  const sectionRadii = [0, 1].map(i => (radii[i] + radii[i + 1]) * Math.cos(halfAngle) / 2);
  const sampled = sectionRadii.map(exact.at);
  // This is the exact flux through either the arc OR its straight chord:
  // integral K*R/(R^2+y^2) dy = K*deltaAngle, where R is the chord midpoint.
  const massFlow = exact.K * (2 * halfAngle);
  const cell = evaluateStreamtubeCell({ lower, upper, densities: sampled.map(s => s.rho),
    massFlow, stagnationEnthalpy: exact.h0, gamma });
  const chordFactor = halfAngle / Math.tan(halfAngle);
  const sectionDistance = Math.abs(sectionRadii[1] - sectionRadii[0]);
  for (let i = 0; i < 2; i++) {
    const state = cell.states[i], reference = sampled[i];
    // Do not accidentally replace actual chord quadrature by arc sampling.
    close(state.q, reference.q * chordFactor);
    close(reference.rho * reference.q * sectionRadii[i], exact.K);
    close(gamma / (gamma - 1) * reference.p / reference.rho + reference.q ** 2 / 2, exact.h0);
    close(reference.p / reference.rho ** gamma, exact.pressureReference);
    assert.ok(state.machSquared > 0 && state.machSquared < 1);
  }
  assert.ok(direction * (sampled[1].q - sampled[0].q) < 0, 'Sink accelerates and source decelerates.');
  assert.ok(Math.abs(sampled[1].q - sampled[0].q) > .001, 'Nonconstant streamwise speed must be resolved.');
  close(cell.pressureCorrection, 0);
  const centerPressure = exact.at(1).p;
  return { h,
    // The production residual is a pressure difference, not a derivative.
    // Divide by ds/r_ref as well as fixed p_ref (r_ref=1) to test consistency.
    normalizedMomentumError: Math.abs(cell.streamwiseResidual) / (exact.pressureReference * sectionDistance),
    sectionPressureError: Math.max(...cell.states.map((s, i) => Math.abs(s.p - sampled[i].p))) / exact.pressureReference,
    // Recovered interface pressure represents the radial side's midpoint
    // r=1; midpoint versus side-average quadrature differs by O(h^2).
    sidePressureError: Math.max(...Object.values(cell.interfacePressure).map(p => Math.abs(p - centerPressure))) / exact.pressureReference,
    rawMomentumResidual: cell.streamwiseResidual,
    entropyJump: cell.entropyJump,
    chordSpeedFactor: chordFactor,
  };
}

test('smooth subsonic radial sink/source momentum and pressure are second-order consistent', t => {
  const fields = ['normalizedMomentumError', 'sectionPressureError', 'sidePressureError'];
  for (const gamma of [1.4, 5 / 3]) for (const referenceMach of [.3, .6]) for (const direction of [-1, 1]) {
    const label = { gamma, referenceMach, flow: direction === -1 ? 'accelerating-sink' : 'decelerating-source' };
    const samples = [.1, .05, .025, .0125].map(h => sampleRadialCell({ ...label, direction, h }));
    const orders = Object.fromEntries(fields.map(field => [field, samples.slice(1).map((sample, i) =>
      Math.log2(samples[i][field] / sample[field]))]));
    for (const field of fields) {
      assert.ok(samples.every(s => Number.isFinite(s[field]) && s[field] > 1e-10), `${field}: nonzero finite error required`);
      for (let i = 1; i < samples.length; i++) assert.ok(samples[i][field] < .35 * samples[i - 1][field],
        `${JSON.stringify(label)} ${field}: ${JSON.stringify(samples)}`);
      for (const order of orders[field].slice(-2)) assert.ok(order > 1.8 && order < 2.2,
        `${JSON.stringify(label)} ${field}: ${JSON.stringify(orders[field])}`);
    }
    t.diagnostic(JSON.stringify({ ...label, samples, orders }));
  }
});
