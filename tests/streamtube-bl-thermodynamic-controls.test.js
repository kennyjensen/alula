import test from 'node:test';
import assert from 'node:assert/strict';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { streamtubeEdgeVelocity } from '../src/euler/streamtube-edge-velocity.js';
import { streamtubeBLEdgeThermodynamics } from '../src/euler/streamtube-bl-edge-thermodynamics.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { referenceGas, isentropeAtSpeed, isentropeAtMach, normalShock, pressureTaylorBound } from './oracles/thermodynamic-collocation.js';

const gas = referenceGas(), conditions = { ...gas, flowModel: 'compressible' };
const tolerance = 512 * Number.EPSILON;
const relativeError = (a, b) => Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b));
const near = (a, b, message) => assert.ok(relativeError(a, b) <= tolerance, `${message}: ${a} != ${b}`);
const station = ue => ({ s: 1, ue, theta: .001, deltaStar: .003, aux: .03 });

function gasAtMatchedStation(cell) {
  const edge = streamtubeEdgeVelocity(cell), physical = station(edge.ue);
  const diagnostic = streamtubeBLEdgeThermodynamics({
    flow: { cells: [[[cell], [cell]]] },
    euler: { conditions, layout: { nx: 2, tubes: [1, 1] } },
    bl: { stations: [{ id: 0, kind: 'surface', body: 0, side: 'upper', i: 1 }] },
    states: [physical],
  }).stations[0];
  const native = createIntegralKernel({ mach: gas.freestreamMach, gamma: gas.gamma, reynolds: 1e6 }).station(physical, 'turbulent');
  const pressure = diagnostic.banks[0].interfacePressure;
  const inferredDensity = gas.gamma * pressure / ((gas.gamma - 1) * (gas.h0 - .5 * edge.ue ** 2));
  const historicalPressure = (gas.gamma - 1) / gas.gamma * native.rho * diagnostic.boundaryLayer.enthalpy;
  near(native.rho, diagnostic.boundaryLayer.rho, 'Native and diagnostic historical density');
  near(native.rho / inferredDensity - 1, historicalPressure / pressure - 1, 'Collocated pressure/density identity');
  return { edge, native, diagnostic, inferredDensity, historicalPressure,
    excess: native.rho / inferredDensity - 1 };
}

test('uniform postshock cell isolates historical BL density excess from all interpolation and speed-bias effects', () => {
  const measurements = [];
  for (const upstreamMach of [1, 1.2, 2]) {
    const shock = normalShock(upstreamMach, gas), a = shock.upstream, b = shock.downstream;
    near(a.rho * a.q, b.rho * b.q, 'Normal-shock mass');
    near(a.p + a.rho * a.q ** 2, b.p + b.rho * b.q ** 2, 'Normal-shock momentum');
    near(b.enthalpy + .5 * b.q ** 2, gas.h0, 'Normal-shock total enthalpy');
    near(b.rho * b.q ** 2 / (gas.gamma * b.p), b.machSquared, 'Independent downstream Mach relation');
    const cellInput = { lower: [{ x: -1, y: 0 }, { x: 0, y: 0 }, { x: 1, y: 0 }],
      upper: [{ x: -1, y: 1 }, { x: 0, y: 1 }, { x: 1, y: 1 }],
      densities: [b.rho, b.rho], massFlow: b.rho * b.q, stagnationEnthalpy: gas.h0,
      gamma: gas.gamma, transportSpeeds: [b.q, b.q] };
    const before = structuredClone(cellInput), cell = evaluateStreamtubeCell(cellInput), matched = gasAtMatchedStation(cell);
    assert.deepEqual(cellInput, before);
    assert.equal(cell.geometry.pressureCurvature, 0); assert.ok(matched.edge.correction === 0);
    near(cell.streamwiseResidual, 0, 'Uniform streamwise momentum');
    near(cell.isentropicResidual, 0, 'Uniform entropy row');
    near(cell.entropyJump, 0, 'Uniform downstream entropy difference');
    near(cell.artificialEntropyJump, 0, 'Uniform artificial entropy difference');
    near(matched.edge.ue, b.q, 'Collocated physical speed');
    near(matched.diagnostic.banks[0].interfacePressure, b.p, 'Collocated physical pressure');
    near(matched.inferredDensity, b.rho, 'Local pressure/speed downstream density');
    const expectedExcess = Math.expm1(shock.entropyOverR);
    near(matched.excess, expectedExcess, 'Historical postshock excess');
    for (const section of matched.diagnostic.banks[0].sections)
      near(section.entropyOverR, shock.entropyOverR, 'Independent section entropy');
    if (upstreamMach === 1) near(matched.excess, 0, 'Zero-entropy control');
    else assert.ok(shock.entropyOverR > 0 && matched.excess > 0);
    measurements.push({ upstreamMach, downstreamMach: Math.sqrt(b.machSquared),
      entropyOverR: shock.entropyOverR, expectedExcess, measuredExcess: matched.excess,
      excessAbsoluteError: Math.abs(matched.excess - expectedExcess),
      physicalDensity: b.rho, collocatedDensity: matched.inferredDensity, historicalDensity: matched.native.rho,
      relativeCollocatedDensityError: matched.inferredDensity / b.rho - 1 });
  }
  console.log('THERMODYNAMIC_CONTROL ' + JSON.stringify({ kind: 'uniform-postshock', measurements,
    tolerance, localCellEvaluations: 3, localNativeStationEvaluations: 3, localEdgeDiagnosticEvaluations: 3 }));
});

test('smooth isentropic collocation has the independently bounded quadratic pressure defect and factor-four refinement', () => {
  const measurements = [];
  for (const centerMach of [.65, 1.3]) {
    const center = isentropeAtMach(centerMach, gas), massFlow = center.rho * center.q;
    const levels = [];
    for (const fraction of [.08, .04, .02, .01]) {
      const deltaQ = fraction * center.q;
      const a = isentropeAtSpeed(center.q - .5 * deltaQ, gas), b = isentropeAtSpeed(center.q + .5 * deltaQ, gas);
      const areaA = massFlow / (a.rho * a.q), areaB = massFlow / (b.rho * b.q);
      const widthCenter = .5 * (areaA + areaB), widths = [2 * areaA - widthCenter, widthCenter, 2 * areaB - widthCenter];
      assert.ok(widths.every(w => w > 0));
      // Symmetric straight, linearly tapered banks give zero normal inertia
      // and zero curvature correction. Section normal areas match mass flow.
      const cell = evaluateStreamtubeCell({ lower: widths.map((w, i) => ({ x: i - 1, y: -.5 * w })),
        upper: widths.map((w, i) => ({ x: i - 1, y: .5 * w })),
        densities: [a.rho, b.rho], massFlow, stagnationEnthalpy: gas.h0, gamma: gas.gamma,
        transportSpeeds: [a.q, b.q] });
      const matched = gasAtMatchedStation(cell), meanPressure = .5 * (a.p + b.p);
      near(cell.geometry.pressureCurvature, 0, 'Straight taper curvature');
      near(matched.edge.correction, 0, 'No sawtooth correction');
      near(matched.edge.ue, center.q, 'Symmetric speed average');
      near(matched.diagnostic.banks[0].interfacePressure, meanPressure, 'Actual mean interface pressure');
      near(matched.historicalPressure, center.p, 'Common-isentrope pressure at mean speed');
      for (const s of matched.diagnostic.banks[0].sections) near(s.entropyOverR, 0, 'Zero endpoint entropy');
      near(cell.entropyJump, 0, 'No physical entropy jump');
      const difference = matched.diagnostic.banks[0].interfacePressure - matched.historicalPressure;
      const analytic = pressureTaylorBound(center.q, deltaQ, gas);
      const roundoffBound = tolerance * Math.max(1, center.p, meanPressure);
      const remainder = Math.abs(difference - analytic.leadingDifference);
      assert.ok(remainder <= analytic.remainderBound + roundoffBound, `Taylor remainder ${remainder} exceeds its independent bound`);
      assert.equal(Math.sign(difference), Math.sign(analytic.pSecond));
      near(matched.excess, center.p / meanPressure - 1, 'Pure sampling density discrepancy');
      levels.push({ deltaQ, difference, collocatedHistoricalDensityExcess: matched.excess,
        pSecond: analytic.pSecond, leadingDifference: analytic.leadingDifference,
        remainder, fourthOrderRemainderBound: analytic.remainderBound, roundoffBound });
    }
    const ratios = levels.slice(1).map((fine, i) => {
      const coarse = levels[i], ratio = coarse.difference / fine.difference;
      const permittedDeparture = (coarse.fourthOrderRemainderBound + coarse.roundoffBound
        + 4 * (fine.fourthOrderRemainderBound + fine.roundoffBound)) / Math.abs(fine.difference);
      assert.ok(permittedDeparture < .04, 'Chosen interval must resolve a factor-four limit to better than one percent');
      assert.ok(Math.abs(ratio - 4) <= permittedDeparture + tolerance, 'Refinement ratio exceeds Taylor/roundoff bound');
      return { ratio, permittedDepartureFromFour: permittedDeparture, order: Math.log2(Math.abs(ratio)) };
    });
    measurements.push({ centerMach, centerSpeed: center.q, levels, ratios });
  }
  console.log('THERMODYNAMIC_CONTROL ' + JSON.stringify({ kind: 'smooth-isentropic-collocation', measurements,
    tolerance, localCellEvaluations: 8, localNativeStationEvaluations: 8, localEdgeDiagnosticEvaluations: 8,
    scope: 'Consistency/sampling control, not a globally solved nozzle or a curved-interface pressure validation.' }));
});
