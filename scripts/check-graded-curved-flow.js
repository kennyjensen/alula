// SPDX-License-Identifier: GPL-2.0-or-later
// Cheap analytic local consistency checks; no global Euler/BL solve or LU.
import fs from 'node:fs';
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { gradedRinglebTube, gradedVortexTube, exactEulerPolygonFlux, exactVortexMomentumFace, simpsonIntegral } from '../tests/oracles/graded-curved-flow.js';
import { directPolygonGeometry } from '../tests/oracles/streamtube-control-volume-geometry.js';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { prepareStreamtubeTransportChain } from '../src/euler/streamtube-transport-chain.js';

const args = process.argv.slice(2);
assert(args.length === 1 && args[0].startsWith('--output='), 'Use --output=<fresh-directory>.');
const output = args[0].slice(9); assert(output && !fs.existsSync(output), 'Preserve previous receipts.');
fs.mkdirSync(output, { recursive: true });
const write = (file, data) => fs.writeFileSync(`${output}/${file}`, JSON.stringify(data, null, 2) + '\n');
const paths = ['scripts/check-graded-curved-flow.js', 'tests/oracles/graded-curved-flow.js', 'tests/oracles/ringleb.js',
  'tests/oracles/streamtube-control-volume-geometry.js', 'src/euler/streamtube-cell.js', 'src/euler/streamtube-transport-chain.js',
  'src/euler/streamtube-speed-upwind.js', 'src/geometry/simple-polygon.js', 'src/geometry/simple-quadrilateral.js',
  'node_modules/robust-orientation/orientation.js'];
const hash = p => createHash('sha256').update(fs.readFileSync(p)).digest('hex');
const hashes = Object.fromEntries(paths.map(p => [p, hash(p)])), started = performance.now();
const limits = { smoothFinalOrders: [1.8, Infinity], positiveCellArea: 0, positiveInterfacePressure: 0,
  vortexCurvatureIdentity: 2e-12, vortexFaceAbsoluteError: 2e-10, massQuadratureAbsoluteError: 1e-12 };
const report = { passed: false, physicalAcceptance: false, hashes, limits,
  operations: { globalFlowEvaluations: 0, BLClosureEvaluations: 0, jacobians: 0, linearSolves: 0, newtonUpdates: 0,
    explicitLocalCells: 0, transportChains: 0, polygonQuadratures: 0, inverseHodographMaps: 0 },
  mapping: { maximumResidual: 0, maximumCondition: 0, maximumIterations: 0 }, ringleb: [], vortex: [], failures: [],
  scope: 'Exact-field local consistency under smooth grading, thin tubes and skewness; fixed-ratio spacing is an ungated contrast. No airfoil, BL, shock, or global-root accuracy acceptance.' };
write('report.json', report);
const mid = (a, b) => ({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 });
const cross = (a, b) => a.x * b.y - a.y * b.x;
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const length = p => Math.hypot(p.x, p.y);
const maximum = values => Math.max(...values.map(Math.abs));
const registerMapping = s => {
  if (!s.mapping) return;
  report.operations.inverseHodographMaps++;
  report.mapping.maximumResidual = Math.max(report.mapping.maximumResidual, s.mapping.residual);
  report.mapping.maximumCondition = Math.max(report.mapping.maximumCondition, s.mapping.maxCondition);
  report.mapping.maximumIterations = Math.max(report.mapping.maximumIterations, s.mapping.iterations);
};
const dual = (lower, upper) => [mid(lower[0], lower[1]), lower[1], mid(lower[1], lower[2]),
  mid(upper[1], upper[2]), upper[1], mid(upper[0], upper[1])];
const ordered = levels => Object.fromEntries(Object.keys(levels[0].errors).map(key => [key,
  levels.slice(1).map((b, i) => levels[i].errors[key] > 1e-13 && b.errors[key] > 1e-13
    ? Math.log2(levels[i].errors[key] / b.errors[key]) : null)]));

function measure(tube, index, { filtered = false, polygon = false, vortex = false } = {}) {
  tube.exactSections.forEach(registerMapping);
  const lower = tube.lower.slice(index - 1, index + 2), upper = tube.upper.slice(index - 1, index + 2);
  const local = { lower, upper, densities: tube.densities.slice(index - 1, index + 1),
    massFlow: tube.massFlow, stagnationEnthalpy: tube.stagnationEnthalpy, gamma: tube.gamma };
  const plain = evaluateStreamtubeCell(local); report.operations.explicitLocalCells++;
  let biased;
  if (filtered) {
    const chain = prepareStreamtubeTransportChain({ ...tube, upwind: { mucon: 1, mcrit: .99, boundary: { kind: 'unfiltered-first-two' } } });
    report.operations.transportChains++;
    biased = evaluateStreamtubeCell({ ...local, transportSpeeds: chain.transportSpeeds.slice(index - 1, index + 1) });
    report.operations.explicitLocalCells++;
  }
  const ref = tube.exactSections.slice(index - 1, index + 1), pRef = (ref[0].p + ref[1].p) / 2;
  const ds = length(plain.geometry.streamwise), scale = pRef * ds, shape = directPolygonGeometry(dual(lower, upper));
  const quads = [0, 1].map(i => directPolygonGeometry([lower[i], lower[i + 1], upper[i + 1], upper[i]]));
  assert(shape.valid && quads.every(q => q.valid && q.minCornerSine > 0));
  const errors = { R1Density: Math.abs(plain.streamwiseResidual / scale), R2Density: Math.abs(plain.isentropicResidual / scale),
    R1MinusR2Density: Math.abs((plain.streamwiseResidual - plain.isentropicResidual) / scale),
    speedRelative: maximum(plain.states.map((s, i) => s.q / ref[i].q - 1)),
    pressureRelative: maximum(plain.states.map((s, i) => s.p / ref[i].p - 1)),
    direction: maximum(ref.map((s, i) => (s.u * plain.geometry.directions[i].y - s.v * plain.geometry.directions[i].x) / s.q)) };
  if (biased) Object.assign(errors, { biasedR1Density: Math.abs(biased.streamwiseResidual / scale),
    biasedR2Density: Math.abs(biased.isentropicResidual / scale) });
  assert(Object.values(errors).every(Number.isFinite));
  const record = { controls: tube.controls, errors, signed: { R1Density: plain.streamwiseResidual / scale,
    R2Density: plain.isentropicResidual / scale, Pc: plain.pressureCorrection },
    state: { minPressure: Math.min(...plain.states.map(s => s.p)), interfacePressure: plain.interfacePressure,
      exactMach: ref.map(s => s.mach), reconstructedMach: plain.states.map(s => Math.sqrt(s.machSquared)) },
    geometry: { minimumCornerSine: Math.min(...quads.map(q => q.minCornerSine)), minimumArea: Math.min(...quads.map(q => q.area)),
      streamwiseLength: ds, normalAreas: plain.geometry.normalAreas, curvature: plain.geometry.pressureCurvature,
      normalizedWidth: Math.min(...plain.geometry.normalAreas) / ds,
      rawNodeAdjacentLengthRatio: length(sub(lower[2], lower[1])) / length(sub(lower[1], lower[0])) },
    sampling: { points: ref.map(s => ({ x: s.x, y: s.y })), exactDensities: ref.map(s => s.rho),
      exactMass: tube.massFlow, densitiesUnmodified: true } };
  if (vortex) {
    const [a, b] = [tube.angles[1] - tube.angles[0], tube.angles[2] - tube.angles[1]];
    const expected = 4 * Math.tan(a / 2) * Math.tan(b / 2);
    record.vortex = { analyticCurvature: expected, curvatureError: plain.geometry.pressureCurvature - expected,
      pcDirectStreamwise: -plain.pressureCorrection * cross(plain.geometry.sides.lower, plain.geometry.sides.upper) / plain.geometry.area,
      massQuadratureDifference: tube.massQuadratureDifference,
      exactChordMass: [0, 1].map(i => {
        const inner = length(mid(upper[i], upper[i + 1])), outer = length(mid(lower[i], lower[i + 1]));
        return simpsonIntegral(r => { const s = tube.exactAt({ x: r, y: 0 }); return s.rho * s.q; }, inner, outer, 128);
      }) };
    record.vortex.chordMassRelativeMismatch = record.vortex.exactChordMass.map(m => m / tube.massFlow - 1);
    assert(Math.abs(record.vortex.curvatureError) < limits.vortexCurvatureIdentity);
    assert(tube.massQuadratureDifference < limits.massQuadratureAbsoluteError);
  }
  if (polygon) {
    const points = dual(lower, upper), exactAt = p => { const s = tube.exactAt(p); registerMapping(s); return s; };
    const coarse = exactEulerPolygonFlux(points, exactAt, 32), fine = exactEulerPolygonFlux(points, exactAt, 64);
    report.operations.polygonQuadratures += 2;
    const difference = maximum(fine.faces.flatMap((f, i) => f.map((v, j) => v - coarse.faces[i][j])));
    record.exactPolygon = { points, faces: fine.faces, total: fine.total, maximumQuadratureDifference: difference,
      normalizedMomentumClosure: Math.hypot(fine.total[1], fine.total[2]) / (pRef * shape.area),
      absoluteLateralMassFlux: maximum([fine.faces[0][0], fine.faces[1][0], fine.faces[3][0], fine.faces[4][0]]) };
    if (vortex) {
      record.exactPolygon.endpointMomentum = points.map((a, i) => exactVortexMomentumFace(a, points[(i + 1) % points.length]));
      const error = maximum(record.exactPolygon.endpointMomentum.flatMap((v, i) => [v.x - fine.faces[i][1], v.y - fine.faces[i][2]]));
      record.exactPolygon.maximumFaceMomentumOracleError = error;
      assert(error < limits.vortexFaceAbsoluteError);
    }
  }
  return record;
}

try {
  // These finite controls and expected smooth-family order are set before any
  // run. Fixed-ratio spacing deliberately does not become smoother on refinement.
  for (const mach of [.9, 1.04]) for (const shape of ['regular', 'thin', 'thin-skew']) for (const direction of [1, -1]) {
    const levels = [0, 1, 2, 3].map(level => measure(gradedRinglebTube({ mach, h: .002 / 2 ** level,
      widthRatio: shape === 'regular' ? 2 : .02, skew: shape === 'thin-skew', direction }), 3,
    { filtered: true, polygon: level === 0 }));
    const item = { mach, shape, direction, levels, orders: ordered(levels) }; report.ringleb.push(item);
    for (const key of ['R1Density', 'R2Density', 'biasedR1Density', 'biasedR2Density', 'speedRelative']) {
      const order = item.orders[key].at(-1);
      if (order === null || order < limits.smoothFinalOrders[0]) report.failures.push({ kind: 'smooth-order', mach, shape, direction, key, order });
    }
  }
  for (const direction of [1, -1]) {
    const levels = [0, 1, 2, 3].map(level => measure(gradedRinglebTube({ mach: .9, h: .002 / 2 ** level,
      widthRatio: .02, skew: true, direction, spacing: 'fixed-ratio' }), 3));
    report.ringleb.push({ mach: .9, shape: 'thin-skew-fixed-ratio', direction, gated: false, levels, orders: ordered(levels) });
  }
  for (const widthRatio of [1, .02]) {
    const levels = [0, 1, 2, 3].map(level => measure(gradedVortexTube({ h: .1 / 2 ** level, widthRatio }), 1,
      { polygon: true, vortex: true }));
    const item = { widthRatio, levels, orders: ordered(levels) }; report.vortex.push(item);
    for (const key of ['R1Density', 'R2Density', 'speedRelative']) {
      const order = item.orders[key].at(-1);
      if (order === null || order < limits.smoothFinalOrders[0]) report.failures.push({ kind: 'smooth-vortex-order', widthRatio, key, order });
    }
  }
  report.passed = report.failures.length === 0;
} catch (error) { report.failures.push({ message: error.message, stack: error.stack }); }
report.seconds = (performance.now() - started) / 1000;
report.sourceChanges = paths.filter(p => hash(p) !== hashes[p]);
if (report.sourceChanges.length) report.passed = false;
write('report.json', report);
write('receipt.json', { passed: report.passed, seconds: report.seconds, hashes, sourceChanges: report.sourceChanges,
  failures: report.failures, operations: report.operations, physicalAcceptance: false, reportSha256: hash(`${output}/report.json`) });
console.log(JSON.stringify({ passed: report.passed, seconds: report.seconds, operations: report.operations, mapping: report.mapping,
  ringleb: report.ringleb.map(c => ({ mach: c.mach, shape: c.shape, direction: c.direction, orders: c.orders,
    finestErrors: c.levels.at(-1).errors, finestGeometry: c.levels.at(-1).geometry })),
  vortex: report.vortex.map(c => ({ widthRatio: c.widthRatio, orders: c.orders })), failures: report.failures }));
if (!report.passed) process.exitCode = 1;
