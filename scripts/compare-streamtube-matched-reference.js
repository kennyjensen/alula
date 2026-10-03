// SPDX-License-Identifier: GPL-2.0-or-later
// Independent native profiles against a saved, re-evaluated coupled root.
// No nonlinear iteration or matching/tuning of acceptance thresholds.
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync } from 'node:fs';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { createContourCurve } from '../src/geometry/contour-curve.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';
const source = process.argv[2] ?? 'docs/streamtube-matched-coupled-solve.json', saved = JSON.parse(readFileSync(source)), r = saved.result;
const input = saved.restart?.input ?? saved.input, options = saved.restart?.options ?? saved.options;
const nativePath = process.argv.find(s => s.startsWith('--native='))?.slice('--native='.length) ?? 'tests/fixtures/fortran/compressible-flow.json';
const native = JSON.parse(readFileSync(nativePath));
const referenceName = process.argv.find(s => s.startsWith('--native-case='))?.slice('--native-case='.length);
const reference = native.cases.find(c => c.options.mach === input.mach && (!referenceName || c.name === referenceName));
const report = { source, sourceHash: sha256(source), nativePath, nativeHash: sha256(nativePath), startedAt: new Date().toISOString(),
  referenceName, sourceHashes: numericalSourceHashes(['scripts/compare-streamtube-matched-reference.js']),
  scope: 'Single-grid attached-flow comparison with native XFOIL. Every Euler surface pressure node, including stagnation, and every resolved BL station. Euler versus Karman–Tsien outer flow; not roundoff parity, physical force or refinement acceptance.',
  gates: { cpMax: .04, ueMax: .02, thetaRelativeMax: .05, deltaStarRelativeMax: .08 }, failures: [] };
const start = performance.now();
const interpolate = (rows, s, key) => {
  if (s < rows[0].parameter || s > rows.at(-1).parameter) throw new Error('Comparison would extrapolate beyond the native profile.');
  let i = 1; while (i < rows.length - 1 && rows[i].parameter < s) i++;
  const a = rows[i - 1], b = rows[i], t = (s - a.parameter) / (b.parameter - a.parameter);
  return (1 - t) * a[key] + t * b[key];
};
try {
  for (const [path, hash] of Object.entries(native.provenance.sha256)) assert.equal(sha256(path), hash);
  assert.ok(reference, 'No native reference at the requested Mach/case.');
  assert.equal(input.alpha, reference.options.alpha, 'Native incidence must match.');
  assert.equal(reference.expected.converged, true); assert.equal(r.converged, true);
  assert.equal(input.bodies.length, 1, 'This native fixture supplies one airfoil only.');
  assert.equal((saved.failures ?? []).length, 0); assert.equal(saved.error, undefined);
  if (saved.restart) assert.deepEqual(changedSources(saved.sourceHashes), []);
  const coupled = createCoupledStreamtubeBody(input, { ...options,
    initialBL: saved.restart ? Float64Array.from(saved.restart.initialBL) : r.x.slice(saved.initial.ne),
    initialEuler: saved.restart?.initialEuler ?? r.flow });
  const checked = coupled.evaluate(coupled.initial);
  report.checkedResidual = Math.max(...Array.from(checked.residual, Math.abs));
  assert.ok(report.checkedResidual < 1e-10);
  const lengthScale = coupled.euler.conditions.lengthScale, curve = createContourCurve(input.bodies[0].points);
  report.lengthScale = lengthScale;
  // Native XFOIL's Re is based on one input-coordinate length. The coupled
  // kernel uses distances divided by its solver length; compare viscosities,
  // not the numerically equal (but differently normalized) input values.
  report.reynoldsNormalization={kernelReynolds:options.reynolds,solverLength:lengthScale,
    effectiveReynoldsPerCoordinateUnit:options.reynolds/lengthScale,nativeReynoldsPerCoordinateUnit:reference.options.reynolds};
  assert.ok(Math.abs(options.reynolds/lengthScale/reference.options.reynolds-1)<1e-12,
    'Native/coupled viscosity mismatch: kernel Reynolds must equal native Reynolds per coordinate unit times solver length.');
  const nativeKernel = createIntegralKernel({ ...reference.options, velocityConvention: 'xfoil' });
  const parameters = reference.geometry?.parameters ?? curve.knots;
  assert.equal(parameters.length, reference.panels + 1);
  assert.equal(parameters[0], 0); assert.equal(parameters.at(-1), curve.length);
  assert.ok(parameters.every((s, i) => Number.isFinite(s) && (!i || s > parameters[i - 1])));
  const cp = reference.expected.cp.filter(q => q.index <= reference.panels + 1).map(q => {
    const parameter = parameters[q.index - 1], p = curve.evaluate(parameter).point;
    assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 1e-11 * lengthScale, 'Native samples must lie on the same solid contour.');
    return { ...q, parameter };
  });
  const { gamma, mach, pInf } = coupled.euler.conditions;
  // Euler supplies a pressure at stagnation, where there is no BL unknown.
  // Restricting Cp validation to BL stations silently omits that node.
  const pressureBranches = input.bodies.flatMap((_, body) => ['lower', 'upper'].map(side => ({ body, side })));
  report.pressureProfiles = pressureBranches.map(branch => {
    const fractions = coupled.euler.fractions[branch.body][branch.side];
    assert.equal(fractions.length, input.bodies[branch.body].trailingIndex - input.bodies[branch.body].leadingIndex + 1);
    return { body: branch.body, side: branch.side, points: fractions.map((_, k) => {
      const solid = curve.branch(branch.side, fractions[k], checked.outer.stagnation[branch.body]);
      const expected = interpolate(cp, solid.parameter, 'cp');
      const pressureCp = 2 * (checked.outer.bodyPressure(branch.body, input.bodies[branch.body].leadingIndex + k, branch.side) - pInf);
      return { k, parameter: solid.parameter, solid: solid.point, cp: pressureCp, expected, error: pressureCp - expected };
    }) };
  });
  const stagnationCp = 2 * pInf * Math.expm1(gamma / (gamma - 1) * Math.log1p(.5 * (gamma - 1) * mach ** 2));
  report.stagnationPressure = { exactCp: stagnationCp,
    points: report.pressureProfiles.map(b => ({ body: b.body, side: b.side, cp: b.points[0].cp,
      cpExcess: b.points[0].cp - stagnationCp })),
    scope: 'Independent smooth perfect-gas stagnation pressure. Positive wall-pressure excess and its refinement trend remain physical diagnostics, not a new fitted acceptance tolerance.' };
  report.uncoveredStations = [];
  report.profiles = coupled.bl.surfaces.map(branch => {
    const side = branch.side === 'upper' ? 1 : 2;
    const rows = reference.expected.bl.filter(q => q.side === side && q.index <= reference.panels + 1).map(q => ({ ...q,
      physicalUe: nativeKernel.station({ s: q.s, ue: q.ue, theta: q.theta, deltaStar: q.deltaStar, aux: q.ctau }).ue,
      parameter: parameters[q.index - 1] })).sort((a, b) => a.parameter - b.parameter);
    const fractions = coupled.euler.fractions[branch.body][branch.side];
    const points = branch.ids.map(id => {
      const raw = { ...coupled.bl.stations[id], ...checked.layers.states[id] };
      const solid = curve.branch(branch.side, fractions[raw.k], checked.outer.stagnation[branch.body]);
      const q = { ...raw, theta: raw.theta * lengthScale, deltaStar: raw.deltaStar * lengthScale };
      const covered = solid.parameter >= rows[0].parameter && solid.parameter <= rows.at(-1).parameter;
      if (!covered) report.uncoveredStations.push({ body: branch.body, side: branch.side, id, k: q.k,
        parameter: solid.parameter, nativeRange: [rows[0].parameter, rows.at(-1).parameter], solid: solid.point });
      // Retain the independently available Cp and the other covered stations
      // even if one BL point lies between native stagnation and its first
      // reported station. Missing coverage fails acceptance; never extrapolate.
      const expected = Object.fromEntries(['theta', 'deltaStar', 'physicalUe'].map(key => [key, covered ? interpolate(rows, solid.parameter, key) : null]));
      expected.cp = interpolate(cp, solid.parameter, 'cp');
      const pressure = 2 * (checked.outer.bodyPressure(branch.body, raw.i, branch.side) - pInf);
      return { ...q, parameter: solid.parameter, solid: solid.point, cp: pressure, expected,
        errors: { cp: pressure - expected.cp, ue: covered ? q.ue - expected.physicalUe : null,
          thetaRelative: covered ? q.theta / expected.theta - 1 : null, deltaStarRelative: covered ? q.deltaStar / expected.deltaStar - 1 : null } };
    });
    return { side: branch.side, points };
  });
  report.errors = {};
  report.worst = {};
  report.blStationCpMax = Math.max(...report.profiles.flatMap(b => b.points.map(p => Math.abs(p.errors.cp))));
  const pressures = report.pressureProfiles.flatMap(b => b.points.map(q => ({ ...q, side: b.side })));
  const worstPressure = pressures.reduce((a, b) => Math.abs(a.error) > Math.abs(b.error) ? a : b);
  report.errors.cpMax = Math.abs(worstPressure.error);
  report.worst.cp = { side: worstPressure.side, x: worstPressure.solid.x, k: worstPressure.k,
    cp: worstPressure.cp, expected: worstPressure.expected, error: worstPressure.error };
  for (const key of ['ue', 'thetaRelative', 'deltaStarRelative']) {
    const points = report.profiles.flatMap(p => p.points.map(q => ({ ...q, side: p.side }))).filter(q => Number.isFinite(q.errors[key]));
    if (!points.length) throw new Error(`No native comparison coverage for ${key}.`);
    const q = points.reduce((a, b) => Math.abs(a.errors[key]) > Math.abs(b.errors[key]) ? a : b);
    report.errors[key + 'Max'] = Math.abs(q.errors[key]);
    report.worst[key] = { side: q.side, x: q.solid.x, k: q.k, regime: q.regime, errors: q.errors };
  }
  if (report.uncoveredStations.length) report.failures.push(`Native BL profile does not cover ${report.uncoveredStations.length} stations; no extrapolation was performed.`);
  for (const [key, bound] of Object.entries(report.gates)) if (!(report.errors[key] <= bound)) report.failures.push(`${key}: ${report.errors[key]} exceeds ${bound}.`);
  report.status = report.failures.length ? 'native profile gates failed' : 'single-grid native profile gates passed; refinement/force acceptance remains open';
} catch (error) { report.status = 'unresolved native comparison'; report.error = error.stack ?? error.message; }
report.elapsedSeconds = (performance.now() - start) / 1000; report.changedSources = changedSources(report.sourceHashes);
if (report.changedSources.length) report.failures.push('Numerical sources changed during comparison.');
report.finishedAt = new Date().toISOString();
writeFileSync(process.argv[3] ?? 'docs/streamtube-matched-native-comparison.json', JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ status: report.status, elapsedSeconds: report.elapsedSeconds, checkedResidual: report.checkedResidual, errors: report.errors, worst: report.worst, failures: report.failures, error: report.error }));
if (report.failures.length || report.error) process.exitCode = 1;
