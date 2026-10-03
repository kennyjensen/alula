// SPDX-License-Identifier: GPL-2.0-or-later
// Recheck immutable strict evidence under the separate temporary policy.
import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { createIntegralKernel } from '../src/viscous/integral.js';
import { createPhysicalBLDomain } from '../src/viscous/physical-domain.js';
import { streamtubeMeshSnapshot } from '../src/euler/streamtube-mesh-preview.js';
import { directStreamtubeVolumeGeometry } from '../tests/oracles/streamtube-control-volume-geometry.js';
import { evaluateCoupledAcceptance, coupledAcceptancePolicy, originalGates } from './validation/coupled-acceptance.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';

const producer = 'scripts/check-current-coupled-acceptance.js';
const args = process.argv.slice(2), option = (key, fallback) => args.find(s => s.startsWith(`--${key}=`))?.slice(key.length + 3) ?? fallback;
assert.ok(args.every(s => /^--(comparison|output)=.+/.test(s)), 'Only --comparison and --output are supported.');
const comparisonPath = option('comparison', 'docs/current-natural-outer-refined-stagnation-native.json');
const output = option('output', 'docs/current-coupled-acceptance.json');
assert.ok(output.endsWith('.json'), 'Output must end in .json.');
const logPath = output.replace(/\.json$/, '.log');
// These are explicit source transitions, not a general stale-source exception.
// Old bytes must match the preserved archive and new bytes must match this
// reviewed version. Every saved state/profile value is still replayed below.
const sourceArchives = [
  { manifest: 'docs/stagnation-density-startup/before.json', files: [
    { path: 'src/euler/streamtube-result.js', archive: 'docs/stagnation-density-startup/streamtube-result.js.txt',
      sha256: 'ef229ee18ed28a365ae18a040ece274440f9e021101bb6aaac223964843cf337',
      currentSha256: 'f88439f12d15e4c26a5b58343348f23ef4d3e84a34580aa2b5f647a4d440a072' },
    { path: 'src/euler/streamtube-coupled-assembly.js', archive: 'docs/stagnation-density-startup/streamtube-coupled-assembly.js.txt',
      sha256: 'f49a202baae72a0ba1fa54fcb1b1f50fbc3f90a6033a20f90d2b21809dbbbe71',
      currentSha256: '58b56e8198166b991b83b3b88da81e5a9a2df5077f97e25bae63521217e7a9a4' },
  ] },
  { manifest: 'docs/transition-coordinate-scale/before.json', files: [
    { path: 'src/viscous/transition-interval.js', archive: 'docs/transition-coordinate-scale/transition-interval-before.js.txt',
      sha256: 'f82cac66b5b70a8406e25b9c3da8a9f37b292bae70c2046b62cb9ce61a816b96',
      currentSha256: '29433007a830f05399228f2940a8132c7001c458ed02ae7b7a224d60e918c638' },
  ] },
];
const producerArchive = { manifest: 'docs/transition-coordinate-scale/acceptance-producer-before.json',
  path: producer, archive: 'docs/transition-coordinate-scale/check-current-coupled-acceptance-before.js.txt',
  sha256: 'c882c4545a6d3a1337f9b6485259b265f169746359e209839264edb6d3b5fc6a' };
// Repeating this new policy check is supported. Its input and original strict
// producer/report/log can never be overwritten by the selected output paths.
const protectedPaths = new Set([comparisonPath, comparisonPath.replace(/\.json$/, '.log'),
  producer, 'scripts/validation/coupled-acceptance.js', 'scripts/compare-streamtube-matched-reference.js',
  producerArchive.manifest, producerArchive.archive,
  ...sourceArchives.flatMap(a => [a.manifest, ...a.files.flatMap(f => [f.path, f.archive])])].map(p => path.resolve(p)));
const savedComparison = JSON.parse(fs.readFileSync(comparisonPath));
for (const p of [savedComparison.source, savedComparison.nativePath]) if (typeof p === 'string') protectedPaths.add(path.resolve(p));
for (const p of [output, logPath]) assert.ok(!protectedPaths.has(path.resolve(p)), `Output overlaps protected input ${p}.`);
fs.writeFileSync(logPath, '');
const start = performance.now(), inputs = {}, sourceHashes = numericalSourceHashes([producer,
  'scripts/validation/coupled-acceptance.js', 'scripts/compare-streamtube-matched-reference.js',
  'tests/oracles/streamtube-control-volume-geometry.js', 'tests/coupled-acceptance.test.js', 'package.json']);
const report = { date: new Date().toISOString(), inProgress: true, passed: false, accepted: false,
  physicalAcceptance: false, fullSolverComplete: false, runtimeChanged: false,
  scope: coupledAcceptancePolicy.scope, policy: coupledAcceptancePolicy,
  originalGates: { ...coupledAcceptancePolicy.originalGates }, currentEffective: { ...coupledAcceptancePolicy.currentEffective },
  sourceHashes, inputs, verifiedSourceArchives: [], historicalSourceChanges: [], operations: { coupledEvaluations: 0, geometryAudits: 0, nativeRuns: 0,
    jacobianAssemblies: 0, factorizations: 0, newtonUpdates: 0 }, failures: [] };
const save = () => fs.writeFileSync(output, JSON.stringify(report, null, 2) + '\n');
const log = value => { const line = JSON.stringify({ time: new Date().toISOString(), ...value });
  fs.appendFileSync(logPath, line + '\n'); console.log(line); };
const track = (p, hash) => { const current = sha256(p); if (hash !== undefined) assert.equal(current, hash, `Stale linked input: ${p}`);
  inputs[p] = current; return current; };
const verifiedTransitions = new Map();
const verifyArchives = () => {
  for (const a of sourceArchives) {
    const manifestHash = track(a.manifest), manifest = JSON.parse(fs.readFileSync(a.manifest));
    assert.deepEqual(manifest.files, a.files.map(({ currentSha256, ...old }) => old), `Unexpected archive entries: ${a.manifest}`);
    for (const f of a.files) {
      track(f.archive, f.sha256); track(f.path, f.currentSha256);
      assert.ok(!verifiedTransitions.has(f.path), `Duplicate archived transition: ${f.path}`);
      verifiedTransitions.set(f.path, { ...f, manifest: a.manifest, manifestSha256: manifestHash });
    }
    if (manifest.baseline) track(manifest.baseline.path, manifest.baseline.sha256);
    report.verifiedSourceArchives.push({ path: a.manifest, sha256: manifestHash, files: a.files });
  }
  const { manifest: p, ...old } = producerArchive;
  const manifestHash = track(p); assert.deepEqual(JSON.parse(fs.readFileSync(p)).files, [old]);
  track(old.archive, old.sha256);
  report.producerSourceChange = { manifest: { path: p, sha256: manifestHash }, ...old,
    currentSha256: sourceHashes[producer] };
};
const verifyHashes = (hashes, evidence) => {
  assert.ok(hashes && Object.keys(hashes).length > 0, 'Missing source provenance.');
  for (const [p, h] of Object.entries(hashes)) {
    const current = sha256(p);
    if (current === h) { track(p, h); continue; }
    const archived = verifiedTransitions.get(p);
    assert.ok(archived, `Unarchived source change: ${p}`);
    assert.equal(h, archived.sha256, `Historical source does not match its archive: ${p}`);
    track(p, archived.currentSha256);
    let change = report.historicalSourceChanges.find(c => c.path === p);
    if (!change) {
      change = { path: p, historicalSha256: h, currentSha256: current,
        archive: { path: archived.archive, sha256: archived.sha256 },
        manifest: { path: archived.manifest, sha256: archived.manifestSha256 }, evidence: [] };
      report.historicalSourceChanges.push(change);
    }
    if (!change.evidence.includes(evidence)) change.evidence.push(evidence);
  }
};
save(); log({ event: 'start', comparisonPath, policy: report.policy });
try {
  verifyArchives();
  const c = savedComparison; track(comparisonPath); track(comparisonPath.replace(/\.json$/, '.log'));
  assert.equal(c.error, undefined); assert.notEqual(c.inProgress, true);
  assert.ok(typeof c.finishedAt === 'string'); assert.deepEqual(c.gates, originalGates);
  assert.deepEqual(c.changedSources, []); verifyHashes(c.sourceHashes, comparisonPath);
  track(c.source, c.sourceHash); track(c.nativePath, c.nativeHash);
  const saved = JSON.parse(fs.readFileSync(c.source)), native = JSON.parse(fs.readFileSync(c.nativePath));
  assert.equal(saved.inProgress, false); assert.equal(saved.error, undefined); assert.deepEqual(saved.failures ?? [], []);
  assert.equal(native.inProgress, false); assert.equal(native.passed, true); assert.equal(native.error, undefined);
  assert.deepEqual(native.failures, []); verifyHashes(saved.sourceHashes, c.source); verifyHashes(native.sourceHashes, c.nativePath);
  // Native executable/data provenance must remain exact, without source-transition exceptions.
  assert.ok(native.provenance?.sha256 && Object.keys(native.provenance.sha256).length > 0, 'Missing native provenance.');
  for (const [p, h] of Object.entries(native.provenance.sha256)) track(p, h);
  report.comparison = { path: comparisonPath, sha256: inputs[comparisonPath], originalStatus: c.status,
    originalFailures: structuredClone(c.failures), source: { path: c.source, sha256: c.sourceHash },
    native: { path: c.nativePath, sha256: c.nativeHash, case: c.referenceName }, originalEvidencePreserved: true };
  const f = saved.checkpoint?.restart;
  assert.ok(f?.initialEuler?.x && f.initialEuler.nodes && f.initialBL, 'Complete supplied Euler/grid/BL restart required; no cold solve.');
  assert.deepEqual(f, saved.restart); assert.deepEqual(saved.checkpoint.families, saved.result.families);
  assert.equal(saved.result.converged, true, 'The temporary policy requires a converged coupled result.');
  if (f.options.transitionMode === 'automatic') assert.ok(Array.isArray(f.options.transitionState), 'Explicit phase map required.');
  const immutableCheckpoint = JSON.stringify(saved.checkpoint);
  const system = createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  const x = system.initial, phase = system.bl.snapshotActive(), value = system.evaluate(x); report.operations.coupledEvaluations++;
  assert.deepEqual(value.families, saved.result.families); assert.deepEqual(value.outer.nodes, f.initialEuler.nodes);
  assert.deepEqual(Array.from(x.slice(system.ne)), f.initialBL); assert.deepEqual(phase, f.options.transitionState);
  const residual = Math.max(...Array.from(value.residual, Math.abs)); assert.equal(residual, c.checkedResidual);
  assert.ok(value.residual.every(Number.isFinite));
  const domain = createPhysicalBLDomain(system.bl.kernel.parameters);
  let minimumShape = Infinity, minimumEnthalpy = Infinity;
  for (const s of value.layers.states) { const d = domain(s); minimumShape = Math.min(minimumShape, d.shape);
    minimumEnthalpy = Math.min(minimumEnthalpy, d.enthalpy); assert.ok(d.shape > 0 && d.enthalpy > 0); }
  const mesh = streamtubeMeshSnapshot({ system: system.euler, nodes: value.outer.nodes }); assert.equal(mesh.quality.valid, true);
  const geometry = directStreamtubeVolumeGeometry(value.outer.nodes); report.operations.geometryAudits++;
  assert.equal(geometry.valid, true); assert.deepEqual(geometry.concavePrimal, []);
  report.admissibility = { accepted: true, finiteStateAndResidual: true, minimumShape, minimumEnthalpy,
    productionMeshQuality: mesh.quality, independentGeometry: geometry };
  report.exactReplay = { families: true, physicalNodes: true, packedBL: true, transitionState: true,
    checkedResidual: true, fullSavedResidualVectorAvailable: false };
  report.state = { path: c.source, sha256: c.sourceHash, unknowns: system.n,
    converged: saved.result.converged, reason: saved.result.reason, families: value.families, maximumResidual: residual };

  // Reconstruct every strict comparison value from the actual state and native
  // rows. The error maxima and coverage are therefore not accepted on trust.
  const reference = native.cases.find(q => q.name === c.referenceName);
  assert.ok(reference); assert.equal(reference.panels, 473); assert.equal(reference.expected.converged, true);
  assert.equal(reference.options.mach, f.input.mach); assert.equal(reference.options.alpha, f.input.alpha);
  assert.equal(reference.options.ncrit, f.options.ncrit); assert.deepEqual(reference.options.trips, f.options.tripFractions[0]);
  assert.equal(f.input.bodies.length, 1, 'Current native evidence is a single-element comparison.');
  const L = system.euler.conditions.lengthScale, curve = system.euler.curves[0], parameters = reference.geometry.parameters;
  assert.equal(L, c.lengthScale); assert.ok(Math.abs(f.options.reynolds / L / reference.options.reynolds - 1) < 1e-12);
  assert.equal(parameters.length, reference.panels + 1); assert.equal(parameters[0], 0); assert.equal(parameters.at(-1), curve.length);
  assert.ok(parameters.every((p, k) => Number.isFinite(p) && (!k || p > parameters[k - 1])));
  const kernel = createIntegralKernel({ ...reference.options, velocityConvention: 'xfoil' });
  const nativeCp = reference.expected.cp.filter(q => q.index <= reference.panels + 1).map(q => {
    const parameter = parameters[q.index - 1], p = curve.evaluate(parameter).point;
    assert.ok(Math.hypot(p.x - q.x, p.y - q.y) < 1e-11 * L); return { ...q, parameter };
  });
  const interpolate = (rows, parameter, key) => {
    assert.ok(parameter >= rows[0].parameter && parameter <= rows.at(-1).parameter, 'Native coverage is incomplete; no extrapolation.');
    let i = 1; while (i < rows.length - 1 && rows[i].parameter < parameter) i++;
    const a = rows[i - 1], b = rows[i], t = (parameter - a.parameter) / (b.parameter - a.parameter);
    return (1 - t) * a[key] + t * b[key];
  };
  assert.equal(c.pressureProfiles.length, 2); assert.equal(c.profiles.length, 2); assert.deepEqual(c.uncoveredStations, []);
  const errors = { cpMax: 0, ueMax: 0, thetaRelativeMax: 0, deltaStarRelativeMax: 0 };
  const coverage = { requestedBL: 0, coveredBL: 0, requestedPressure: 0, coveredPressure: 0, uncoveredStations: 0 };
  for (const side of ['upper', 'lower']) {
    const branch = system.bl.surfaces.find(b => b.body === 0 && b.side === side), body = system.euler.layout.bodies[0];
    const fractions = system.euler.fractions[0][side], pressures = c.pressureProfiles.find(b => b.body === 0 && b.side === side);
    const profile = c.profiles.find(b => b.side === side); assert.ok(branch && pressures && profile);
    assert.equal(pressures.points.length, fractions.length); assert.equal(profile.points.length, branch.ids.length);
    assert.equal(pressures.points.length, profile.points.length + 1);
    coverage.requestedPressure += fractions.length; coverage.requestedBL += branch.ids.length;
    for (let k = 0; k < fractions.length; k++) {
      const solid = curve.branch(side, fractions[k], value.outer.stagnation[0]), expected = interpolate(nativeCp, solid.parameter, 'cp');
      const cp = 2 * (value.outer.bodyPressure(0, body.leadingIndex + k, side) - system.euler.conditions.pInf);
      assert.deepEqual(pressures.points[k], { k, parameter: solid.parameter, solid: solid.point, cp, expected, error: cp - expected });
      errors.cpMax = Math.max(errors.cpMax, Math.abs(cp - expected)); coverage.coveredPressure++;
    }
    const rows = reference.expected.bl.filter(q => q.side === (side === 'upper' ? 1 : 2) && q.index <= reference.panels + 1)
      .map(q => ({ ...q, parameter: parameters[q.index - 1],
        physicalUe: kernel.station({ s: q.s, ue: q.ue, theta: q.theta, deltaStar: q.deltaStar, aux: q.ctau }).ue }))
      .sort((a, b) => a.parameter - b.parameter);
    assert.ok(rows.length > 1);
    for (let k = 0; k < branch.ids.length; k++) {
      const id = branch.ids[k], raw = { ...system.bl.stations[id], ...value.layers.states[id] };
      const solid = curve.branch(side, fractions[raw.k], value.outer.stagnation[0]);
      const q = { ...raw, theta: raw.theta * L, deltaStar: raw.deltaStar * L };
      const expected = Object.fromEntries(['theta', 'deltaStar', 'physicalUe'].map(key => [key, interpolate(rows, solid.parameter, key)]));
      expected.cp = interpolate(nativeCp, solid.parameter, 'cp');
      const cp = 2 * (value.outer.bodyPressure(0, raw.i, side) - system.euler.conditions.pInf);
      const pointErrors = { cp: cp - expected.cp, ue: q.ue - expected.physicalUe,
        thetaRelative: q.theta / expected.theta - 1, deltaStarRelative: q.deltaStar / expected.deltaStar - 1 };
      assert.deepEqual(profile.points[k], { ...q, parameter: solid.parameter, solid: solid.point, cp, expected, errors: pointErrors });
      for (const key of ['ue', 'thetaRelative', 'deltaStarRelative']) errors[key + 'Max'] = Math.max(errors[key + 'Max'], Math.abs(pointErrors[key]));
      coverage.coveredBL++;
    }
  }
  assert.deepEqual(errors, c.errors);
  assert.deepEqual(c.failures, Object.entries(originalGates).filter(([key, limit]) => errors[key] > limit)
    .map(([key, limit]) => `${key}: ${errors[key]} exceeds ${limit}.`));
  report.exactReplay.allComparisonPoints = true; report.exactReplay.nativeInterpolants = true; report.exactReplay.errorMaxima = true;
  assert.deepEqual(system.bl.snapshotActive(), phase); assert.equal(JSON.stringify(saved.checkpoint), immutableCheckpoint);
  const acceptance = evaluateCoupledAcceptance({ complete: true, converged: saved.result.converged, admissible: true,
    provenanceValid: true, exactReplay: true, residual, errors, coverage });
  Object.assign(report, acceptance); report.worst = structuredClone(c.worst);
  report.errorsUnaltered = true; report.passed = acceptance.accepted;
} catch (error) { report.passed = false; report.accepted = false; report.failures.push(error.message); report.error = error.stack; }
report.changedSources = changedSources(sourceHashes);
report.changedInputs = Object.entries(inputs).filter(([p, h]) => sha256(p) !== h).map(([p]) => p);
if (report.changedSources.length || report.changedInputs.length) {
  report.failures.push('Source or immutable evidence changed during the check.'); report.passed = false; report.accepted = false;
}
report.inProgress = false; report.seconds = (performance.now() - start) / 1000; save();
log({ event: 'terminal', output, passed: report.passed, accepted: report.accepted, errors: report.errors,
  residual: report.residual, coverage: report.coverage, failures: report.failures, seconds: report.seconds });
if (!report.passed) process.exitCode = 1;
