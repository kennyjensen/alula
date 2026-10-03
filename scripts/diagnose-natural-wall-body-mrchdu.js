// SPDX-License-Identifier: GPL-2.0-or-later
// Native verification of complete per-body local preparation; no flow solve.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { spawnSync } from 'node:child_process';
import assert from 'node:assert/strict';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { relaxXfoilBody } from './validation/mrchdu-body.js';
import { numericalSourceHashes, changedSources, sha256 } from './validation/provenance.js';

const output = 'docs/current-natural-wall-body-mrchdu.json', fixture = 'tests/fixtures/fortran/mrchdu-body.json';
assert.ok(!fs.existsSync(output) && !fs.existsSync(fixture));
const started = performance.now(), parents = [], cases = [], replays = [];
const sourceHashes = numericalSourceHashes(['scripts/diagnose-natural-wall-body-mrchdu.js', 'scripts/validation/mrchdu-body.js', 'scripts/reference/mrchdu-body.f']);
const read = (p, historical = false) => {
  const r = JSON.parse(fs.readFileSync(p)), changes = changedSources(r.sourceHashes);
  if (!historical) assert.deepEqual(changes, []);
  else assert.deepEqual(changes, [
    'src/euler/streamtube-boundary-layers.js', 'src/euler/streamtube-coupled-initializer.js',
    'src/euler/tests/streamtube-coupled-redistribution.js', 'src/euler/streamtube-coupled-refinement.js', 'src/viscous/transition-selection.js',
  ]);
  parents.push({ path: p, sha256: sha256(p), historicalSourceChanges: changes }); return r;
};
const add = (system, value, name, expected) => {
  assert.deepEqual(value.families, expected);
  replays.push({ name, unknowns: system.n, families: value.families, exactFamilyReplay: true });
  for (const w of system.bl.wakes) {
    const surfaces = system.bl.surfaces.filter(s => s.body === w.body);
    assert.deepEqual(surfaces.map(s => s.side), ['upper', 'lower']);
    cases.push({ name: `${name}/body-${w.body}`, body: w.body, parameters: system.bl.kernel.parameters,
      input: { surfaces: surfaces.map(s => s.ids.map(id => value.layers.states[id])),
        phases: surfaces.map(s => s.transition), wake: w.ids.map(id => value.layers.states[id]) } });
  }
};
const make = f => createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
const events = read('docs/current-natural-wall-newton-events.json'), s = make(events.checkpoint.restart);
add(s, s.evaluate(s.initial), 'before-newton', events.initial.families);
const ray = events.trials.find(t => t.originalAcceptedStepComparison);
assert.equal(ray.guardPassed, true); assert.equal(ray.event.changed, false);
const point = s.initial.map((v, i) => v + (i < s.euler.layout.densityCount
  ? Math.log1p(ray.step * events.linear.direction[i]) : ray.step * events.linear.direction[i]));
add(s, s.evaluate(point), 'first-accepted-ray', ray.beforeMaintenance.families);
for (const [p, key, name] of [
  ['docs/current-natural-graded-wall-first.json', 'restart', 'eight-update-endpoint'],
  ['docs/current-natural-wall-predicted-descent.json', 'restart', 'Euler-predictor-and-descent'],
  ['docs/current-natural-wall-mrchdu-euler-equilibrated.json', 'seed', 'surface-MRCHDU-and-converged-Euler'],
  ['docs/current-quad-automatic-mrchue-public.json', 'restart', 'default-two-element-root'],
]) {
  const r = read(p, name === 'default-two-element-root'); assert.equal(r.inProgress, false);
  const system = make(r[key]), value = system.evaluate(system.initial);
  add(system, value, name, key === 'seed' ? r.initial.families : r.result.families);
  assert.deepEqual(value.outer.nodes, r[key].initialEuler.nodes);
}
const nativeFiles = ['scripts/reference/mrchdu-body.f', ...['xblsys.f', 'xbl.f', 'xsolve.f', 'spline.f'].map(p => `third_party/Xfoil/src/${p}`)];
const nativeHashes = Object.fromEntries([...nativeFiles, ...['XFOIL.INC', 'XBL.INC', 'BLPAR.INC', 'PINDEX.INC', 'xpanel.f'].map(p => `third_party/Xfoil/src/${p}`)].map(p => [p, sha256(p)]));
const compiler = process.env.FC ?? 'gfortran', folder = fs.mkdtempSync(path.join(os.tmpdir(), 'mses-body-mrchdu-'));
const executable = path.join(folder, 'mrchdu-body');
const flags = ['-O2', '-std=legacy', '-fdefault-real-8', '-fallow-argument-mismatch', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections'];
const report = { date: new Date().toISOString(), physicalAcceptance: false, runtimeChanged: false, inProgress: true,
  sourceHashes, parents, replays, cases, scope: 'Unmodified original MRCHDU on differing upper/lower surfaces and a merged sharp-TE wake, using native IBLPAN/XICALC indexing. Prescribed retained states only; no coupled solve or physical acceptance. The historical default root is used only after exact current-equation family and grid replay; its changed source paths remain explicit.' };
const save = () => fs.writeFileSync(output, JSON.stringify(report) + '\n'); save();
try {
  const built = spawnSync(compiler, [...flags, `-I${path.resolve('third_party/Xfoil/src')}`, ...nativeFiles, '-o', executable], { encoding: 'utf8', timeout: 60000 });
  assert.equal(built.status, 0, built.error?.message ?? built.stderr);
  const compilerVersion = spawnSync(compiler, ['--version'], { encoding: 'utf8' }).stdout.split('\n')[0];
  const lines = [cases.length];
  for (const c of cases) {
    assert.equal(c.parameters.gamma, 1.4); assert.equal(c.parameters.velocityConvention, 'physical');
    const { surfaces, wake, phases } = c.input;
    lines.push([c.parameters.reynolds, c.parameters.mach, c.parameters.ncrit, ...surfaces.map(s => s.length), wake.length, ...phases].join(' '));
    [...surfaces.flat(), ...wake].forEach(p => lines.push([p.s, p.ue, p.aux, p.theta, p.deltaStar].join(' ')));
  }
  const run = spawnSync(executable, [], { input: lines.join('\n') + '\n', encoding: 'utf8', timeout: 20000, maxBuffer: 8 * 1024 * 1024 });
  assert.equal(run.status, 0, run.error?.message ?? run.stderr); report.nativeStdout = run.stdout;
  let current;
  for (const line of run.stdout.trim().split('\n')) {
    const fields = line.trim().split(/\s+/), kind = fields[0];
    if (kind === 'CASE') { current = cases[Number(fields[1]) - 1]; current.nativeWarnings = []; }
    if (/Convergence failed/.test(line)) current.nativeWarnings.push(line.trim());
    if (!['STATION', 'TRANSITION'].includes(kind)) continue;
    const c = cases[Number(fields[1]) - 1]; c.native ??= { surfaces: [{ states: [] }, { states: [] }], wake: [] };
    const side = Number(fields[2]) - 1;
    if (kind === 'TRANSITION') Object.assign(c.native.surfaces[side], { transition: Number(fields[3]), s: Number(fields[4]), forced: fields[5] === 'T' });
    else {
      const i = Number(fields[3]), [s, ue, aux, theta, deltaStar] = fields.slice(4).map(Number), point = { s, ue, aux, theta, deltaStar };
      if (side === 1 && i >= c.input.surfaces[1].length) c.native.wake[i - c.input.surfaces[1].length] = point;
      else c.native.surfaces[side].states[i] = point;
    }
  }
  let stations = 0, maximumError = 0, changedIntervals = 0;
  for (const c of cases) {
    const before = structuredClone(c.input), result = relaxXfoilBody(c.input, c.parameters);
    assert.deepEqual(c.input, before); c.javascript = result;
    const warningSites = lines => lines.map(l => l.match(/failed at\s*(\d+)\s+side\s*(\d+)/)?.slice(1));
    assert.deepEqual(warningSites(result.localConvergenceWarnings), warningSites(c.nativeWarnings));
    let worst = { error: 0 }, maxWakeThetaChange = 0;
    c.native.surfaces.forEach((expected, side) => {
      assert.equal(result.surfaces[side].transition, expected.transition); assert.equal(result.surfaces[side].forced, expected.forced);
      assert.ok(Math.abs(result.surfaces[side].s - expected.s) < 1e-12);
      if (expected.transition !== c.input.phases[side]) changedIntervals++;
    });
    for (const [region, actual, expected] of [
      ['upper', result.surfaces[0].states, c.native.surfaces[0].states],
      ['lower', result.surfaces[1].states, c.native.surfaces[1].states], ['wake', result.wake, c.native.wake],
    ]) {
      assert.equal(actual.length, expected.length);
      actual.forEach((p, i) => {
        assert.equal(p.s, expected[i].s);
        for (const key of ['ue', 'aux', 'theta', 'deltaStar']) {
          const a = p[key], b = expected[i][key], error = Math.abs(a - b) / (key === 'aux' ? Math.max(.01, Math.abs(b)) : Math.abs(b));
          assert.ok(Number.isFinite(error)); if (error > worst.error) worst = { region, i, key, error, actual: a, expected: b };
        }
        stations++;
      });
    }
    result.wake.forEach((p, i) => { maxWakeThetaChange = Math.max(maxWakeThetaChange, Math.abs(p.theta / c.input.wake[i].theta - 1)); });
    const [a, b] = result.surfaces.map(s => s.states.at(-1)), w = result.wake[0];
    const teDefects = { theta: w.theta - a.theta - b.theta, deltaStar: w.deltaStar - a.deltaStar - b.deltaStar,
      shearMomentum: w.aux * w.theta - a.aux * a.theta - b.aux * b.theta };
    c.comparison = { worst, beforePhase: c.input.phases, afterPhase: result.surfaces.map(s => s.transition),
      warningCount: result.localConvergenceWarnings.length, maxWakeThetaRelativeChange: maxWakeThetaChange, teDefects };
    maximumError = Math.max(maximumError, worst.error); assert.ok(worst.error < 1e-10, `${c.name}: ${JSON.stringify(worst)}`);
    console.log(JSON.stringify({ name: c.name, ...c.comparison }));
  }
  report.summary = { cases: cases.length, stations, maximumError, changedIntervals,
    localWarnings: cases.reduce((sum, c) => sum + c.javascript.localConvergenceWarnings.length, 0) };
  report.provenance = { compiler, compilerVersion, flags, sha256: nativeHashes };
  report.passed = true;
  fs.writeFileSync(fixture, JSON.stringify({ scope: report.scope, provenance: report.provenance,
    parent: { path: output, sourceHashes }, cases: cases.map(c => ({ name: c.name, parameters: c.parameters, input: c.input,
      native: c.native, nativeWarnings: c.nativeWarnings })) }) + '\n');
} catch (error) { report.passed = false; report.error = error.stack; }
finally { fs.rmSync(folder, { recursive: true, force: true }); }
assert.deepEqual(changedSources(sourceHashes), []); assert.deepEqual(changedSources(nativeHashes), []);
parents.forEach(p => assert.equal(sha256(p.path), p.sha256));
report.inProgress = false; report.seconds = (performance.now() - started) / 1000; save();
console.log(JSON.stringify({ output, passed: report.passed, seconds: report.seconds, summary: report.summary, error: report.error }));
if (!report.passed) process.exitCode = 1;
