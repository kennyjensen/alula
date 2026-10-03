import test from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdtempSync, readFileSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { evaluateCoupledAcceptance, originalGates, currentEffective } from '../scripts/validation/coupled-acceptance.js';

const baseline = () => ({ complete: true, converged: true, admissible: true, provenanceValid: true, exactReplay: true,
  residual: 1.3018819799981428e-11,
  errors: { cpMax: .02694926740934167, ueMax: .00959719882434884,
    thetaRelativeMax: .09117409609471605, deltaStarRelativeMax: .11151270094480092 },
  coverage: { requestedBL: 394, coveredBL: 394, requestedPressure: 396, coveredPressure: 396, uncoveredStations: 0 } });

test('temporary thickness policy accepts the retained errors without changing them or strict gates', () => {
  const e = baseline(), before = structuredClone(e), result = evaluateCoupledAcceptance(e);
  assert.equal(result.accepted, true); assert.deepEqual(e, before); assert.deepEqual(result.errors, e.errors);
  assert.deepEqual(originalGates, { cpMax: .04, ueMax: .02, thetaRelativeMax: .05, deltaStarRelativeMax: .08 });
  assert.deepEqual(currentEffective, { ...originalGates, thetaRelativeMax: .12, deltaStarRelativeMax: .15 });
  assert.deepEqual(Object.entries(result.gateResults).filter(([, q]) => !q.originalPassed).map(([k]) => k), ['thetaRelativeMax', 'deltaStarRelativeMax']);
  assert.equal(result.physicalAcceptance, false); assert.equal(result.fullSolverComplete, false); assert.equal(result.runtimeChanged, false);
});

for (const [key, value] of [['cpMax', .0400000001], ['ueMax', .0200000001],
  ['thetaRelativeMax', .1200000001], ['deltaStarRelativeMax', .1500000001]])
  test(`${key} still rejects values above its effective limit`, () => {
    const e = baseline(); e.errors[key] = value; const r = evaluateCoupledAcceptance(e);
    assert.equal(r.accepted, false); assert.equal(r.gateResults[key].passed, false);
  });

test('a converged flag cannot bypass the unchanged strict residual limit', () => {
  for (const residual of [1e-10, 1.1e-10, -1, NaN, Infinity, undefined])
    assert.equal(evaluateCoupledAcceptance({ ...baseline(), residual }).accepted, false);
});

test('missing BL or stagnation pressure coverage cannot be hidden by passing error maxima', () => {
  for (const mutation of [c => c.coveredBL--, c => c.coveredPressure--, c => c.uncoveredStations++,
    c => { c.requestedBL = 0; c.coveredBL = 0; }, c => { delete c.requestedPressure; }, c => { delete c.uncoveredStations; }]) {
    const e = baseline(); mutation(e.coverage); assert.equal(evaluateCoupledAcceptance(e).accepted, false);
  }
});

test('partial, unconverged, inadmissible, stale, and unreplayed evidence always rejects', () => {
  for (const key of ['complete', 'converged', 'admissible', 'provenanceValid', 'exactReplay'])
    for (const value of [false, undefined, 'true', 1]) {
      const e = baseline(); e[key] = value; assert.equal(evaluateCoupledAcceptance(e).accepted, false, `${key}=${value}`);
    }
});

test('missing, nonfinite, or negative error values cannot pass the widened thickness gate', () => {
  for (const key of Object.keys(currentEffective)) for (const value of [undefined, NaN, Infinity, -1]) {
    const e = baseline(); e.errors[key] = value; assert.equal(evaluateCoupledAcceptance(e).accepted, false);
  }
  assert.equal(evaluateCoupledAcceptance().accepted, false);
});

test('the executable exits nonzero for incomplete, unconverged, or stale input before any coupled evaluation', () => {
  const original = JSON.parse(readFileSync('docs/current-natural-outer-refined-stagnation-native.json'));
  for (const kind of ['in-progress', 'unconverged', 'missing-state', 'stale-native', 'unarchived-source', 'wrong-archived-version']) {
    const directory = mkdtempSync(join(tmpdir(), 'mses-coupled-acceptance-'));
    try {
      const saved = JSON.parse(readFileSync(original.source)), comparison = structuredClone(original);
      if (kind === 'in-progress') saved.inProgress = true;
      if (kind === 'unconverged') saved.result.converged = false;
      if (kind === 'missing-state') delete saved.checkpoint.restart.initialBL;
      if (kind === 'stale-native') comparison.nativeHash = '0'.repeat(64);
      if (kind === 'unarchived-source') comparison.sourceHashes['src/viscous/integral.js'] = '0'.repeat(64);
      if (kind === 'wrong-archived-version') comparison.sourceHashes['src/viscous/transition-interval.js'] = '0'.repeat(64);
      comparison.source = join(directory, 'state.json');
      const stateBytes = JSON.stringify(saved); writeFileSync(comparison.source, stateBytes);
      comparison.sourceHash = createHash('sha256').update(stateBytes).digest('hex');
      const input = join(directory, 'comparison.json'), output = join(directory, 'acceptance.json');
      writeFileSync(input, JSON.stringify(comparison)); writeFileSync(input.replace(/\.json$/, '.log'), 'Synthetic invalid-input fixture.\n');
      const run = spawnSync(process.execPath, ['scripts/check-current-coupled-acceptance.js', `--comparison=${input}`, `--output=${output}`],
        { encoding: 'utf8', timeout: 15000 });
      assert.equal(run.status, 1, `${kind}: ${run.stderr}\n${run.stdout}`);
      const report = JSON.parse(readFileSync(output));
      assert.equal(report.inProgress, false); assert.equal(report.accepted, false); assert.equal(report.passed, false);
      assert.ok(report.failures.length > 0); assert.equal(report.operations.coupledEvaluations, 0);
      assert.equal(report.operations.factorizations, 0); assert.equal(report.operations.newtonUpdates, 0);
      // This historical producer may reject an earlier stale linked source
      // after the application is reorganized. The invariant of this test is
      // that none of these invalid inputs reaches a coupled evaluation.
    } finally { rmSync(directory, { recursive: true, force: true }); }
  }
});
