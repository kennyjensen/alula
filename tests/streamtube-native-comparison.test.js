import test from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { join } from 'node:path';
import { tmpdir } from 'node:os';

test('native Cp validation includes the stagnation node omitted from the BL station list', () => {
  const directory = mkdtempSync(join(tmpdir(), 'mses-native-pressure-'));
  try {
    const path = join(directory, 'comparison.json');
    const r = spawnSync(process.execPath, ['scripts/compare-streamtube-matched-reference.js',
      'docs/coupled-nested-isentropic-root.json', path], { encoding: 'utf8', timeout: 10000 });
    assert.equal(r.status, 1, r.stderr); const report = JSON.parse(readFileSync(path));
    assert.equal(report.error, undefined); assert.equal(report.checkedResidual < 1e-10, true);
    assert.equal(report.pressureProfiles.length, 2);
    for (const branch of report.pressureProfiles) {
      const bl = report.profiles.find(b => b.side === branch.side);
      assert.equal(branch.points.length, bl.points.length + 1); assert.equal(branch.points[0].k, 0);
    }
    // The previously reported BL-only Cp maximum would pass. Full surface
    // coverage must expose the independently unphysical stagnation pressure.
    assert.ok(report.blStationCpMax < report.gates.cpMax);
    assert.equal(report.worst.cp.k, 0); assert.ok(report.errors.cpMax > .08);
    assert.ok(report.failures.some(s => s.startsWith('cpMax:')));
    assert.equal(report.gates.cpMax, .04); assert.equal(report.uncoveredStations.length, 1);
    assert.ok(report.stagnationPressure.points.every(p => p.cpExcess > .08));
    assert.equal(report.changedSources.length, 0);
  } finally { rmSync(directory, { recursive: true, force: true }); }
});

test('native leading-edge refinement uses exact contour parameters, closes coverage and retains the failed BL gates', () => {
  const directory = mkdtempSync(join(tmpdir(), 'mses-native-leading-'));
  try {
    const path = join(directory, 'comparison.json');
    const r = spawnSync(process.execPath, ['scripts/compare-streamtube-matched-reference.js',
      'docs/coupled-stagnation-normal-refined-root.json', path,
      '--native=docs/native-leading-reference-refinement.json', '--native-case=leading-4'], { encoding: 'utf8', timeout: 10000 });
    assert.equal(r.status, 1, r.stderr); const report = JSON.parse(readFileSync(path));
    assert.equal(report.error, undefined); assert.equal(report.uncoveredStations.length, 0);
    assert.ok(report.errors.cpMax < report.gates.cpMax);
    assert.ok(report.errors.ueMax > report.gates.ueMax);
    assert.ok(report.errors.thetaRelativeMax > report.gates.thetaRelativeMax);
    assert.ok(report.errors.deltaStarRelativeMax > report.gates.deltaStarRelativeMax);
    assert.equal(report.failures.length, 3); assert.equal(report.changedSources.length, 0);
  } finally { rmSync(directory, { recursive: true, force: true }); }
});
