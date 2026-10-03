// SPDX-License-Identifier: GPL-2.0-or-later
// Compile and call only the unchanged original XICALC routine on three tiny
// prescribed geometries. No panel solution, BL interval or nonlinear solve.
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { createXfoilDeadAirGap } from '../src/viscous/xfoil-dead-air-gap.js';

const source = 'third_party/Xfoil/src/xpanel.f', fixture = 'tests/fixtures/fortran/dead-air-gap.json';
const output = 'docs/nlr-finite-base/dead-air-gap-native.json';
assert.ok(!fs.existsSync(fixture) && !fs.existsSync(output), 'Preserve original native evidence.');
const tracked = [source, 'third_party/Xfoil/src/XFOIL.INC', 'third_party/Xfoil/src/PINDEX.INC',
  'src/viscous/xfoil-dead-air-gap.js', 'scripts/check-xfoil-dead-air-gap-native.js'];
const sha = data => createHash('sha256').update(data).digest('hex');
const hashes = () => Object.fromEntries(tracked.map(p => [p, sha(fs.readFileSync(p))]));
const sourceHashes = hashes(), start = performance.now(), directory = fs.mkdtempSync(path.join(os.tmpdir(), 'mses-dead-air-gap-'));
const text = fs.readFileSync(source, 'utf8'), begin = text.indexOf('      SUBROUTINE XICALC'), end = text.indexOf('      SUBROUTINE UICALC', begin);
assert.ok(begin > 0 && end > begin);
const extracted = text.slice(begin, end); fs.writeFileSync(path.join(directory, 'xicalc.f'), extracted);
const driver = `      PROGRAM GAPCONTROL
      INCLUDE 'XFOIL.INC'
      INTEGER NC, IC, K, M
      REAL DIST(IWX)
      READ(*,*) NC
      DO 100 IC=1,NC
        READ(*,*) ANTE,XP(1),YP(1),XP(4),YP(4),NW
        READ(*,*) (DIST(K),K=1,NW)
        N=4
        SST=0.3
        DO 10 K=1,N
          S(K)=0.2*(K-1)
   10   CONTINUE
        SHARP=.FALSE.
        IBLTE(1)=2
        IBLTE(2)=2
        IPAN(2,1)=1
        IPAN(2,2)=N
        NBL(1)=2+NW
        NBL(2)=2+NW
        DO 20 K=1,NW
          X(N+K)=DIST(K)
          Y(N+K)=0.0
          IPAN(2+K,1)=N+K
          IPAN(2+K,2)=N+K
   20   CONTINUE
        CALL XICALC
        DO 30 K=1,NW
          WRITE(*,'(I3,1X,I3,1X,ES26.17,1X,ES26.17)')
     &      IC,K,XSSI(2+K,2)-XSSI(2,2),WGAP(K)
   30   CONTINUE
  100 CONTINUE
      END
`;
fs.writeFileSync(path.join(directory, 'driver.f'), driver);
const compiler = process.env.FC || (fs.existsSync('/tmp/mses-fortran-toolchain/gfortran') ? '/tmp/mses-fortran-toolchain/gfortran' : 'gfortran');
const executable = path.join(directory, 'gap-control');
const flags = ['-O2', '-std=legacy', '-fdefault-real-8', `-I${path.resolve('third_party/Xfoil/src')}`];
const c = 1 / Math.sqrt(5), cases = [
  { name: 'parallel-endpoint-derivatives', input: { normalGap: .002, upperDerivative: { x: -1, y: 0 }, lowerDerivative: { x: 1, y: 0 } } },
  { name: 'negative-source-slope-clamp', input: { normalGap: .00089, upperDerivative: { x: -1, y: 0 }, lowerDerivative: { x: c, y: 2 * c } } },
  { name: 'positive-source-slope-clamp', input: { normalGap: .001150281917618465, upperDerivative: { x: -1, y: 0 }, lowerDerivative: { x: c, y: -2 * c } } },
];
const fractions = [0, .04, .2, .6, .95, 1, 1.1];
const report = { date: new Date().toISOString(), inProgress: true, passed: false, sourceHashes, directory,
  scope: 'Three prescribed ANTE/endpoint-derivative XICALC controls, 21 scalar widths. Exact extracted original routine; no TECALC, spline, airfoil, BL residual, Euler or nonlinear solve.',
  source: { path: source, firstLine: text.slice(0, begin).split('\n').length,
    extractedSHA256: sha(extracted), exactByteExtraction: true },
  native: { compiler, flags, driver, driverSHA256: sha(driver) }, cases, maximumGapError: null };
const save = () => fs.writeFileSync(output, JSON.stringify(report, null, 2) + '\n'); save();
try {
  const built = spawnSync(compiler, [...flags, path.join(directory, 'driver.f'), path.join(directory, 'xicalc.f'), '-o', executable],
    { encoding: 'utf8', timeout: 20000 });
  fs.writeFileSync(path.join(directory, 'compiler.log'), (built.stdout ?? '') + (built.stderr ?? ''));
  assert.ok(built.status === 0, built.error?.message ?? built.stderr);
  report.native.version = spawnSync(compiler, ['--version'], { encoding: 'utf8' }).stdout.split('\n')[0];
  report.native.executableSHA256 = sha(fs.readFileSync(executable));
  const lines = [cases.length];
  for (const row of cases) {
    const g = row.input;
    row.requestedDistances = fractions.map(f => f * 2.5 * g.normalGap);
    lines.push([g.normalGap, g.upperDerivative.x, g.upperDerivative.y, g.lowerDerivative.x, g.lowerDerivative.y, fractions.length].join(' '));
    lines.push(row.requestedDistances.join(' ')); row.samples = [];
  }
  const stdin = lines.join('\n') + '\n'; fs.writeFileSync(path.join(directory, 'stdin.txt'), stdin);
  const run = spawnSync(executable, [], { input: stdin, encoding: 'utf8', timeout: 10000 });
  fs.writeFileSync(path.join(directory, 'stdout.log'), run.stdout ?? ''); fs.writeFileSync(path.join(directory, 'stderr.log'), run.stderr ?? '');
  assert.ok(run.status === 0, run.error?.message ?? run.stderr);
  report.native.stdout = run.stdout; report.native.stdoutSHA256 = sha(run.stdout);
  let maximumGapError = 0;
  for (const line of run.stdout.trim().split('\n')) {
    const [ic, k, distance, expectedGap] = line.trim().split(/\s+/).map(Number);
    assert.ok(Number.isInteger(ic) && ic >= 1 && ic <= cases.length && k >= 1 && k <= fractions.length
      && Number.isFinite(distance) && Number.isFinite(expectedGap));
    const row = cases[ic - 1], value = createXfoilDeadAirGap(row.input).at(distance), error = Math.abs(value.gap - expectedGap);
    maximumGapError = Math.max(maximumGapError, error);
    row.samples.push({ distance, expectedGap, javascript: value, error });
    assert.ok(error <= 2e-14 * row.input.normalGap, 'Original XICALC width parity failed.');
  }
  assert.ok(cases.every(c => c.samples.length === fractions.length));
  report.maximumGapError = maximumGapError; report.passed = true;
} catch (error) { report.error = { message: error.message, stack: error.stack }; }
report.sourceHashesAfter = hashes();
if (JSON.stringify(sourceHashes) !== JSON.stringify(report.sourceHashesAfter)) report.passed = false;
report.inProgress = false; report.seconds = (performance.now() - start) / 1000; save();
if (report.passed) fs.writeFileSync(fixture, JSON.stringify({ passed: true, sourceHashes, source: report.source,
  native: report.native, cases }, null, 2) + '\n');
console.log(JSON.stringify({ passed: report.passed, error: report.error, output, fixture: report.passed ? fixture : null,
  cases: cases.length, samples: cases.reduce((n, c) => n + (c.samples?.length ?? 0), 0), maximumGapError: report.maximumGapError,
  seconds: report.seconds }));
if (!report.passed) process.exitCode = 1;
