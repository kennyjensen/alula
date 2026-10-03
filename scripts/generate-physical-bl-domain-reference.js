// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import { mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { spawnSync } from 'node:child_process';
import { sha256 } from './validation/provenance.js';

const cases = [];
for (const mach of [0, .2, .75]) for (const gamma of [1.4, 1.3]) for (const h of [1.01, 1.2, 2.6]) {
  const theta = .001, wakeGap = h === 1.2 ? .0003 : 0;
  cases.push({ parameters: { mach, gamma }, station: { theta, deltaStar: theta * h + wakeGap, wakeGap, ue: 1.3 } });
}
const directory = await mkdtemp(join(tmpdir(), 'mses-bl-domain-native-'));
const executable = join(directory, 'bl-domain'), compiler = process.env.FC ?? 'gfortran';
// Compile the original file intact. Section garbage collection discards
// unrelated routines; no full XFOIL run or iterative profile is needed.
const files = ['scripts/reference/physical-bl-domain.f', 'third_party/Xfoil/src/xblsys.f'];
const flags = ['-O2', '-std=legacy', '-fdefault-real-8', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections'];
const build = spawnSync(compiler, [...flags, `-I${resolve('third_party/Xfoil/src')}`, ...files, '-o', executable], { encoding: 'utf8', timeout: 60000 });
if (build.error || build.status !== 0) throw new Error(`BL domain oracle build failed: ${build.error?.message ?? build.stderr}`);
const lines = [cases.length, ...cases.map(({ parameters: p, station: s }) => [p.mach, p.gamma, s.theta, s.deltaStar, s.wakeGap, s.ue].join(' '))];
const run = spawnSync(executable, [], { input: lines.join('\n') + '\n', encoding: 'utf8', timeout: 10000 });
if (run.error || run.status !== 0) throw new Error(`BL domain oracle failed: ${run.error?.message ?? run.stderr}`);
for (const line of run.stdout.trim().split('\n')) {
  const [id, machSquared, machSquaredUe, hk, ...hkGradient] = line.trim().split(/\s+/).map(Number);
  cases[id - 1].expected = { machSquared, machSquaredUe, hk, hkGradient };
}
if (cases.some(c => !c.expected || !Object.values(c.expected).flat().every(Number.isFinite))) throw new Error('Incomplete native domain fixture.');
const sourceFiles = [...files, 'third_party/Xfoil/src/XBL.INC', 'scripts/generate-physical-bl-domain-reference.js'];
const report = { provenance: { compiler: spawnSync(compiler, ['--version'], { encoding: 'utf8' }).stdout.split('\n')[0], flags,
  sourceHashes: Object.fromEntries(sourceFiles.map(p => [p, sha256(p)])) }, cases, stdout: run.stdout };
fs.writeFileSync('tests/fixtures/fortran/physical-bl-domain.json', JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ directory, compiler: report.provenance.compiler, cases: cases.length }));
