// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import { mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { spawnSync } from 'node:child_process';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { sha256 } from './validation/provenance.js';

const inputFile = 'docs/accepted-physical-domain-coupled-continuation.json', seed = JSON.parse(fs.readFileSync(inputFile)).restart;
const system = createCoupledStreamtubeBody(seed.input, { ...seed.options, initialEuler: seed.initialEuler, initialBL: Float64Array.from(seed.initialBL) });
const states = system.evaluate(system.initial).layers.states;
const intervals = [124, 125, 126, 127, 128].map(id => ({ id, parameters: system.bl.kernel.parameters,
  input: { upstream: states[id - 1], downstream: states[id], regime: 'wake' }, expected: { residual: [], upstream: [[], [], []], downstream: [[], [], []] } }));
const primitives = [1.02, 1.2, 1.99, 2.6, 4.5, 8].flatMap(hk => [5, 150, 1000].map(rt => ({ hk, rt })));
const directory = await mkdtemp(join(tmpdir(), 'mses-wake-jacobian-native-')), compiler = process.env.FC ?? 'gfortran';
const flags = ['-O2', '-std=legacy', '-fdefault-real-8', '-fallow-argument-mismatch', '-ffunction-sections', '-fdata-sections'];
const source = resolve('third_party/Xfoil/src');
const run = (executable, args, input, timeout = 60000) => {
  const r = spawnSync(executable, args, { input, encoding: 'utf8', timeout });
  if (r.error || r.status !== 0) throw new Error(`${executable} failed: ${r.error?.message ?? r.stderr}`);
  return r.stdout;
};
for (const name of ['xblsys', 'xbl']) run(compiler, [...flags, `-I${source}`, '-c', `${source}/${name}.f`, '-o', join(directory, `${name}.o`)]);
const execute = (driver, input) => {
  const executable = join(directory, driver);
  run(compiler, [...flags, `-I${source}`, '-Wl,--gc-sections', `scripts/reference/${driver}.f`,
    ...['xblsys', 'xbl'].map(n => join(directory, `${n}.o`)), '-o', executable]);
  return run(executable, [], input, 10000);
};
const primitiveStdout = execute('wake-dissipation', [primitives.length, ...primitives.map(c => `${c.hk} ${c.rt}`)].join('\n') + '\n');
for (const line of primitiveStdout.trim().split('\n')) {
  const [id, di, diHk, diRt, ...differences] = line.trim().split(/\s+/).map(Number);
  primitives[id - 1].expected = { di, diHk, diRt, differences };
}
const lines = [intervals.length];
for (const c of intervals) {
  lines.push([3, c.parameters.reynolds, c.parameters.mach, 0, c.parameters.ncrit, 1e10, 1].join(' '));
  for (const s of [c.input.upstream, c.input.downstream]) lines.push([s.s, s.aux, s.aux, s.theta, s.deltaStar, 0, s.ue].join(' '));
  lines.push('0.03 0.003 0.005');
}
const intervalStdout = execute('integral', lines.join('\n') + '\n');
for (const line of intervalStdout.trim().split('\n')) {
  const [tag, ...parts] = line.trim().split(/\s+/), [id, row, ...v] = parts.map(Number);
  if (tag === 'R') intervals[id - 1].expected.residual[row - 1] = v[0];
  if (tag === 'J') { intervals[id - 1].expected.upstream[row - 1][v[0] - 1] = v[1]; intervals[id - 1].expected.downstream[row - 1][v[0] - 1] = v[2]; }
}
if (primitives.some(c => !c.expected || !Object.values(c.expected).flat().every(Number.isFinite))
  || intervals.some(c => c.expected.residual.length !== 3 || c.expected.upstream.some(r => r.length !== 5)
    || c.expected.downstream.some(r => r.length !== 5) || !Object.values(c.expected).flat(2).every(Number.isFinite))) throw new Error('Incomplete wake oracle output.');
const files = ['scripts/generate-wake-jacobian-reference.js', 'scripts/reference/wake-dissipation.f', 'scripts/reference/integral.f',
  'third_party/Xfoil/src/xblsys.f', 'third_party/Xfoil/src/xbl.f', 'third_party/Xfoil/src/XBL.INC', 'third_party/Xfoil/src/BLPAR.INC'];
const report = { provenance: { compiler: run(compiler, ['--version']).split('\n')[0], flags,
  input: { path: inputFile, sha256: sha256(inputFile) }, sourceHashes: Object.fromEntries(files.map(p => [p, sha256(p)])) }, primitives, intervals, primitiveStdout, intervalStdout };
fs.writeFileSync('tests/fixtures/fortran/wake-jacobian.json', JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ directory, primitives: primitives.length, intervals: intervals.length, compiler: report.provenance.compiler }));
