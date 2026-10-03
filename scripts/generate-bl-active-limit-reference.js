// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import { mkdtemp } from 'node:fs/promises';
import { join } from 'node:path';
import { tmpdir } from 'node:os';
import { spawnSync } from 'node:child_process';
import { buildBLReference } from './reference/build-bl.js';
import { sha256 } from './validation/provenance.js';

const input = 'tests/fixtures/bl-active-limit-cases.json', cases = JSON.parse(fs.readFileSync(input)).cases;
const directory = await mkdtemp(join(tmpdir(), 'mses-bl-active-limits-')), native = await buildBLReference(directory, 'integral');
const lines = [cases.length];
for (const c of cases) {
  lines.push([c.input.regime === 'wake' ? 3 : 1, c.parameters.reynolds, c.parameters.mach, 0, c.parameters.ncrit ?? 9, 1e10, 1].join(' '));
  for (const s of [c.input.upstream, c.input.downstream]) lines.push([s.s, s.aux, s.aux, s.theta, s.deltaStar, 0, s.ue].join(' '));
  lines.push('0.03 0.003 0.005');
  c.expected = { residual: [], upstream: [[], [], []], downstream: [[], [], []] };
}
const r = spawnSync(native.executable, [], { input: lines.join('\n') + '\n', encoding: 'utf8', timeout: 10000 });
if (r.error || r.status !== 0) throw new Error(`Native limit oracle failed: ${r.error?.message ?? r.stderr}`);
for (const line of r.stdout.trim().split('\n')) {
  const [tag, ...words] = line.trim().split(/\s+/), [id, row, ...v] = words.map(Number);
  if (tag === 'R') cases[id - 1].expected.residual[row - 1] = v[0];
  if (tag === 'J') { cases[id - 1].expected.upstream[row - 1][v[0] - 1] = v[1]; cases[id - 1].expected.downstream[row - 1][v[0] - 1] = v[2]; }
}
if (cases.some(c => c.expected.residual.length !== 3 || c.expected.upstream.some(row => row.length !== 5)
  || c.expected.downstream.some(row => row.length !== 5) || !Object.values(c.expected).flat(2).every(Number.isFinite))) throw new Error('Incomplete native active-limit output.');
const files = [...native.files, input, 'third_party/Xfoil/src/XBL.INC', 'third_party/Xfoil/src/BLPAR.INC', 'scripts/reference/build-bl.js', 'scripts/generate-bl-active-limit-reference.js'];
const report = { provenance: { compiler: native.compiler, flags: native.flags, sourceHashes: Object.fromEntries(files.map(p => [p, sha256(p)])) }, cases, stdout: r.stdout };
fs.writeFileSync('tests/fixtures/fortran/bl-active-limits.json', JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ directory, cases: cases.length, compiler: native.compiler }));
