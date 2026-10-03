// SPDX-License-Identifier: GPL-2.0-or-later
// Original native transition at exactly the TE, compared to the laminar
// interval at identical endpoint states. No airfoil solve or fitted output.
import fs from 'node:fs';
import { mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { spawnSync } from 'node:child_process';
import assert from 'node:assert/strict';
import { buildBLReference } from './reference/build-bl.js';
import { sha256 } from './validation/provenance.js';
const directory = await mkdtemp(join(tmpdir(), 'mses-terminal-transition-'));
const check = await buildBLReference(directory, 'automatic-transition'), integral = await buildBLReference(directory, 'integral');
const run = (executable, lines) => {
  const r = spawnSync(executable, [], { input: lines.join('\n') + '\n', encoding: 'utf8', timeout: 10000 });
  assert.equal(r.status, 0, r.error?.message ?? r.stderr); return r.stdout.trim().split('\n').map(l => l.trim().split(/\s+/));
};
const cases = [];
for (const mach of [0, .2, .3]) {
  const parameters = { reynolds: 1e5, mach, ncrit: 9 }, tripS = 1;
  const upstream = { s: .8, ue: 1, theta: .002, deltaStar: .005, aux: .7 };
  const downstream = { s: 1, ue: .98, theta: .0022, deltaStar: .0057, aux: .03 };
  const primitive = s => [s.s, s.aux, s.theta, s.deltaStar, s.ue].join(' ');
  const raw = run(check.executable, [1, [parameters.reynolds, mach, 9, tripS].join(' '), primitive(upstream), primitive(downstream)]);
  const checked = raw.find(r => r[0] === 'CHECK'), shear = Number(raw.find(r => r[0] === 'SHEAR')[2]);
  assert.equal(checked[2], 'T'); assert.equal(checked[3], 'T'); assert.equal(Number(checked[4]), tripS);
  const amplification = Number(checked[5]); assert.ok(amplification < 9 && shear > 0);
  downstream.aux = shear;
  const lines = [2];
  for (const type of [4, 1]) {
    lines.push([type, parameters.reynolds, mach, 0, 9, tripS, 1].join(' '));
    lines.push([upstream.s, upstream.aux, upstream.aux, upstream.theta, upstream.deltaStar, 0, upstream.ue].join(' '));
    lines.push([downstream.s, type === 4 ? upstream.aux : amplification, type === 4 ? shear : amplification,
      downstream.theta, downstream.deltaStar, 0, downstream.ue].join(' '));
    lines.push('0 0 0');
  }
  const output = run(integral.executable, lines), residuals = [[], []];
  output.filter(r => r[0] === 'R').forEach(r => { residuals[Number(r[1]) - 1][Number(r[2]) - 1] = Number(r[3]); });
  assert.ok(residuals.flat().every(Number.isFinite)); assert.ok(Math.abs(residuals[0][0]) < 1e-12);
  for (const row of [1, 2]) assert.ok(Math.abs(residuals[0][row] - residuals[1][row]) < 1e-12);
  cases.push({ parameters, upstream, downstream, tripS, expected: { shear, amplification, transition: residuals[0], laminar: residuals[1] }, stdout: { check: raw, integral: output } });
}
const paths = [...new Set([...check.files, ...integral.files, 'third_party/Xfoil/src/XBL.INC', 'third_party/Xfoil/src/BLPAR.INC',
  'scripts/reference/build-bl.js', 'scripts/generate-terminal-transition-reference.js'])];
const report = { provenance: { date: new Date().toISOString(), compiler: check.compiler, flags: check.flags,
  scope: 'Executed unmodified XFOIL TRCHEK, BLSYS, BLVAR. A zero-length turbulent part at the trailing edge retains laminar momentum/energy equations and supplies turbulent shear for wake matching.',
  sha256: Object.fromEntries(paths.map(p => [p, sha256(p)])) }, cases };
const output = 'tests/fixtures/fortran/terminal-transition.json';
fs.writeFileSync(output, JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ output, compiler: check.compiler, cases: cases.length }));
