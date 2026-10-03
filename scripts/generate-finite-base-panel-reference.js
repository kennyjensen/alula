// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import { mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, relative } from 'node:path';
import { spawnSync } from 'node:child_process';
import { naca4, transform } from '../src/geometry/airfoil.js';
import { buildReference } from './reference/build.js';
import { sha256 } from './validation/provenance.js';

const output = 'tests/fixtures/fortran/finite-base-panel.json';
if (fs.existsSync(output)) throw new Error('Preserve the original native fixture.');
const started = performance.now(), directory = await mkdtemp(join(tmpdir(), 'mses-finite-base-panel-'));
const built = await buildReference(directory, 'driver-finite-base-panel');
const surface = (gap, skew = 0) => naca4('0012', 40).map((p, i) => ({ x: p.x + (i > 20 ? skew * p.x ** 2 : 0),
  y: p.y + (i <= 20 ? 1 : -1) * .5 * gap * p.x ** 2 }));
const cases = [
  { name: 'symmetric-small-base', points: surface(.002), alpha: 4 },
  { name: 'thicker-skew-base', points: surface(.008, .003), alpha: 7 },
  { name: 'rotated-skew-base', points: transform(surface(.008, .003), { angle: 17 }), alpha: 24 },
];
const fields = [{ x: -.2, y: .13 }, { x: .4, y: .2 }, { x: .4, y: -.2 }, { x: 1.2, y: .12 }, { x: 1.2, y: -.12 }];
const report = { date: new Date().toISOString(), passed: false, physicalAcceptance: false,
  scope: 'Three tiny inviscid original-Fortran finite-TE controls. No viscous/Newton solve or NLR assembly solve.',
  provenance: { compiler: built.compiler, flags: built.flags,
    sourceHashes: Object.fromEntries([...built.files, 'scripts/reference/driver-finite-base-panel.f', 'scripts/reference/build.js',
      'scripts/generate-finite-base-panel-reference.js', 'src/geometry/airfoil.js'].map(p => [relative(process.cwd(), p), sha256(p)])) }, cases: [] };
for (const item of cases) {
  const queries = item.name.startsWith('rotated') ? transform(fields, { angle: 17 }) : fields;
  const input = `${item.points.length} ${item.alpha}\n${item.points.map(p => `${p.x} ${p.y}`).join('\n')}\n${queries.length}\n${queries.map(p => `${p.x} ${p.y}`).join('\n')}\n`;
  const result = spawnSync(built.executable, [], { input, encoding: 'utf8', timeout: 15000 });
  if (result.error || result.status !== 0) throw new Error(`Native control failed: ${result.error?.message ?? result.stderr}`);
  const rows = result.stdout.split(/\r?\n/).map(line => line.trim().split(/\s+/)), get = label => rows.filter(row => row[0] === label);
  const values = get('RESULT')[0]?.slice(1);
  if (values?.length !== 10) throw new Error('Missing complete native inviscid result.');
  const keys = ['cl', 'cm', 'pressureDrag', 'chord', 'normalGap', 'tangentialGap', 'gap', 'sourceStrength', 'vortexStrength'];
  const expected = Object.fromEntries(keys.map((key, i) => [key, Number(values[i])]));
  expected.sharp = values[9] === 'T';
  expected.nodes = get('NODE').map(row => ({ index: Number(row[1]), gamma: Number(row[2]), dxds: Number(row[3]), dyds: Number(row[4]) }));
  expected.fields = get('FIELD').map(row => ({ index: Number(row[1]), x: Number(row[2]), y: Number(row[3]),
    psi: Number(row[4]), u: Number(row[5]), v: Number(row[6]) }));
  if (expected.sharp || expected.nodes.length !== item.points.length || expected.fields.length !== queries.length)
    throw new Error('Native finite-base output shape failed.');
  report.cases.push({ ...item, queries, expected });
  console.log(JSON.stringify({ case: item.name, nodes: item.points.length, expected: { ...expected, nodes: undefined, fields: undefined } }));
}
report.passed = true; report.seconds = (performance.now() - started) / 1000;
fs.writeFileSync(output, JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({ passed: true, output, seconds: report.seconds, cases: report.cases.length, nativeViscousSolves: 0, nlrAssemblySolves: 0 }));
