// SPDX-License-Identifier: GPL-2.0-or-later
// Compare a saved source tree against this checkout on the same restart.
// node scripts/validation/benchmark-coupled-kernels.js BEFORE_ROOT CHECKPOINT OUT [MEMBER]
// The source tree must include src/ and an ES-module package.json.
import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import { pathToFileURL, fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { createHash } from 'node:crypto';
import { createCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';

const [beforeRoot, source, output, member] = process.argv.slice(2);
if (!beforeRoot || !source || !output) throw new Error('Supply baseline source root, checkpoint and unused output filename.');
if (fs.existsSync(output)) throw new Error('Output already exists.');
const bytes = fs.readFileSync(source), hash = bytes => createHash('sha256').update(bytes).digest('hex');
let checkpoint = JSON.parse(source.endsWith('.gz') ? gunzipSync(bytes) : bytes);
if (member) for (const key of member.split('.')) checkpoint = checkpoint[key];
const r = checkpoint.restart;
if (!r) throw new Error('Select a coupled checkpoint with restart data.');
const options = { ...r.options, initialEuler: r.initialEuler, initialBL: r.initialBL };
const { createCoupledStreamtubeBody: createBefore } = await import(pathToFileURL(path.resolve(beforeRoot, 'src/euler/streamtube-coupled.js')));
const before = createBefore(structuredClone(r.input), structuredClone(options));
const after = createCoupledStreamtubeBody(structuredClone(r.input), structuredClone(options));
assert.deepEqual(after.initial, before.initial);
const samples = [], checks = {};
for (const kernel of ['evaluate', 'admissibleValue', 'jacobian']) {
  const run = system => kernel === 'evaluate'
    ? system.euler.evaluate(system.initial.subarray(0, system.ne))
    : system[kernel](system.initial);
  const serializable = value => JSON.parse(JSON.stringify(value));
  const a = serializable(run(before)), b = serializable(run(after));
  assert.ok(a && b, `${kernel}: checkpoint must be admissible`);
  assert.deepEqual(b, a, `${kernel}: values, geometry, diagnostics and sparse entries must match exactly`);
  checks[kernel] = hash(JSON.stringify(b));
  const repeats = kernel === 'jacobian' ? 3 : 30;
  for (const [label, system] of [['before', before], ['after', after], ['after', after], ['before', before]]) {
    for (let i = 0; i < 3; i++) run(system);
    const start = performance.now();
    for (let i = 0; i < repeats; i++) run(system);
    const sample = { kernel, label, repeats, milliseconds: performance.now() - start };
    samples.push(sample); console.log(JSON.stringify(sample));
  }
}
const sourceHashes = root => Object.fromEntries(fs.readdirSync(path.join(root, 'src'), { recursive: true })
  .map(p => path.join('src', p)).filter(p => fs.statSync(path.join(root, p)).isFile()).sort()
  .map(p => [p, hash(fs.readFileSync(path.join(root, p)))]));
fs.writeFileSync(output, JSON.stringify({ runtime: process.version, source, member, checkpointHash: hash(bytes),
  scope: 'Warm ABBA kernel timings, excluding construction and equality checks. No cold-start or end-to-end speed claim.',
  beforeHashes: sourceHashes(beforeRoot), afterHashes: sourceHashes(fileURLToPath(new URL('../..', import.meta.url))),
  checks, samples }, null, 2), { flag: 'wx' });
