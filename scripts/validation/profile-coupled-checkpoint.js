// SPDX-License-Identifier: GPL-2.0-or-later
// Profile a fixed number of real coupled updates from an unchanged checkpoint.
// node scripts/validation/profile-coupled-checkpoint.js CHECKPOINT OUT [UPDATES=3] [MEMBER]
// Example: tests/fixtures/rae64x24-early-trip-extension.json.gz /tmp/rae-profile 6 source.checkpoint
// No browser rendering, worker message transport, cold grid generation or file
// parsing is included. Use benchmark-grid-smoothing.js for cold grid preparation.
import fs from 'node:fs';
import path from 'node:path';
import { Session } from 'node:inspector';
import { createHash } from 'node:crypto';
import { gunzipSync } from 'node:zlib';
import { solveCoupledStreamtubeIses } from '../../src/euler/streamtube-coupled-ises.js';
import { summarizeCpuProfile } from './summarize-cpu-profile.js';

const [source, output, count = '3', member] = process.argv.slice(2), updates = Number(count);
if (!source || !output || !Number.isSafeInteger(updates) || updates < 1)
  throw new Error('Supply checkpoint, unused output directory, positive update count and optional dotted member path.');
const bytes = fs.readFileSync(source);
let checkpoint = JSON.parse(source.endsWith('.gz') ? gunzipSync(bytes) : bytes);
if (member) for (const key of member.split('.')) checkpoint = checkpoint[key];
if (checkpoint?.version !== 1 || !checkpoint.restart || !checkpoint.continuation)
  throw new Error('Select a complete coupled checkpoint.');
fs.mkdirSync(output); // Preserve earlier evidence instead of overwriting it.
const session = new Session(); session.connect();
const post = (method, params = {}) => {
  let result, error, completed = false;
  session.post(method, params, (e, r) => { error = e; result = r; completed = true; });
  if (!completed) throw new Error('Expected synchronous in-process inspector command.');
  if (error) throw error;
  return result;
};
post('Profiler.enable'); post('Profiler.setSamplingInterval', { interval: 1000 });
const hash = value => createHash('sha256').update(value).digest('hex');
const sourceHashes = Object.fromEntries(fs.readdirSync('src', { recursive: true })
  .map(p => path.join('src', p)).filter(p => fs.statSync(p).isFile()).sort()
  .map(p => [p, hash(fs.readFileSync(p))]));
post('Profiler.start');
const iterations = [], started = performance.now(); let last = started, result, milliseconds;
try {
  result = solveCoupledStreamtubeIses(undefined, { ...checkpoint.continuation, resume: checkpoint,
    maxIterations: updates, onIteration: h => {
      const now = performance.now();
      iterations.push({ iteration: h.iteration, milliseconds: now - last, residual: h.residual,
        step: h.step, backtracks: h.backtracks }); last = now;
    } });
} finally {
  milliseconds = performance.now() - started;
  const { profile } = post('Profiler.stop'); session.disconnect();
  fs.writeFileSync(path.join(output, 'cpu.cpuprofile'), JSON.stringify(profile));
  fs.writeFileSync(path.join(output, 'cpu-summary.json'), JSON.stringify(summarizeCpuProfile(profile), null, 2));
}
const resultHash = createHash('sha256');
for (const values of [result.x, result.residual]) resultHash.update(Buffer.from(Float64Array.from(values).buffer));
resultHash.update(JSON.stringify(result.flow.nodes));
const report = { source, member, checkpointFileHash: hash(bytes), sourceHashes, runtime: process.version,
  scope: 'Checkpoint replay and coupled updates; excludes imports, input parsing, cold initialization, browser rendering and worker messaging.',
  requestedUpdates: updates, milliseconds, iterations,
  unknowns: result.x.length, converged: result.converged, reason: result.reason, families: result.families,
  quality: result.mesh.quality, resultHash: resultHash.digest('hex') };
fs.writeFileSync(path.join(output, 'result.json'), JSON.stringify(report, null, 2));
console.log(JSON.stringify({ milliseconds: report.milliseconds, iterations, resultHash: report.resultHash }));
