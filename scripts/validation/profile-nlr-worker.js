// SPDX-License-Identifier: GPL-2.0-or-later
// node scripts/validation/profile-nlr-worker.js output-directory [source-root]
// Runs the real Worker handler in Node, including progress observables and
// structured-clone cost. Browser rendering/IPC are outside this CPU profile.
import fs from 'node:fs';
import path from 'node:path';
import { Session } from 'node:inspector';
import { fileURLToPath, pathToFileURL } from 'node:url';

const repository = fileURLToPath(new URL('../../', import.meta.url));
const output = process.argv[2], root = path.resolve(process.argv[3] ?? repository);
if (!output) throw new Error('Supply an output directory for the CPU profiles and result.');
fs.mkdirSync(output, { recursive: true });
const receipt = JSON.parse(fs.readFileSync(path.join(repository, 'docs/solver-reliability/nlr64x24-subsonic/browser.json')));
const session = new Session(); session.connect();
// The in-process inspector dispatches these commands synchronously. Keeping
// rotation synchronous prevents timing boundaries from drifting during solve.
const post = (method, params = {}) => {
  let result, error, completed = false;
  session.post(method, params, (e, r) => { error = e; result = r; completed = true; });
  if (!completed) throw new Error('Profiler requires synchronous in-process inspector commands.');
  if (error) throw error;
  return result;
};
post('Profiler.enable'); post('Profiler.setSamplingInterval', { interval: 2000 });
let segment = 0, segmentStart = performance.now(), started = segmentStart, outputMilliseconds = 0;
const history = [], messages = {}, profiles = [];
post('Profiler.start');
function rotate(label) {
  const now = performance.now(), { profile } = post('Profiler.stop');
  const name = `${String(segment++).padStart(2, '0')}-${label}.cpuprofile`;
  fs.writeFileSync(path.join(output, name), JSON.stringify(profile));
  const entry = { profile: name, milliseconds: now - segmentStart, elapsed: now - started };
  profiles.push(entry); console.log(JSON.stringify(entry));
  segmentStart = performance.now(); post('Profiler.start');
}
let terminal;
globalThis.self = { postMessage(message) {
  const clock = performance.now(); structuredClone(message); outputMilliseconds += performance.now() - clock;
  messages[message.type] = (messages[message.type] ?? 0) + 1;
  if (message.type === 'flow-stage') rotate(message.stage);
  if (message.type === 'iteration') {
    const h = message.iteration; history.push({ ...h, elapsed: performance.now() - started });
    console.log(JSON.stringify({ stage: h.stage, iteration: h.iteration, residual: h.residual, elapsed: performance.now() - started }));
    if (h.stage === 'coupled' && h.iteration > 0 && h.iteration % 6 === 0) rotate(`coupled-${h.iteration}`);
  }
  if (message.type === 'result' || message.type === 'error') {
    rotate('terminal'); const r = message.result;
    terminal = { input: receipt.input, sourceRoot: root, runtime: process.version, route: 'node-worker-handler',
      elapsed: performance.now() - started, outputMilliseconds, messages, profiles, history,
      result: r ? { status: r.status, diagnostics: r.diagnostics, cl: r.cl, cd: r.cd, cm: r.cm,
        initialization: r.initialization } : message };
    fs.writeFileSync(path.join(output, 'result.json'), JSON.stringify(terminal));
    console.log(JSON.stringify({ terminal: true, status: r?.status, residual: r?.diagnostics?.equationResidual,
      elapsed: terminal.elapsed, outputMilliseconds }));
  }
} };
try {
  await import(pathToFileURL(path.join(root, 'src/worker/solver.js')).href);
  await self.onmessage({ data: { id: 1, task: 'solve', caseData: receipt.input } });
  if (terminal?.result.status !== 'research-coupled-equations-converged') process.exitCode = 1;
} finally { post('Profiler.stop'); session.disconnect(); }
