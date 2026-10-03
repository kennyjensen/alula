// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import path from 'node:path';
import { createHash } from 'node:crypto';
import { solverReliabilityCases, buildReliabilityCase } from './validation/solver-reliability-cases.js';
import { assessSolverResult } from './validation/solver-reliability-acceptance.js';
import { runBoundedWorker } from './validation/bounded-worker.js';
const args = Object.fromEntries(process.argv.slice(2).map(arg => {
  const match = /^--([a-z-]+)=(.+)$/.exec(arg);
  if (!match) throw new Error(`Use --key=value arguments; received ${arg}`);
  return [match[1], match[2]];
}));
for (const key of Object.keys(args)) if (!['stage', 'tier', 'case', 'seconds', 'out'].includes(key)) throw new Error(`Unknown option: ${key}`);
const stage = args.stage ?? 'geometry';
if (!['geometry', 'mesh', 'solve'].includes(stage)) throw new Error('Stage must be geometry, mesh or solve.');
const timeoutMs = 1000 * Number(args.seconds ?? (stage === 'geometry' ? 15 : 180));
if (!(Number.isFinite(timeoutMs) && timeoutMs > 0 && timeoutMs <= 3600000)) throw new Error('Time limit must be in (0, 3600] seconds per case.');
const selectedIds = args.case?.split(',');
if (selectedIds?.some(id => !solverReliabilityCases.some(c => c.id === id))) throw new Error('Unknown case ID; see scripts/validation/solver-reliability-cases.js.');
const cases = solverReliabilityCases.filter(c => (!selectedIds || selectedIds.includes(c.id)) && (!args.tier || args.tier === 'all' || c.tier === args.tier));
if (!cases.length) throw new Error('No matching cases.');
// Full solves are explicitly requested. Default runs only cheap geometry
// gates, never an unbounded Cartesian sweep of all visible settings.
const out = args.out ?? `docs/solver-reliability/${new Date().toISOString().replace(/[:.]/g, '-')}-${stage}`;
if (fs.existsSync(out)) throw new Error(`Preserve existing report directory: ${out}`);
fs.mkdirSync(out, { recursive: true });
const hash = value => createHash('sha256').update(value).digest('hex');
const serialize = value => JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v, 2) + '\n';
const write = (name, value) => fs.writeFileSync(path.join(out, name), serialize(value));
function hashes() {
  const paths = fs.readdirSync('src', { recursive: true }).map(p => path.join('src', p)).filter(p => fs.statSync(p).isFile());
  paths.push('index.html', 'scripts/validate-solver-reliability.js', ...fs.readdirSync('scripts/validation').map(p => path.join('scripts/validation', p)).filter(p => fs.statSync(p).isFile()));
  return Object.fromEntries(paths.sort().map(p => [p, hash(fs.readFileSync(p))]));
}
const sourceHashes = hashes(); write('sources.json', sourceHashes);
const report = { startedAt: new Date().toISOString(), stage, timeoutSeconds: timeoutMs / 1000,
  scope: 'Current source, actual browser Worker message handler in isolated Node threads; browser UI and panel streamlines emitted after the flow result are tested separately. Geometry passes are not flow convergence or physical validation.',
  physicalValidation: 'not evaluated', selected: cases.map(c => c.id), notRun: solverReliabilityCases.filter(c => !cases.includes(c)).map(c => c.id), cases: [], inProgress: true };
write('report.json', report);
for (const spec of cases) {
  const entry = { id: spec.id, preset: spec.preset, mode: spec.mode, tier: spec.tier, passed: false };
  console.log(JSON.stringify({ event: 'start', case: spec.id, stage }));
  try {
    const { caseData, guiRestriction } = buildReliabilityCase(spec);
    entry.guiRestriction = guiRestriction || null;
    write(`${spec.id}.input.json`, { spec, caseData, guiRestriction });
    entry.inputHash = hash(serialize(caseData));
    if (stage === 'mesh' && caseData.flowModel !== 'streamtube-grid') {
      entry.status = 'not-applicable'; entry.reason = 'Panel modes have no volume mesh.';
    } else {
      const outcome = await runBoundedWorker(new URL('./validation/bounded-solver-worker.js', import.meta.url), { caseData, stage }, {
        timeoutMs, onProgress: message => {
          // Persist compact progress even if the orchestrator itself stops.
          fs.appendFileSync(path.join(out, `${spec.id}.progress.jsonl`), JSON.stringify(message) + '\n');
        },
      });
      write(`${spec.id}.state.json`, outcome);
      entry.seconds = outcome.seconds; entry.status = outcome.type;
      if (stage === 'solve' && outcome.type === 'result') {
        entry.acceptance = assessSolverResult(caseData, outcome.result); entry.passed = entry.acceptance.passed;
      } else if (stage === 'geometry') {
        entry.passed = outcome.type === 'geometry-ready'; entry.geometry = outcome.geometry;
      } else if (stage === 'mesh') entry.passed = outcome.type === 'mesh-ready' && outcome.retained.mesh?.quality?.valid === true;
      if (!entry.passed) entry.reason = outcome.message ?? outcome.result?.diagnostics?.reason ?? outcome.result?.reason ?? entry.acceptance?.failures?.join(', ');
      entry.lastStage = outcome.stage ?? outcome.retained.stage ?? stage;
    }
  } catch (error) { entry.status = 'error'; entry.reason = error.message; entry.code = error.code; }
  report.cases.push(entry); write('report.json', report);
  console.log(JSON.stringify({ event: 'terminal', ...entry }));
}
const finalHashes = hashes();
report.sourceChanges = [...new Set([...Object.keys(sourceHashes), ...Object.keys(finalHashes)])].filter(p => sourceHashes[p] !== finalHashes[p]);
report.inProgress = false; report.finishedAt = new Date().toISOString();
report.passed = report.cases.every(c => c.passed || c.status === 'not-applicable') && !report.sourceChanges.length;
report.publicModeNumericalAcceptance = stage === 'solve' && report.passed && report.cases.every(c => !c.guiRestriction);
write('report.json', report);
console.log(JSON.stringify({ report: path.join(out, 'report.json'), passed: report.passed, publicModeNumericalAcceptance: report.publicModeNumericalAcceptance, sourceChanges: report.sourceChanges }));
if (!report.passed) process.exitCode = 1;
