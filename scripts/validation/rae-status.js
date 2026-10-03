// The table is a view of cold browser-form receipts, never experiment labels.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
import { isDeepStrictEqual } from 'node:util';
import { buildReliabilityCase } from './solver-reliability-cases.js';

const root = fileURLToPath(new URL('../../', import.meta.url));
const reportPath = path.join(root, 'docs/rae-status.md');
const evidencePath = path.join(root, 'docs/rae-status-results.json');
export const raeGrids = [8, 16, 32, 64, 128].flatMap(n => (n === 128 ? [7, 9, 11, 24] : [7, 9, 11]).map(t => [n, t]));
export function raeSourceHash() {
  const files = ['index.html', 'scripts/validation/solver-reliability-cases.js',
    'scripts/validation/rae-status.js', 'tests/browser/rae-default-startup.spec.js'];
  function visit(dir) {
    for (const entry of fs.readdirSync(path.join(root, dir), { withFileTypes: true })) {
      const name = `${dir}/${entry.name}`;
      if (entry.isDirectory()) { if (entry.name !== 'tests') visit(name); }
      else files.push(name);
    }
  }
  visit('src');
  const hash = createHash('sha256');
  for (const name of files.sort()) hash.update(name).update('\0').update(fs.readFileSync(path.join(root, name))).update('\0');
  return hash.digest('hex');
}
export function verifiedRaeReceipt(receipt, sourceHash) {
  if (!receipt || receipt.sourceHash !== sourceHash || receipt.route !== 'browser-form'
    || !receipt.browser || !receipt.recordedAt) return false;
  const input = receipt.input, r = receipt.result;
  if (!input || !raeGrids.some(([n, t]) => n === input.gridIntervals && t === input.gridTubes)) return false;
  const viscous = input.quadBoundaryLayers === true;
  const { caseData: expected } = buildReliabilityCase({ preset: 'rae2822-mses', mode: viscous ? 'streamtube-bl' : 'streamtube-grid',
    changes: { mach: .74, alpha: 2.68, gridIntervals: input.gridIntervals, gridTubes: input.gridTubes,
      ...(viscous ? { reynolds: 2.7e6, ncrit: 4 } : { eulerStartup: 'standard' }) } });
  // Reject hidden overrides, different geometry, startup modes and tolerances.
  if (!isDeepStrictEqual(input, expected)) return false;
  return !receipt.error && r?.status === (viscous ? 'research-coupled-equations-converged' : 'research-converged')
    && r.mach === .74 && r.alpha === 2.68 && r.gridValid === true
    && Number.isFinite(r.residual) && r.residual <= 1e-10
    && r.upwind?.mucon === 1 && r.upwind?.mcrit === .99;
}
export function renderRaeStatus(receipts, sourceHash) {
  const status = (n, t, viscous) => receipts.some(r => r.input?.gridIntervals === n && r.input?.gridTubes === t
    && r.input?.quadBoundaryLayers === viscous && verifiedRaeReceipt(r, sourceHash)) ? 'works in app' : 'unresolved';
  return `Cold browser-form solves: MIT RAE 2822, Mach 0.74, α 2.68°, Standard (automatic) startup, automatic spacing, smoothing on, inlet/outlet auto. Viscous: Re 2.7e6, Ncrit 4, automatic transition.\n\nOnly a passing browser receipt for the current source and these exact settings earns “works in app”. “Unresolved” includes cases not reverified after code changes; historical experiments are not app verification. Regenerate with \`npm run status:rae\`; the tests reject a stale table. Receipts: [rae-status-results.json](rae-status-results.json). Verification applies to the recorded build/browser, not an older deployed or cached app.\n\n| Grid | Inviscid solve | Viscous solve |\n| --- | --- | --- |\n`
    + raeGrids.map(([n, t]) => `| rae${n}x${t} | ${status(n, t, false)} | ${status(n, t, true)} |`).join('\n') + '\n';
}
export function readRaeReceipts() {
  return fs.existsSync(evidencePath) ? JSON.parse(fs.readFileSync(evidencePath, 'utf8')) : [];
}
if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  fs.mkdirSync(path.dirname(reportPath), { recursive: true });
  const sourceHash = raeSourceHash();
  let receipts = readRaeReceipts();
  const record = process.argv.find(a => a.startsWith('--record='));
  if (record) {
    const receipt = JSON.parse(fs.readFileSync(record.slice('--record='.length), 'utf8'));
    if (!verifiedRaeReceipt(receipt, sourceHash)) throw new Error('Receipt is stale, failed, or differs from the documented app settings.');
    receipts = receipts.filter(r => r.input?.gridIntervals !== receipt.input.gridIntervals
      || r.input?.gridTubes !== receipt.input.gridTubes || r.input?.quadBoundaryLayers !== receipt.input.quadBoundaryLayers);
    // Detailed iteration logs remain in the test artifact, not the status file.
    const { progress, ...summary } = receipt;
    receipts.push(summary);
    fs.writeFileSync(evidencePath, JSON.stringify(receipts, null, 2) + '\n');
  }
  const report = renderRaeStatus(receipts, sourceHash);
  if (process.argv.includes('--check')) {
    if (fs.readFileSync(reportPath, 'utf8') !== report) throw new Error('RAE status is stale. Run npm run status:rae; rerun browser verification to restore success claims.');
  } else fs.writeFileSync(reportPath, report);
}
