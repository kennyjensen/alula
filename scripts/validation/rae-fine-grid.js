// Bounded cold-start comparisons at fixed physical conditions. All grid
// overrides are public input controls. Receipts are experiments, not app status.
// Usage: node scripts/validation/rae-fine-grid.js /tmp/rae128x7 '{"gridTubes":7}'
import fs from 'node:fs';
import { buildReliabilityCase } from './solver-reliability-cases.js';
import { raeSourceHash } from './rae-status.js';
import { solveStreamtubeAssembly } from '../../src/euler/streamtube-result.js';
import { createStreamtubeBodySystem } from '../../src/euler/streamtube-body.js';
import { capturePreparedStreamtubeAssembly } from '../../src/euler/streamtube-prepared-assembly.js';
const prefix = process.argv[2];
if (!prefix) throw new Error('Supply an output prefix.');
const changes = JSON.parse(process.argv[3] ?? '{}');
const { caseData: input } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-grid',
  changes: { mach: .74, alpha: 2.68, gridIntervals: 128, gridTubes: 7, eulerStartup: 'standard', ...changes } });
const receipt = { recordedAt: new Date().toISOString(), sourceHash: raeSourceHash(), route: 'experiment', input, history: [] };
receipt.preparedSource = process.env.RAE_PREPARED;
receipt.experimentalMucon = process.env.RAE_TEMP_MUCON;
let preparedEuler;
if (process.env.RAE_PREPARED) {
  const saved = JSON.parse(fs.readFileSync(process.env.RAE_PREPARED));
  const p = saved.prepared, system = createStreamtubeBodySystem(p.input);
  const initial = system.adoptGeometry(Float64Array.from(Object.values(p.initial)), p.nodes);
  preparedEuler = capturePreparedStreamtubeAssembly({ ...p, system, initial, initialEuler: { x: initial, nodes: p.nodes } }, saved.caseData);
}
const start = performance.now();
try {
  const result = solveStreamtubeAssembly(input, {
    maxIterations: Number(process.env.RAE_MAX_ITERATIONS ?? 80),
    ...(preparedEuler ? { preparedEuler } : {}),
    onEulerPrepared: p => fs.writeFileSync(`${prefix}-prepared.json`, JSON.stringify(p)),
    onIteration: h => {
      receipt.history.push(h);
      console.log(JSON.stringify({ iteration: h.iteration, residual: h.residual, step: h.step,
        backtracks: h.backtracks, mucon: h.dissipation?.mucon, mcrit: h.dissipation?.mcrit }));
    },
  });
  receipt.result = { status: result.status, diagnostics: result.diagnostics,
    settings: result.solverSettings, quality: result.mesh.quality, initialization: result.mesh.initialization,
    cells: result.mesh.cells.length, law: result.flow.solverInput.upwind,
    lastRejectedStep: result.flow.lastRejectedStep };
  fs.writeFileSync(`${prefix}-checkpoint.json`, JSON.stringify(result.flow.checkpoint));
} catch (e) { receipt.error = { message: e.message, code: e.code, diagnostics: e.diagnostics }; }
receipt.seconds = (performance.now() - start) / 1000;
fs.writeFileSync(`${prefix}.json`, JSON.stringify(receipt, null, 2) + '\n');
console.log(JSON.stringify({ status: receipt.result?.status, residual: receipt.result?.diagnostics?.equationResidual,
  reason: receipt.result?.diagnostics?.reason, error: receipt.error, seconds: receipt.seconds }));
