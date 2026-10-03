// Full app mesh preparation, including panel tracing, SLOR and wall refinement.
// Usage: node scripts/validation/benchmark-grid-smoothing.js [intervals=64] [tubes=11]
import { createHash } from 'node:crypto';
import { buildReliabilityCase } from './solver-reliability-cases.js';
import { prepareStreamtubeAssembly } from '../../src/euler/streamtube-result.js';

const gridIntervals = Number(process.argv[2] ?? 64), gridTubes = Number(process.argv[3] ?? 11);
const { caseData } = buildReliabilityCase({ preset: 'rae2822-mses', mode: 'streamtube-grid',
  changes: { mach: .74, alpha: 2.68, gridIntervals, gridTubes } });
const start = performance.now(), result = prepareStreamtubeAssembly(caseData, { meshOnly: true });
const milliseconds = performance.now() - start, smoothing = result.diagnostics.gridSmoothing;
console.log(JSON.stringify({ gridIntervals, gridTubes, milliseconds, gridValid: result.mesh.quality.valid,
  converged: smoothing.converged,
  regions: smoothing.regions.map(region => ({ reason: region.reason, historyEntries: region.history.length,
    residual: region.history.at(-1)?.residual, quality: region.quality })),
  nodesHash: createHash('sha256').update(JSON.stringify(result.nodes)).digest('hex'),
  smoothingHash: createHash('sha256').update(JSON.stringify(smoothing)).digest('hex'),
}, null, 2));
if (!result.mesh.quality.valid || !smoothing.converged) process.exitCode = 1;
