// Controlled same-process JS/WASM comparison, independent of panel tracing.
// Usage: node scripts/validation/benchmark-smoothing-kernel.js [stations=128] [tubes=24] [sweeps=20]
import assert from 'node:assert/strict';
import { createEllipticStreamtubeGrid } from '../../src/geometry/elliptic-streamtube-grid.js';
import { smoothPairedBoundaryGrid } from '../../src/geometry/paired-boundary-slor.js';
const nx = Number(process.argv[2] ?? 128), nt = Number(process.argv[3] ?? 24), maxSweeps = Number(process.argv[4] ?? 20);
const xi = Array.from({ length: nx + 1 }, (_, i) => i / nx), eta = Array.from({ length: nt + 1 }, (_, j) => j / nt);
const nodes = xi.map(u => eta.map(e => ({ x: u + .04 * Math.sin(Math.PI * u) * e * (1 - e), y: e + .025 * Math.sin(Math.PI * u) * Math.sin(Math.PI * e) })));
const config = { nodes, massFlows: eta.slice(1).map((v, j) => v - eta[j]), streamwiseCoordinates: xi,
  discretization: 'giles-1985', lineLinearization: 'full-metrics', lineGrouping: 'boundary-pairs', lineSearch: 'armijo',
  orthogonalBoundaryControl: { sourceForm: 'metric-stretch', background: nodes.map(row => row.map(() => 0)), decay: { lower: 4, upper: 4 } } };
let reference, maxCoordinateDifference = 0;
const samples = { javascript: [], wasm: [] };
for (let round = 0; round < 4; round++) for (const backend of round % 2 ? ['wasm', 'javascript'] : ['javascript', 'wasm']) {
  const system = createEllipticStreamtubeGrid(config), start = performance.now();
  const result = smoothPairedBoundaryGrid(system, { backend, maxSweeps });
  const elapsed = performance.now() - start;
  assert.equal(result.backend, backend); assert.equal(result.referenceReplays, 0);
  assert.ok(result.quality.valid); assert.equal(result.reason, 'sweep limit');
  if (!reference) reference = result;
  result.nodes.forEach((row, i) => row.forEach((p, j) => {
    for (const key of ['x', 'y']) maxCoordinateDifference = Math.max(maxCoordinateDifference, Math.abs(p[key] - reference.nodes[i][j][key]));
  }));
  assert.ok(maxCoordinateDifference < 1e-11);
  assert.deepEqual(result.history.map(row => row.rowSteps?.map(step => step.fraction)), reference.history.map(row => row.rowSteps?.map(step => step.fraction)));
  if (round) samples[backend].push(elapsed);
}
const median = values => [...values].sort((a, b) => a - b)[1];
console.log(JSON.stringify({ nx, nt, maxSweeps, milliseconds: samples, speedup: median(samples.javascript) / median(samples.wasm), maxCoordinateDifference }, null, 2));
