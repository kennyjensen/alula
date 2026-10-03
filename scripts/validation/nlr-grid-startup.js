// SPDX-License-Identifier: GPL-2.0-or-later
// node scripts/validation/nlr-grid-startup.js [receipt.json] [alpha]
// Cold production-worker replay of the reported automatic-transition case.
import fs from 'node:fs';
import { getBenchmarkAirfoil } from '../../src/geometry/benchmark-airfoils.js';

const benchmark = getBenchmarkAirfoil('nlr7301');
const input = {
  elements: benchmark.elements, alpha: Number(process.argv[3] ?? 6),
  referenceChord: 1, momentReference: { x: .25, y: 0 },
  flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .185,
  reynolds: 2510000, ncrit: 9, transitionMode: 'automatic', materialTrips: [[1, 1], [1, 1]],
  gridIntervals: 32, gridTubes: 24, gridCrosslinePlacement: 'potential', eulerIsmom: 4,
  gridChordExponent: 0, gridSurfaceSpacing: 'automatic', gridEllipticSmoothing: true,
  gridSmoothingMethod: 'elliptic', geometrySource: benchmark.provenance,
};
const started = performance.now(), history = [];
let terminal;
globalThis.self = { postMessage(message) {
  if (message.type === 'iteration') {
    const h = message.iteration;
    const row = Object.fromEntries(['stage', 'coarseStage', 'iteration', 'residual', 'euler',
      'boundaryLayer', 'edgeMatching', 'step', 'backtracks'].filter(k => h[k] !== undefined).map(k => [k, h[k]]));
    history.push(row); console.log(JSON.stringify(row));
  }
  if (message.type === 'result' || message.type === 'error') {
    const r = message.result;
    terminal = { input, seconds: (performance.now() - started) / 1000, history,
      result: r ? { status: r.status, cl: r.cl, cd: r.cd, cm: r.cm,
        converged:r.converged, ncrit:r.restart?.options.ncrit, materialTrips:r.materialTrips,
        gridSequence:r.gridSequence, diagnostics: r.diagnostics,
        mesh: { cellCount: r.mesh?.cells?.length, quality: r.mesh?.quality },
        initialization: r.initialization } : message };
  }
} };
await import('../../src/worker/solver.js');
await self.onmessage({ data: { id: 1, task: 'solve', caseData: input } });
if (process.argv[2]) fs.writeFileSync(process.argv[2], JSON.stringify(terminal, null, 2));
console.log(JSON.stringify({ seconds: terminal?.seconds, status: terminal?.result.status,
  residual: terminal?.result.diagnostics?.equationResidual }));
if (!terminal?.result.converged || terminal.result.gridSequence?.reachedTarget !== true
  || terminal.result.ncrit !== input.ncrit
  || ![terminal.result.cl, terminal.result.cd, terminal.result.cm].every(Number.isFinite)) process.exitCode = 1;
