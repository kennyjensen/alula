// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledNcrit, quadCoupledNcritResult } from '../src/ui/quad-coupled-ncrit.js';
import { quadCoupledResultForDisplay } from '../src/ui/quad-coupled-result.js';
import { quadCoupledIterationPressure } from '../src/ui/quad-coupled-iteration-pressure.js';
import { createQuadMeshProgress } from '../src/ui/quad-mesh-progress.js';

test('actual worker, pressure, mesh and terminal adapters keep a retained grid paired with its own BL/Cp', async () => {
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/current-multielement-automatic-16x9-slor-browser.json', import.meta.url)));
  const request = { ...saved.input, gridIntervals: 128 }, messages = [];
  // Existing physical fixture; only sequence metadata is synthetic. No flow solve.
  const raw = { ...structuredClone(saved.result), boundaryLayer: structuredClone(saved.result.numericalBoundaryLayer),
    sourceCase: structuredClone(request), reason: 'Controlled finer-grid failure', stateConverged: true,
    gridSequence: { kind: 'coarse-to-fine', reachedTarget: false, actualGridIntervals: 16, requestedGridIntervals: 128 } };
  raw.automaticRefinement = raw.gridSequence;
  // Intermediate grid-level results have transfer metadata but no Euler
  // precursor record (the shape which crashed the reported 128/11 solve).
  raw.initialization = { method: 'Complete-state nested grid sequencing' };
  const key = '__gridWorkerFixture', self = { postMessage: message => messages.push(structuredClone(message)) };
  const modules = {
    '../euler/streamtube-coupled-assembly.js': { solveCoupledStreamtubeAssembly(input, options) {
      assert.deepEqual(input, request);
      const labels = { stage: 'coupled-grid-refinement', gridLevel: 16, requestedGridIntervals: 128, retained: true,
        actualNcrit: 9, targetNcrit: 9 };
      options.onStage(labels);
      options.onMesh(raw.mesh, 'solving', labels);
      options.onFlow({ ...labels, checkpoint: raw.checkpoint, flow: raw.flow, bl: raw.boundaryLayer,
        bodies: raw.solverInput.bodies, referenceChord: raw.referenceChord,
        iteration: raw.history.at(-1), mach: raw.conditions.mach });
      return raw;
    } },
    '../ui/quad-coupled-result.js': { quadCoupledResultForDisplay },
    '../ui/quad-mesh-progress.js': { createQuadMeshProgress },
  };
  globalThis[key] = { self, modules, quadCoupledNcrit, quadCoupledNcritResult, quadCoupledIterationPressure,
    quadCoupledIterationCoefficients: () => structuredClone(saved.result.coefficients) };
  const text = fs.readFileSync(new URL('../src/worker/solver.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names.replace(/\s+as\s+/g, ':')}} = globalThis[${JSON.stringify(key)}];`)
    .replace(/await import\('([^']+)'\)/g, (_, path) => `globalThis[${JSON.stringify(key)}].modules[${JSON.stringify(path)}]`);
  try {
    await import('data:text/javascript;base64,' + Buffer.from(`const self=globalThis.${key}.self;\n${text}`).toString('base64'));
    await self.onmessage({ data: { id: 17, caseData: request } });
    assert.deepEqual(messages.map(m => m.type), ['flow-stage', 'mesh', 'coefficients', 'result']);
    const pressure = messages[2], result = messages[3].result;
    assert.equal(pressure.gridLevel, 16); assert.equal(pressure.requestedGridIntervals, 128); assert.equal(pressure.retained, true);
    assert.deepEqual(pressure.pressure.elements.map(e => e.cp), result.elements.map(e => e.cp));
    assert.deepEqual(result.mesh.vertices, raw.mesh.vertices);
    assert.equal(result.sourceCase.gridIntervals, 16); assert.equal(result.requestedCase.gridIntervals, 128);
    assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
    assert.equal(result.coefficientStatus, 'unconverged');
    assert.equal(raw.sourceCase.gridIntervals, 128); assert.equal(raw.converged, true, 'Adapters do not mutate the saved numerical state.');
  } finally { delete globalThis[key]; }
});
