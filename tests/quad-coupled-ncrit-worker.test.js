import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledNcrit, quadCoupledNcritResult } from '../src/ui/quad-coupled-ncrit.js';

test('actual worker pressure and terminal source-case messages preserve intermediate checkpoint Ncrit', async () => {
  const messages = [], request = { flowModel: 'streamtube-grid', quadBoundaryLayers: true, transitionMode: 'automatic',
    mach: .185, ncrit: 9, referenceChord: 1, elements: [{ name: 'fixture', points: [] }] };
  const source = { converged: true, conditions: { ncrit: 8, mach: .185 }, mesh: { quality: { valid: true } },
    checkpoint: { restart: { options: { ncrit: 8 } } }, sourceCase: { ...request },
    actualNcrit: 8, targetNcrit: 9, ncritContinuation: { actualNcrit: 8, targetNcrit: 9, reachedTarget: false } };
  const key = '__ncritWorkerFixture', self = { postMessage: message => messages.push(structuredClone(message)) };
  const modules = {
    '../euler/streamtube-coupled-assembly.js': { solveCoupledStreamtubeAssembly(input, options) {
      assert.deepEqual(input, request);
      options.onStage({ stage: 'coupled', actualNcrit: 7.5, targetNcrit: 9 });
      options.onFlow({ checkpoint: { version: 1, restart: { options: { ncrit: 7.5 }, input: { mach: .185 } } },
        actualNcrit: 9, targetNcrit: 9, iteration: { iteration: 2 }, stage: 'coupled' });
      return source;
    } },
    '../ui/quad-coupled-result.js': { quadCoupledResultForDisplay: result => result },
    '../ui/quad-mesh-progress.js': { createQuadMeshProgress: () => ({ stage() {}, mesh: value => value }) },
  };
  globalThis[key] = { self, modules, quadCoupledNcrit, quadCoupledNcritResult,
    quadCoupledIterationPressure: frame => ({ elements: [], actualMach: frame.checkpoint.restart.input.mach,
      fixtureActualNcrit: frame.checkpoint.restart.options.ncrit }) };
  const text = fs.readFileSync(new URL('../src/worker/solver.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`)
    .replace(/await import\('([^']+)'\)/g, (_, path) => `globalThis[${JSON.stringify(key)}].modules[${JSON.stringify(path)}]`);
  try {
    await import('data:text/javascript;base64,' + Buffer.from(`const self=globalThis.${key}.self;\n${text}`).toString('base64'));
    await self.onmessage({ data: { id: 5, caseData: request } });
    assert.deepEqual(messages.map(m => m.type), ['flow-stage', 'pressure', 'result']);
    const pressure = messages[1], result = messages[2].result;
    assert.equal(pressure.actualNcrit, 7.5); assert.equal(pressure.targetNcrit, 9);
    assert.equal(pressure.pressure.fixtureActualNcrit, 7.5);
    assert.equal(result.sourceCase.ncrit, 8); assert.equal(result.requestedCase.ncrit, 9);
    assert.equal(result.converged, false); assert.equal(result.stateConverged, true);
    assert.equal(source.sourceCase.ncrit, 9); assert.equal(source.converged, true, 'Worker labels must not mutate numerical result.');
  } finally { delete globalThis[key]; }
});
