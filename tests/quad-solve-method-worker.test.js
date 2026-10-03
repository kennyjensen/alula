import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadCoupledNcrit, quadCoupledNcritResult } from '../src/ui/quad-coupled-ncrit.js';

for (const route of ['cold', 'mach', 'refinement']) test(`${route} Worker progress labels use the current physical checkpoint`, async () => {
  const key = `__solveMethodWorker_${route}`, messages = [];
  const request = { flowModel: 'streamtube-grid', quadBoundaryLayers: true, mach: .2, ncrit: 9,
    transitionMode: 'automatic', referenceChord: 1, elements: [] };
  const checkpoint = { version: 1, restart: { input: { mach: .2 }, options: { ncrit: 4, hkFloorLinearization: 'native' } },
    continuation: { shearCoordinate: 'logarithmic' } };
  const raw = { flow: {}, boundaryLayer: {}, checkpoint, conditions: { mach: .2, ncrit: 4 },
    mesh: { quality: { valid: true } }, sourceCase: request };
  const frame = { checkpoint, stage: route === 'mach' ? 'coupled-mach' : 'coupled', iteration: { iteration: 41 },
    shearCoordinate: 'linear', hkFloorLinearization: 'exact' };
  const solve = (input, options) => { options.onFlow(frame); return raw; };
  const modules = {
    '../euler/streamtube-coupled-assembly.js': { solveCoupledStreamtubeAssembly: solve },
    '../euler/streamtube-coupled-mach-assembly.js': { solveCoupledStreamtubeMach: solve },
    '../euler/streamtube-coupled-refinement-assembly.js': { solveCoupledStreamtubeRefinement(input, parent, options) {
      options.onPrepared({ checkpoint, system: { bl: {}, euler: { layout: { bodies: [] }, conditions: { gamma: 1.4 } } },
        settings: { normalization: { solverLength: 1, referenceChord: 1 } } });
      options.onMesh(raw.mesh, 'solving', { flow: {}, iteration: { iteration: 41 } });
      return raw;
    } },
    '../ui/quad-coupled-result.js': { quadCoupledResultForDisplay: result => result },
    '../ui/quad-mesh-progress.js': { createQuadMeshProgress: () => ({ stage() {}, mesh: value => value }) },
    '../ui/quad-coupled-transonic-coefficients.js': { quadCoupledTransonicCoefficients: () => ({}) },
    '../ui/quad-coupled-iteration-coefficients.js': { quadCoupledIterationCoefficients: () => ({}) },
  };
  const self = { postMessage: message => messages.push(structuredClone(message)) };
  globalThis[key] = { self, modules, quadCoupledNcrit, quadCoupledNcritResult, quadCoupledIterationPressure: () => ({ elements: [] }) };
  const text = fs.readFileSync(new URL('../src/worker/solver.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`)
    .replace(/await import\('([^']+)'\)/g, (_, path) => `globalThis[${JSON.stringify(key)}].modules[${JSON.stringify(path)}]`);
  try {
    await import('data:text/javascript;base64,' + Buffer.from(`const self=globalThis.${key}.self;\n${text}`).toString('base64'));
    await self.onmessage({ data: { id: 7, caseData: request,
      ...(route === 'mach' ? { task: 'continue-coupled' } : route === 'refinement' ? { task: 'refine-coupled', parentResult: raw } : {}) } });
    assert.equal(messages.some(message => message.type === 'error'), false);
    const progress = messages.find(message => ['pressure', 'coefficients'].includes(message.type));
    assert.ok(progress);
    assert.equal(progress.shearCoordinate, 'logarithmic');
    assert.equal(progress.hkFloorLinearization, 'native');
    assert.equal(progress.actualNcrit, 4);
    assert.equal(progress.iteration, 41);
    assert.equal(frame.shearCoordinate, 'linear', 'Presentation must not modify numerical observations.');
  } finally { delete globalThis[key]; }
});
