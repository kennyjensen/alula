import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { observableFlow, observableBL } from '../src/euler/streamtube-flow-preview.js';

test('Euler precursor and standalone Euler forward accepted loads and recover from unavailable loads', async () => {
  const messages = [], previous = globalThis.self;
  const coefficients = { cl: .5, cd: .02, cm: -.1 };
  const solve = (_input, options) => {
    options.onMesh({}, 'initial'); // Geometry alone cannot supply loads.
    for (let iteration = 0; iteration < 3; iteration++) {
      options.onIteration({ iteration, residual: 1 / (iteration + 1) });
      options.onMesh({ iteration: { iteration, stage: 'euler', startupAttempt: 1 },
        coefficientProgress: { kind: 'euler-pressure', iteration, mach: .47, actualAlpha: 2,
          gridLevel: 16, ...(iteration === 1 ? { message: 'Missing wall pressure' } : { coefficients }) } }, 'solving');
    }
    return {};
  };
  globalThis.__eulerCoefficientTest = { solveCoupledStreamtubeAssembly: solve,
    solveStreamtubeAssembly: solve, quadCoupledResultForDisplay: x => x };
  globalThis.self = { postMessage: x => messages.push(structuredClone(x)) };
  try {
    const location = new URL('../src/worker/solver.js', import.meta.url);
    let source = fs.readFileSync(location, 'utf8').replace(/await import\('\.\.\/(?:euler\/streamtube-(?:coupled-assembly|result)|ui\/quad-coupled-result)\.js'\)/g,
      'globalThis.__eulerCoefficientTest');
    source = source.replace(/from '(\.[^']+)'/g, (_m, relative) => `from '${new URL(relative, location).href}'`)
      .replace(/import\('(\.[^']+)'\)/g, (_m, relative) => `import('${new URL(relative, location).href}')`);
    await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
    for (const quadBoundaryLayers of [true, false]) {
      messages.length = 0;
      await self.onmessage({ data: { id: 17, caseData: { flowModel: 'streamtube-grid', quadBoundaryLayers,
        mach: .74, alpha: 2.68, elements: [], referenceChord: 1 } } });
      const loads = messages.filter(m => m.type.startsWith('coefficients'));
      assert.deepEqual(loads.map(m => [m.type, m.iteration]),
        [['coefficients', 0], ['coefficients-unavailable', 1], ['coefficients', 2]]);
      for (const message of loads) {
        assert.equal(message.id, 17);
        assert.equal(message.kind, 'euler-pressure');
        assert.equal(message.mach, .47);
        assert.equal(message.targetMach, .74);
        assert.equal(message.actualAlpha, 2);
        assert.equal(message.targetAlpha, 2.68);
        assert.equal(message.gridLevel, 16);
        assert.equal(message.stage, 'euler');
        if (message.type === 'coefficients') assert.deepEqual(message.coefficients, coefficients);
      }
      assert.equal(messages.at(-1).type, 'result');
      assert.equal(messages.some(m => m.type === 'error'), false);
    }
  } finally { globalThis.self = previous; delete globalThis.__eulerCoefficientTest; }
});

// Exercise the real Worker reporting code and numerical coefficient routines.
// Substitute the flow solve with two already captured physical snapshots.
test('ordinary startup and cold Mach startup publish provisional loads as well as Cp', async () => {
  const saved = JSON.parse(fs.readFileSync(new URL('../docs/current-multielement-automatic-16x9-slor-browser.json', import.meta.url)));
  const r = saved.result, cp = saved.checkpoint, messages = [], previous = globalThis.self;
  const frame = { checkpoint: cp, flow: observableFlow(r.flow),
    bl: observableBL({ ...r.numericalBoundaryLayer, scale: 1 / Math.sqrt(r.kernelReynolds) }),
    bodies: cp.restart.input.bodies, normalization: { solverLength: r.solverLength, referenceChord: r.referenceChord },
    mach: cp.restart.input.mach };
  const solve = (_input, options) => {
    for (const [iteration, stage] of [[0, 'coupled'], [1, 'coupled-ncrit-startup']])
      options.onFlow({ ...frame, stage, iteration: { iteration }, startupAttempt: 1 });
    return {};
  };
  globalThis.__liveCoefficientTest = { solveCoupledStreamtubeAssembly: solve, solveCoupledStreamtubeMach: solve,
    quadCoupledResultForDisplay: x => x };
  globalThis.self = { postMessage: x => messages.push(structuredClone(x)) };
  try {
    const location = new URL('../src/worker/solver.js', import.meta.url);
    let source = fs.readFileSync(location, 'utf8').replace(/await import\('\.\.\/(?:euler\/streamtube-coupled-(?:assembly|mach-assembly)|ui\/quad-coupled-result)\.js'\)/g,
      'globalThis.__liveCoefficientTest');
    source = source.replace(/from '(\.[^']+)'/g, (_m, relative) => `from '${new URL(relative, location).href}'`)
      .replace(/import\('(\.[^']+)'\)/g, (_m, relative) => `import('${new URL(relative, location).href}')`);
    await import('data:text/javascript;base64,' + Buffer.from(source).toString('base64'));
    for (const mach of [.2, .74]) {
      messages.length = 0;
      await self.onmessage({ data: { id: 1, caseData: { flowModel: 'streamtube-grid', quadBoundaryLayers: true,
        mach, ncrit: cp.restart.options.ncrit, referenceChord: r.referenceChord, momentReference: r.momentReference,
        elements: cp.restart.input.bodies.map((b, i) => ({ name: `Element ${i}` })) } } });
      const loads = messages.filter(m => m.type === 'coefficients');
      assert.equal(loads.length, 2, JSON.stringify(messages.map(m => ({ type: m.type, message: m.message }))));
      assert.equal(messages.some(m => m.type === 'coefficients-unavailable'), false);
      for (const message of loads) {
        assert.deepEqual(message.coefficients, r.coefficients);
        assert(message.pressure.elements.length > 0);
        assert.equal(message.mach, cp.restart.input.mach);
        assert.equal(message.targetMach, mach);
      }
    }
  } finally { globalThis.self = previous; delete globalThis.__liveCoefficientTest; }
});
