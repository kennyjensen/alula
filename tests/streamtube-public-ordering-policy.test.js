import * as shearPolicy from '../src/euler/streamtube-coupled-shear-policy.js';
import { coupledCheckpointHkPolicy } from '../src/euler/streamtube-coupled-assembly.js';
// SPDX-License-Identifier: GPL-2.0-or-later
// Public adapter seams only: no numerical construction, Newton or factorization.
import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';

const encode = text => 'data:text/javascript;base64,' + Buffer.from(text).toString('base64');
let serial = 0;
async function stubbed(file, hooks) {
  const key = `__publicOrdering${serial++}`;
  globalThis[key] = { ...shearPolicy, coupledCheckpointHkPolicy, ...hooks };
  const source = fs.readFileSync(new URL(file, import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}}=globalThis[${JSON.stringify(key)}];`);
  try { return await import(encode(source)); } finally { delete globalThis[key]; }
}

test('transition refinement preserves selected ordering and resets only the previous grid fallback history', async () => {
  for (const policy of [undefined, 'auto', 'station', 'station-auto', 'aligned-auto']) {
    const families = { euler: .1, boundaryLayer: .2, edgeMatching: .3 }, calls = [];
    const input = { mach: .2, wakeGeometry: 'independent-banks', streamwiseMode: 'hybrid', hybrid: { ismom: 4 } };
    const continuation = { projectionGeometry: 'boundary-increment', ...(policy === undefined ? {} : { linearOrdering: policy }),
      ...(policy === 'station-auto' ? { stationFallback: true } : {}) };
    const result = { checkpoint: { version: 1, restart: { input,
      options: { transitionMode: 'automatic', ncrit: 4 }, initialEuler: { x: [0] }, initialBL: [0] }, continuation },
      families, mesh: { initialization: { gridSmoothing: { enabled: true } } } };
    const before = structuredClone(result), mapped = { input: structuredClone(input), options: { transitionMode: 'automatic', ncrit: 4 },
      initialEuler: { x: [1] }, initialBL: [2], system: { n: 2 }, diagnostics: { mapped: true } };
    const module = await stubbed('../src/euler/streamtube-transition-recovery.js', {
      createCoupledStreamtubeBody: () => ({ n: 1, initial: [0], evaluate: () => ({ families }) }),
      refineCoupledStreamtubeBody: () => mapped,
      solveCoupledStreamtubeIses: (target, options) => {
        calls.push({ input: target, options }); return { mesh: { initialization: {} } };
      }, streamtubeMeshSnapshot: () => { throw Error('Unexpected mesh evaluation'); },
    });
    module.recoverCoupledTransition(result, { plan: { normalFactor: 1 }, maxIterations: 0 });
    assert.equal(calls.length, 1); assert.equal(calls[0].options.linearOrdering, policy ?? 'auto');
    assert.equal(Object.hasOwn(calls[0].options, 'stationFallback'), false);
    assert.equal(Object.hasOwn(calls[0].options, 'resume'), false);
    assert.equal(calls[0].options.projectionGeometry, continuation.projectionGeometry);
    assert.equal(calls[0].input, mapped.input); assert.deepEqual(calls[0].input.hybrid, input.hybrid);
    assert.deepEqual(result, before);
  }
});

test('public paired refinement resumes the entire transferred ordering and fallback checkpoint without selecting a new policy', async () => {
  const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-refinement-assembly.js', import.meta.url), 'utf8');
  const solveFunction = source.slice(source.indexOf('export function solveCoupledStreamtubeRefinement('));
  assert.ok(solveFunction.startsWith('export function'));
  for (const policy of [undefined, 'station', 'station-auto', 'aligned-auto']) {
    const continuation = { iterationGeometry: 'ises-sampled', stepAcceptance: 'admissible', stagnationLimiter: 'listing',
      ...(policy === undefined ? {} : { linearOrdering: policy }), ...(policy === 'station-auto' ? { stationFallback: true } : {}) };
    const checkpoint = { continuation, restart: { options: {}, input: { bodies: [], hybrid: { ismom: 4 } } } }, before = structuredClone(checkpoint);
    const parent = { solverSettings: {}, initialization: { euler: {} } }, calls = [], key = `__pairedOrdering${serial++}`;
    globalThis[key] = { ...shearPolicy, coupledCheckpointHkPolicy,
      require: (condition, message) => { if (!condition) throw Error(message); },
      prepareCoupledStreamtubeRefinement: () => ({ checkpoint,
        settings: { normalization: {}, materialTrips: [], elementOrder: [] }, refinement: {}, plan: { tolerance: 1e-10 } }),
      solveCoupledStreamtubeIses: (input, options) => {
        calls.push({ input, options }); return { checkpoint, history: [{}], families: {},
          mesh: { quality: { valid: true } }, solverInput: checkpoint.restart.input };
      }, streamtubeMeshSnapshot: () => { throw Error('Unexpected mesh evaluation'); },
    };
    try {
      const module = await import(encode(`const {coupledResultShearCoordinate,coupledCheckpointHkPolicy,require,prepareCoupledStreamtubeRefinement,solveCoupledStreamtubeIses,streamtubeMeshSnapshot}=globalThis['${key}'];\n` + solveFunction));
      const result = module.solveCoupledStreamtubeRefinement({}, parent, { maxIterations: 0 });
      assert.equal(calls.length, 1); assert.equal(calls[0].input, undefined);
      assert.equal(calls[0].options.resume, checkpoint);
      assert.equal(Object.hasOwn(calls[0].options, 'linearOrdering'), false);
      assert.deepEqual(result.checkpoint, before); assert.deepEqual(checkpoint, before);
    } finally { delete globalThis[key]; }
  }
});
