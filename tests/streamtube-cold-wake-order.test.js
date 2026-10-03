import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-initializer.js', import.meta.url), 'utf8');
let serial = 0;
async function order(initialization, options = {}, initialThicknessFactor = 1, recovery = {}) {
  const key = `__mrchueWarningGate${++serial}`, rows = [], input = { source: 'unchanged' };
  const controls = { initialEuler: {}, transitionMode: 'automatic', ...options }, before = structuredClone(controls);
  let constructors = 0;
  const forbidden = () => { throw new Error('No geometry/flow evaluation is allowed in the gate seam.'); };
  const seed = { ne: 0, initial: Float64Array.of(.1, 1, 3, 1, .1, 1, 2, 1), initialization,
    bl: { transitionMode: controls.transitionMode, snapshotActive: () => [], scale: 1, wakes: [{ ids: [0, 1] }] } };
  globalThis[key] = {
    createCoupledStreamtubeBody(_, controls) {
      if (constructors++ === 0) return seed;
      if (controls.blInitialization === 'direct-inverse') {
        recovery.attempts = (recovery.attempts ?? 0) + 1;
        recovery.afterCandidates = rows.length;
        if (recovery.resolve) return { ...seed, initialization: { boundaryLayer: [] } };
      }
      throw new Error('Controlled rejection before geometry: inspect the complete guess order.');
    },
    extendStreamtubeDisplacement: forbidden, streamtubeMeshSnapshot: forbidden,
    initializeStreamtubeWakeCorrespondence: forbidden,
  };
  const code = source.replace(/^import \{([^}]+)\} from '[^']+';/gm, (_, names) => `const {${names}}=globalThis.${key};`);
  try {
    const module = await import('data:text/javascript;base64,' + Buffer.from(code + `\n// ${serial}`).toString('base64'));
    assert.throws(() => module.initializeCoupledStreamtubeBody(input, controls, {
      maximumBacktracks: 1, initialThicknessFactor, onAttempt: row => rows.push(row),
    }), /Coupled initialization failed/);
    assert.deepEqual(controls, before);
    return rows.map(row => ({ thicknessFactor: row.thicknessFactor, wakeInitialization: row.wakeInitialization }));
  } finally { delete globalThis[key]; }
}
const warning = { boundaryLayer: [{ method: 'mrchue', localConvergenceWarnings: ['MRCHUE: Convergence failed at a local surface station'] }] };
const regular = [{ thicknessFactor: 1, wakeInitialization: undefined },
  { thicknessFactor: 1, wakeInitialization: 'iset-linear-shape' }, { thicknessFactor: .5, wakeInitialization: undefined }];

test('only an explicit unresolved MRCHUE surface march moves the existing ISET guess first', async () => {
  assert.deepEqual(await order(warning), [regular[1], regular[0], regular[2]]);
  assert.deepEqual(await order({ boundaryLayer: [{ method: 'mrchue', localConvergenceWarnings: [] },
    warning.boundaryLayer[0]] }), [regular[1], regular[0], regular[2]]);
});

test('an unfinished surface is resolved once after both full-thickness guesses fail, before thinning', async () => {
  const recovery = { resolve: true };
  assert.deepEqual(await order(warning, {}, 1, recovery), [regular[1], regular[0], regular[0], regular[2]]);
  assert.equal(recovery.attempts, 1);
  assert.equal(recovery.afterCandidates, 2);
  const rejected = {};
  assert.deepEqual(await order(warning, {}, 1, rejected), [regular[1], regular[0], regular[2]]);
  assert.equal(rejected.attempts, 1, 'An unresolved alternate must retain the original bounded guesses.');
});

test('complete, missing, non-MRCHUE and malformed warning metadata keep the original order', async () => {
  for (const initialization of [undefined, {}, { boundaryLayer: [] }, { boundaryLayer: {} },
    { boundaryLayer: [null] }, { boundaryLayer: [{ method: 'mrchue' }] },
    { boundaryLayer: [{ method: 'mrchue', localConvergenceWarnings: [] }] },
    { boundaryLayer: [{ method: 'mrchue', localConvergenceWarnings: 'Convergence failed' }] },
    { boundaryLayer: [{ method: 'direct-inverse', localConvergenceWarnings: ['Convergence failed'] }] }])
    assert.deepEqual(await order(initialization), regular);
});

test('fixed transition, supplied BL and explicitly thinner retries retain their previous candidates', async () => {
  const ordinaryOnly = [regular[0], regular[2]];
  assert.deepEqual(await order(warning, { transitionMode: 'fixed-trip' }), ordinaryOnly);
  assert.deepEqual(await order(warning, { initialBL: [1, 2, 3, 4] }), ordinaryOnly);
  assert.deepEqual(await order(warning, {}, .25), [
    { thicknessFactor: .25, wakeInitialization: undefined }, { thicknessFactor: .125, wakeInitialization: undefined },
  ]);
});
