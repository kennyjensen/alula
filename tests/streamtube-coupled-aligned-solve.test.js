// SPDX-License-Identifier: GPL-2.0-or-later
import fs from 'node:fs';
import test from 'node:test';
import assert from 'node:assert/strict';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';
import { createCoupledStreamtubeBody } from '../src/euler/streamtube-coupled.js';
import { solveCoupledAlignedSystem } from '../src/euler/streamtube-coupled-aligned-solve.js';
import { solveCoupledLinearSystem } from '../src/euler/streamtube-coupled-linear-solve.js';
import { solveSparseDirectAligned } from '../src/numerics/klu.js';
import { createStreamtubeStationOrdering } from '../src/euler/streamtube-station-ordering.js';
import { sparseProduct } from '../src/numerics/sparse.js';

function fixture() {
  const system = createCoupledStreamtubeBody(intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 3 }),
    { edgeMatching: 'section-velocity', transitionMode: 'automatic' });
  const matrix = system.jacobian(system.initial);
  const exact = Float64Array.from({ length: matrix.n }, (_, i) => Math.sin(.37 * i));
  return { matrix, exact, rhs: sparseProduct(matrix, exact),
    controls: { layout: system.euler.layout, stations: system.bl.stations, preferredOrdering: 'amd', pivotTolerance: .001 } };
}

test('matched AMD solves the full coupled matrix in original coordinates without changing it', () => {
  const { matrix, rhs, exact, controls } = fixture(), before = structuredClone({ matrix, rhs });
  const result = solveCoupledAlignedSystem(matrix, rhs, controls);
  assert.equal(result.linear.equationOrdering, 'aligned');
  assert.equal(result.stationPolicy.mode, 'aligned-auto');
  assert.equal(result.stationPolicy.fallback, false);
  assert(result.linear.relativeResidual <= 1e-10);
  assert(Math.max(...result.linear.x.map((v, i) => Math.abs(v - exact[i]))) < 1e-6);
  assert.deepEqual({ matrix, rhs }, before);
});

let serial = 0;
async function injected(overrides) {
  const key = `__aligned${serial++}`;
  globalThis[key] = { solveSparseDirectAligned, createStreamtubeStationOrdering, solveCoupledLinearSystem, ...overrides };
  const source = fs.readFileSync(new URL('../src/euler/streamtube-coupled-aligned-solve.js', import.meta.url), 'utf8')
    .replace(/import \{([^}]+)\} from '[^']+';/g, (_, names) => `const {${names}} = globalThis[${JSON.stringify(key)}];`);
  try { return (await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`)).solveCoupledAlignedSystem; }
  finally { delete globalThis[key]; }
}

test('numerical alignment rejection falls back to the established policy on the identical system', async () => {
  const f = fixture(); let calls = 0;
  const solve = await injected({ solveSparseDirectAligned() { throw Object.assign(new Error('accuracy'), { code: 'KLU_RESIDUAL_LIMIT' }); },
    solveCoupledLinearSystem(a, b, options) { calls++; assert.equal(a, f.matrix); assert.equal(b, f.rhs);
      assert.equal(options.mode, 'station-auto'); return solveCoupledLinearSystem(a, b, options); } });
  const result = solve(f.matrix, f.rhs, f.controls);
  assert.equal(calls, 1); assert.equal(result.stationPolicy.fallback, true);
  assert.equal(result.stationPolicy.alignedFailure.code, 'KLU_RESIDUAL_LIMIT');
  assert(result.linear.relativeResidual <= 1e-10);
});

test('cancellation and unexpected errors are not retried', async () => {
  const f = fixture();
  for (const error of [Object.assign(new Error('cancelled'), { name: 'AbortError' }), new Error('unexpected')]) {
    const solve = await injected({ solveSparseDirectAligned() { throw error; },
      solveCoupledLinearSystem() { assert.fail('must not retry'); } });
    assert.throws(() => solve(f.matrix, f.rhs, f.controls), e => e === error);
  }
});
