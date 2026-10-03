import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { requirePositiveSimpleQuad } from '../src/geometry/simple-quadrilateral.js';
import { quadLaplaceMatrix } from '../src/numerics/quad-laplace.js';
import { evaluateStreamtubeCell } from '../src/euler/streamtube-cell.js';
import { linearizeStreamtubeCell } from '../src/euler/streamtube-linearization.js';
import { solveStreamtubeIses } from '../src/euler/streamtube-ises-update.js';
import { createStreamtubeBodySystem } from '../src/euler/streamtube-body.js';
import { sparseProduct } from '../src/numerics/sparse.js';
import { intrinsicBodyFixture } from './fixtures/intrinsic-body.js';

const points = a => a.map(([x, y]) => ({ x, y }));
const concave = points([[0, 0], [1, 0], [.49, .49], [0, 1]]);
const close = (a, b, tolerance = 2e-8) => assert.ok(Math.abs(a - b) <= tolerance * Math.max(1, Math.abs(a), Math.abs(b)), `${a} != ${b}`);

test('explicit intermediate polygons may be concave, but never crossed, reversed or degenerate', () => {
  close(requirePositiveSimpleQuad(concave), .49, 1e-14);
  close(requirePositiveSimpleQuad(concave.map(p => ({ x: p.x + 4, y: p.y - 3 }))), .49, 1e-14);
  for (const vertices of [concave.toReversed(), points([[0, 0], [2, 1], [0, 1], [1, 0]]),
    points([[0, 0], [1, 0], [.5, 0], [0, 1]]), points([[0, 0], [1, 0], [1, 0], [0, 1]]),
    points([[0, 0], [1, 0], [2, 0], [3, 0]]), points([[0, 0], [NaN, 0], [1, 1], [0, 1]])])
    assert.throws(() => requirePositiveSimpleQuad(vertices));
});

test('sampled STIFF experiment checks Gauss determinants without claiming a valid Q1 map', () => {
  assert.throws(() => quadLaplaceMatrix(concave), /Invalid quadrilateral/);
  const matrix = quadLaplaceMatrix(concave, { quadratureOrder: 2, quadratureDomain: 'sampled-positive' });
  // Independent affine-field energy: its discrete integral equals signed
  // polygon area even though the concave corner invalidates the full Q1 map.
  for (const [gx, gy] of [[1, 0], [0, 1], [.3, -.8]]) {
    const u = concave.map(p => gx * p.x + gy * p.y); let energy = 0;
    for (let i = 0; i < 4; i++) {
      close(matrix.slice(4 * i, 4 * i + 4).reduce((a, b) => a + b, 0), 0, 1e-13);
      for (let j = 0; j < 4; j++) { close(matrix[4 * i + j], matrix[4 * j + i], 1e-13); energy += u[i] * matrix[4 * i + j] * u[j]; }
    }
    close(energy, .49 * (gx * gx + gy * gy), 1e-13);
  }
  const deeper = points([[0, 0], [1, 0], [.1, .1], [0, 1]]);
  assert.ok(requirePositiveSimpleQuad(deeper) > 0);
  assert.throws(() => quadLaplaceMatrix(deeper, { quadratureOrder: 2, quadratureDomain: 'sampled-positive' }), /Invalid quadrilateral/);
});

test('concave intermediate Euler cell retains independent flux balance and coordinate derivatives', () => {
  const p = { lower: points([[0, 0], [1, 0], [2, 0]]), upper: points([[0, 1], [.49, .49], [2, 1]]),
    densities: [1.05, 1.02], massFlow: .27, stagnationEnthalpy: 8, geometryDomain: 'positive-simple' };
  assert.throws(() => evaluateStreamtubeCell({ ...p, geometryDomain: 'convex' }), /Folded/);
  const { value: c, apply } = linearizeStreamtubeCell(p);
  const { geometry: g, states: [a, b], interfacePressure: pi } = c, flux = { x: 0, y: 0 };
  const pressure = (value, edge, sign) => { flux.x += sign * value * edge.y; flux.y -= sign * value * edge.x; };
  pressure(a.p, g.sections[0], -1); pressure(b.p, g.sections[1], 1);
  pressure(pi.lower, g.sides.lower, 1); pressure(pi.upper, g.sides.upper, -1);
  for (const key of ['x', 'y']) flux[key] += p.massFlow * (b.q * g.directions[1][key] - a.q * g.directions[0][key]);
  close((flux.x * g.transverse.x + flux.y * g.transverse.y) / g.area, 0, 1e-13);
  close((flux.x * g.streamwise.x + flux.y * g.streamwise.y) / g.area, c.streamwiseResidual, 1e-13);
  const outputs = c => [c.streamwiseResidual, c.isentropicResidual, c.interfacePressure.lower, c.interfacePressure.upper];
  for (const side of ['lower', 'upper']) for (let i = 0; i < 3; i++) for (const key of ['x', 'y']) {
    const tangent = points([[0, 0], [0, 0], [0, 0]]); tangent[i][key] = 1;
    const plus = structuredClone(p), minus = structuredClone(p), h = 1e-6;
    plus[side][i][key] += h; minus[side][i][key] -= h;
    const a = outputs(evaluateStreamtubeCell(plus)), b = outputs(evaluateStreamtubeCell(minus));
    outputs(apply({ [side]: tangent })).forEach((v, k) => close(v, (a[k] - b[k]) / (2 * h)));
  }
});

test('research intermediate-domain option preserves a controlled root and final convexity requirement', () => {
  const input = intrinsicBodyFixture({ elements: 2, bodySegments: 4, tubes: 2 });
  const strict = solveStreamtubeIses(input), sampled = solveStreamtubeIses(input, { iterationGeometry: 'ises-sampled' });
  assert.equal(strict.iterationGeometry, 'convex');
  assert.equal(sampled.residualConverged, true, sampled.reason);
  assert.equal(sampled.finalQuality.valid, true); assert.equal(sampled.converged, true);
  assert.deepEqual(sampled.x, strict.x); assert.deepEqual(sampled.nodes, strict.nodes);
  assert.throws(() => solveStreamtubeIses(input, { iterationGeometry: 'unknown' }), /controls/);
  assert.throws(() => solveStreamtubeIses({ ...input, geometryDomain: 'positive-simple' }), /Conflicting/);
});

test('a retained concave default state differentiates correctly but cannot pass final grid acceptance', () => {
  // This is a saved algebraic root, explicitly rejected for final geometry.
  // Reuse it for fast assembly/acceptance checks; do not rerun its trajectory.
  const row = JSON.parse(readFileSync(new URL('../docs/default-euler-ises-sampled.json', import.meta.url))).cases[0];
  const { input, initialEuler } = row.restart, system = createStreamtubeBodySystem(input);
  const state = system.adoptGeometry(Float64Array.from(initialEuler.x), initialEuler.nodes);
  assert.ok(system.evaluate(state).diagnostics.residual < 1e-10);
  // A five-step difference sweep found truncation above 1e-6 and pressure
  // subtraction roundoff below 1e-8; 1e-7 resolves this direction best.
  const matrix = system.jacobian(state, { sparse: true }), h = 1e-7;
  const direction = state.map((_, i) => Math.sin(i + .4)), product = sparseProduct(matrix, direction);
  const a = system.residual(state.map((v, i) => v + h * direction[i]));
  const b = system.residual(state.map((v, i) => v - h * direction[i]));
  for (let i = 0; i < matrix.n; i++) close(product[i], (a[i] - b[i]) / (2 * h), 2e-6);
  const strict = createStreamtubeBodySystem({ ...input, geometryDomain: 'convex' });
  const strictState = strict.adoptGeometry(Float64Array.from(initialEuler.x), initialEuler.nodes);
  assert.throws(() => strict.evaluate(strictState), /Folded/);
  // Deliberately loose residual tolerance isolates the independent grid gate
  // after the driver's initial redistribution changes this restart's residual.
  const result = solveStreamtubeIses(input, { initialEuler, iterationGeometry: 'ises-sampled', maxIterations: 0, tolerance: 1 });
  assert.equal(result.residualConverged, true);
  assert.equal(result.finalQuality.valid, false);
  assert.equal(result.converged, false);
  assert.match(result.reason, /final grid fails/);
});
