import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { solveChainQuadratic } from '../src/numerics/chain-quadratic-program.js';
import { polishChainQuadratic } from '../src/numerics/chain-quadratic-polish.js';
import { fitConstrainedStations } from '../src/numerics/constrained-stations.js';

test('the exact three-element failed passage recovers a certified optimum without changing spacing tolerances', () => {
  const { input } = JSON.parse(fs.readFileSync(new URL('./fixtures/three-passage-quadratic.json', import.meta.url)));
  const before = structuredClone(input), output = fitConstrainedStations(input);
  assert.ok(output.sweeps < 17);
  for (const value of [output.primalResidual, output.stationarityResidual, output.complementarity]) assert.ok(value <= input.tolerance);
  const p = output.positions, h = p.slice(1).map((s, i) => s - p[i]);
  assert.ok(h.every(v => v > 0));
  assert.equal(h[0], input.firstSpacing);
  assert.ok(Math.abs(h.at(-1) - input.lastSpacing) < 1e-15);
  assert.equal(p.at(-1), 1);
  for (let i = 1; i < h.length; i++) assert.ok(Math.max(h[i] / h[i - 1], h[i - 1] / h[i]) <= input.maximumGrowth + 1e-12);
  assert.deepEqual(input, before);
  const reflected = fitConstrainedStations({ ...input, positions: input.positions.toReversed().map(s => 1 - s),
    firstSpacing: input.lastSpacing, lastSpacing: input.firstSpacing });
  assert.ok(Math.abs(reflected.objective - output.objective) < 1e-10);
  p.forEach((s, i) => assert.ok(Math.abs(s - (1 - reflected.positions.at(-1 - i))) < 1e-12));
});

test('active-chain polishing solves the original equality-constrained objective and retains caller arrays', () => {
  const constraints = [{ equality: true, terms: [[0, 1], [1, 1]], rhs: 1 },
    { terms: [[0, 1]], rhs: .2 }, { terms: [[1, -1]], rhs: 0 }];
  const input = { constraints, slack: [1e-13, .8], dual: [.6, 1e-13], equalityMultiplier: .2, tolerance: 1e-10 };
  const before = structuredClone(input), result = polishChainQuadratic(input);
  assert.ok(result); assert.deepEqual(result.x, [.2, .8]);
  assert.ok(Math.abs(result.equalityMultiplier - .2) < 1e-15);
  assert.ok(Math.abs(result.inequalityMultipliers[0] - .6) < 1e-15);
  assert.equal(result.inequalityMultipliers[1], 0);
  assert.deepEqual(input, before);
});

test('a guessed active set cannot pass with a negative multiplier or an omitted violated inequality', () => {
  const equality = { equality: true, terms: [[0, 1], [1, 1]], rhs: 1 };
  const controls = { slack: [1e-13], dual: [1], equalityMultiplier: 0, tolerance: 1e-10 };
  assert.equal(polishChainQuadratic({ ...controls, constraints: [equality, { terms: [[0, 1]], rhs: .9 }] }), null);
  assert.equal(polishChainQuadratic({ ...controls, slack: [1], constraints: [equality, { terms: [[0, 1]], rhs: .2 }] }), null);
});

test('dependent fixed chain endpoints preserve a valid nonnegative dual certificate', () => {
  const constraints = [{ equality: true, terms: [[0, 1], [1, 1]], rhs: .6 },
    { terms: [[0, 1]], rhs: .2 }, { terms: [[0, -2], [1, 1]], rhs: 0 }, { terms: [[1, -1]], rhs: -.4 }];
  const result = polishChainQuadratic({ constraints, slack: [1e-14, 1e-14, 1e-14], dual: [2, 1, .4], equalityMultiplier: 0, tolerance: 1e-10 });
  assert.ok(result); assert.deepEqual(result.x, [.2, .4]);
  assert.ok(result.dualResidual <= 1e-10); assert.ok(result.inequalityMultipliers.every(v => v >= 0));
});

test('nonhomogeneous adjacent constraints use the unchanged barrier path and exhausted iterations still fail', () => {
  const constraints = [{ equality: true, terms: [[0, 1], [1, 1]], rhs: 1 },
    { terms: [[0, 1], [1, -1]], rhs: -.2 }, { terms: [[0, -1]], rhs: 0 }];
  assert.equal(polishChainQuadratic({ constraints, slack: [1e-12, .4], dual: [1, 1e-12], equalityMultiplier: 0, tolerance: 1e-10 }), null);
  const result = solveChainQuadratic({ constraints });
  assert.ok(result.converged); assert.ok(Math.abs(result.x[0] - .4) < 1e-9);
  assert.equal(solveChainQuadratic({ constraints, maxIterations: 1 }).converged, false);
});
