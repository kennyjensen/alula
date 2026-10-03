import test from 'node:test';
import assert from 'node:assert/strict';
import { solveBlockTridiagonal2 } from '../src/numerics/block-tridiagonal-2.js';
import { solveLinear } from '../src/numerics/linear.js';

test('general block Thomas matches a pivoted dense solve for nonsymmetric lines at several physical scales', () => {
  for (const n of [1, 2, 9]) for (const scale of [1e-10, 1, 1e10]) {
    const blocks = {
      lower: Array.from({ length: n }, (_, k) => [1, .3, -.2, 1.5].map(v => v * scale)),
      diagonal: Array.from({ length: n }, (_, k) => [-8 - k, 1.2, -.7, -6 - k].map(v => v * scale)),
      upper: Array.from({ length: n }, (_, k) => [1.2, -.6, .2, .9].map(v => v * scale)),
      rhs: Array.from({ length: n }, (_, k) => [Math.sin(k + 1), Math.cos(2 * k)].map(v => v * scale))
    };
    const original = structuredClone(blocks), a = new Float64Array(4 * n * n);
    for (let k = 0; k < n; k++) for (const [name, offset] of [['lower', -1], ['diagonal', 0], ['upper', 1]]) {
      const j = k + offset; if (j < 0 || j >= n) continue;
      for (let r = 0; r < 2; r++) for (let c = 0; c < 2; c++) a[(2 * k + r) * 2 * n + 2 * j + c] = blocks[name][k][2 * r + c];
    }
    const expected = solveLinear(a, blocks.rhs.flat()), actual = solveBlockTridiagonal2(blocks).flat();
    actual.forEach((v, i) => assert.ok(Math.abs(v - expected[i]) < 2e-14));
    assert.deepEqual(blocks, original);
  }
});

test('block pivot may be indefinite but unresolved elimination pivots fail explicitly', () => {
  const zero = [0, 0, 0, 0];
  assert.deepEqual(solveBlockTridiagonal2({ lower: [zero], diagonal: [[0, 2, 3, 0]], upper: [zero], rhs: [[4, 9]] }), [[3, 2]]);
  assert.throws(() => solveBlockTridiagonal2({ lower: [zero, [1, 0, 0, 1]], diagonal: [[1, 0, 0, 1], [1, 0, 0, 1]],
    upper: [[1, 0, 0, 1], zero], rhs: [[1, 2], [3, 4]] }), /pivot/);
});
