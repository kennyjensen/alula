// SPDX-License-Identifier: GPL-2.0-or-later
// Exact algebraic BL elimination. Keep every Euler/grid unknown and each
// mapped displacement thickness (TE matching retains wake Ue instead).
import { factorInverseBL } from './inverse-bl-factor.js';
import { sparseMatrix, sparseAdd } from '../../src/numerics/sparse.js';

export function reduceInverseBL(system, matrix, rhs) {
  const { n, ne, bl } = system, factor = factorInverseBL(system, matrix);
  const relations = new Array(n), add = (row, col, v) => { if (v) row.set(col, (row.get(col) ?? 0) + v); };
  let relationNonzeros = 0;
  for (const block of factor.blocks) {
    const work = block.links.map(links => {
      const row = new Map();
      for (const { col, value } of links) {
        if (factor.owner[col] < 0) add(row, col, -value);
        else for (const [i, v] of relations[col]) add(row, i, -value * v);
      }
      return row;
    });
    const maps = block.columns.map(() => new Map()), columns = new Set(work.flatMap(row => [...row.keys()]));
    for (const col of columns) {
      const values = block.solve(work.map(row => row.get(col) ?? 0));
      values.forEach((v, i) => { if (v) maps[i].set(col, v); });
    }
    block.columns.forEach((col, i) => { relations[col] = maps[i]; relationNonzeros += maps[i].size; });
  }
  const rows = [...Array(ne).keys(), ...bl.stations.map(p => ne + 4 * p.id + 3)], columns = factor.retained;
  if (rows.length !== columns.length) throw new Error('Inverse BL reduction is not square.');
  const indices = new Int32Array(n).fill(-1); columns.forEach((col, i) => { indices[col] = i; });
  const constants = factor.lift(new Float64Array(n), rhs), reducedRhs = new Float64Array(rows.length);
  const coefficients = rows.map((row, i) => {
    const result = new Map(); let b = rhs[row];
    for (let k = matrix.rowPtr[row]; k < matrix.rowPtr[row + 1]; k++) {
      const col = matrix.colIndex[k], v = matrix.values[k]; if (!v) continue;
      if (factor.owner[col] < 0) add(result, indices[col], v);
      else {
        b -= v * constants[col];
        for (const [c, weight] of relations[col]) add(result, indices[c], v * weight);
      }
    }
    reducedRhs[i] = b; return result;
  });
  const reduced = sparseMatrix(coefficients.map(row => row.keys()));
  coefficients.forEach((row, i) => row.forEach((v, col) => sparseAdd(reduced, i, col, v)));
  const lift = y => {
    if (y.length !== columns.length) throw new Error('Wrong reduced direction length.');
    const retained = new Float64Array(n); columns.forEach((col, i) => { retained[col] = y[i]; });
    return factor.lift(retained, rhs);
  };
  return { matrix: reduced, rhs: reducedRhs, lift, rows, columns, relationNonzeros };
}
