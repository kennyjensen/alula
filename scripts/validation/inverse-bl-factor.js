// SPDX-License-Identifier: GPL-2.0-or-later
// Exact station-block forward substitution for the conditional inverse BL
// Jacobian. This is a linear factor of the complete simultaneous equations.
import { factorLinear } from '../../src/numerics/linear.js';

export function factorInverseBL(system, matrix) {
  const { bl, ne, n } = system;
  if (matrix.n !== n) throw new Error('Inverse BL matrix dimension mismatch.');
  const te = new Set(bl.wakes.map(w => w.ids[0])), owner = new Int32Array(n).fill(-1);
  const blocks = bl.stations.map(({ id }, order) => {
    const columns = (te.has(id) ? [0, 1, 2] : [0, 1, 3]).map(k => ne + 4 * id + k);
    columns.forEach(col => { owner[col] = order; });
    return { id, columns, rows: [0, 1, 2].map(k => ne + 4 * id + k) };
  });
  blocks.forEach((block, order) => {
    const local = new Float64Array(9), links = [];
    block.rows.forEach((row, i) => {
      const dependencies = []; links.push(dependencies);
      for (let k = matrix.rowPtr[row]; k < matrix.rowPtr[row + 1]; k++) {
        const col = matrix.colIndex[k], v = matrix.values[k]; if (!v) continue;
        const j = block.columns.indexOf(col);
        if (j >= 0) local[3 * i + j] = v;
        else {
          if (owner[col] >= order) throw new Error('Inverse BL dependency is outside its station order.');
          dependencies.push({ col, value: v });
        }
      }
    });
    Object.assign(block, { local, links, solve: factorLinear(local, 3) });
  });
  const lift = (retained, rhs = new Float64Array(n)) => {
    if (retained.length !== n || rhs.length !== n || !retained.every(Number.isFinite) || !rhs.every(Number.isFinite))
      throw new Error('Invalid inverse BL linear right-hand side.');
    const x = Float64Array.from(retained);
    for (const block of blocks) {
      const work = Float64Array.from(block.rows, (row, i) => block.links[i].reduce((r, p) => r - p.value * x[p.col], rhs[row]));
      const a = block.solve(work); block.columns.forEach((col, j) => { x[col] = a[j]; });
    }
    return x;
  };
  return { lift, blocks, owner, retained: Array.from(owner.keys()).filter(i => owner[i] < 0) };
}
