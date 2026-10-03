// SPDX-License-Identifier: GPL-2.0-or-later
// Assemble the complete intrinsic channel Jacobian from local chain rules.
// Adjacent cells contribute opposite interface-pressure derivatives to the
// same normal-momentum row. Inlet thermodynamics and end slopes are included.
import { linearizeStreamtubeCell } from '../streamtube-linearization.js';
import { sparseMatrix, sparseAdd } from '../../numerics/sparse.js';

export function createStreamtubeChannelJacobian(system) {
  const { n, nx, nt, densityCount, positionIndex, decode, conditions } = system;
  const { massFlows, stagnationEnthalpy: h0, gamma, referencePressure, heightScale, pressureCorrectionFactor, streamwiseMode } = conditions;
  const stencils = [], columns = Array.from({ length: n }, () => new Set());
  for (let i = 1; i < nx; i++) for (let j = 0; j < nt; j++) {
    const streamwiseRow = (i - 1) * nt + j;
    const inletRow = i === 1 ? (nx - 1) * nt + j : null;
    const upperRow = j < nt - 1 ? densityCount + (i - 1) * (nt - 1) + j : null;
    const lowerRow = j > 0 ? densityCount + (i - 1) * (nt - 1) + j - 1 : null;
    const variables = [0, 1].map(slot => ({ col: (i - 1 + slot) * nt + j, slot }));
    for (const [side, streamline] of [['lower', j], ['upper', j + 1]]) if (streamline > 0 && streamline < nt)
      for (let node = 0; node < 3; node++) variables.push({ col: positionIndex(i - 1 + node, streamline), side, node });
    const rows = [streamwiseRow, inletRow, upperRow, lowerRow].filter(row => row !== null);
    for (const row of rows) for (const { col } of variables) columns[row].add(col);
    stencils.push({ i, j, streamwiseRow, inletRow, upperRow, lowerRow, variables });
  }
  const endRows = [];
  let row = densityCount + (nx - 1) * (nt - 1);
  for (const i of [0, nx - 1]) for (let j = 1; j < nt; j++) {
    const left = positionIndex(i, j), right = positionIndex(i + 1, j);
    columns[row].add(left); columns[row].add(right); endRows.push({ row: row++, left, right });
  }
  if (row !== n) throw new Error('Intrinsic channel Jacobian row count mismatch.');
  const pattern = sparseMatrix(columns);

  return (state, { sparse = false } = {}) => {
    const { nodes, densities } = decode(state);
    const matrix = sparse ? { ...pattern, values: new Float64Array(pattern.values.length) } : new Float64Array(n * n);
    const add = sparse ? (row, col, v) => sparseAdd(matrix, row, col, v) : (row, col, v) => { matrix[row * n + col] += v; };
    for (const { i, j, streamwiseRow, inletRow, upperRow, lowerRow, variables } of stencils) {
      const cell = linearizeStreamtubeCell({ lower: [nodes[i - 1][j], nodes[i][j], nodes[i + 1][j]],
        upper: [nodes[i - 1][j + 1], nodes[i][j + 1], nodes[i + 1][j + 1]],
        densities: [densities[i - 1][j], densities[i][j]], massFlow: massFlows[j],
        stagnationEnthalpy: h0[j], gamma, pressureCorrectionFactor });
      if (cell.value.states.some(s => s.machSquared >= 1)) throw new Error('Intrinsic channel currently requires subsonic flow.');
      for (const { col, slot, side, node } of variables) {
        const tangent = {};
        if (side) {
          tangent[side] = Array.from({ length: 3 }, () => ({ x: 0, y: 0 }));
          tangent[side][node].y = heightScale;
        } else {
          tangent.densities = [0, 0]; tangent.densities[slot] = densities[i - 1 + slot][j];
        }
        const d = cell.apply(tangent);
        add(streamwiseRow, col, (streamwiseMode === 'isentropic' ? d.isentropicResidual : d.streamwiseResidual) / referencePressure);
        if (upperRow !== null) add(upperRow, col, d.interfacePressure.upper / referencePressure);
        if (lowerRow !== null) add(lowerRow, col, -d.interfacePressure.lower / referencePressure);
        if (inletRow !== null) {
          const s = cell.value.states[0], ds = d.states[0];
          add(inletRow, col, ds.rho / s.rho - ds.enthalpy / ((gamma - 1) * s.enthalpy));
        }
      }
    }
    // Encoded y offsets are divided by heightScale in the residual, so the
    // two end-segment coefficients are exactly +/-1, including on inclined ends.
    for (const { row, left, right } of endRows) { add(row, left, -1); add(row, right, 1); }
    return matrix;
  };
}
