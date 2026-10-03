// SPDX-License-Identifier: GPL-2.0-or-later
// Serialize already accumulated square row Maps directly to FP64 CSR.
// Equivalent on valid input to sparseMatrix(rows.map(r => r.keys())) followed
// by one sparseAdd per Map entry: explicit zeros stay, missing diagonals are
// forced, and each stored value is 0 + d (normalizing negative zero).
export function sparseMatrixFromRows(rows) {
  const n = rows?.length;
  if (!Array.isArray(rows) || !Number.isInteger(n) || n < 1 || n >= 0x80000000)
    throw new Error('Sparse rows require a nonempty square array of row Maps.');
  const rowPtr = new Int32Array(n + 1), sortedRows = new Array(n);
  let nonzeros = 0;
  for (let row = 0; row < n; row++) {
    const values = rows[row];
    if (!(values instanceof Map)) throw new Error(`Sparse row ${row} must be a Map.`);
    const columns = Array.from(values.keys());
    for (const col of columns) {
      if (!Number.isInteger(col) || col < 0 || col >= n)
        throw new Error(`Sparse row ${row} has a column outside the square matrix.`);
      if (!Number.isFinite(values.get(col))) throw new Error(`Sparse row ${row} has a nonfinite or nonnumeric value.`);
    }
    if (!values.has(row)) columns.push(row);
    columns.sort((a, b) => a - b);
    nonzeros += columns.length;
    if (nonzeros >= 0x80000000) throw new Error('Sparse row entries exceed signed 32-bit CSR indexing.');
    rowPtr[row + 1] = nonzeros; sortedRows[row] = columns;
  }
  const colIndex = new Int32Array(nonzeros), values = new Float64Array(nonzeros);
  for (let row = 0; row < n; row++) {
    let index = rowPtr[row];
    for (const col of sortedRows[row]) {
      colIndex[index] = col;
      values[index++] = 0 + (rows[row].get(col) ?? 0);
    }
  }
  return { n, rowPtr, colIndex, values };
}
