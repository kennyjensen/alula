// SPDX-License-Identifier: GPL-2.0-or-later
// Prescribed scalar displacement coordinates only. Geometry and equation
// validity are checked by the caller, independently of this nodal chart.
export function createTransverseNodeChart({ nodes, directions }) {
  const nx = nodes?.length - 1, nt = nodes?.[0]?.length - 1;
  const shape = grid => Array.isArray(grid) && grid.length === nx + 1 && grid.every(row => Array.isArray(row)
    && row.length === nt + 1 && row.every(p => Number.isFinite(p?.x) && Number.isFinite(p?.y)));
  if (!(nx >= 2 && nt >= 2) || !shape(nodes) || !shape(directions)) throw new Error('Invalid transverse node chart.');
  const initial = Object.freeze(nodes.map(row => Object.freeze(row.map(p => Object.freeze({ x: p.x, y: p.y })))));
  const guide = Object.freeze(directions.map(row => Object.freeze(row.map(d => {
    const length = Math.hypot(d.x, d.y);
    if (!(length > 0) || !Number.isFinite(length)) throw new Error('A transverse guide must be nonzero.');
    return Object.freeze({ x: d.x / length, y: d.y / length });
  }))));
  let lengthScale = 0;
  for (const row of initial) for (const p of row) lengthScale = Math.max(lengthScale, Math.hypot(p.x - initial[0][0].x, p.y - initial[0][0].y));
  if (!(lengthScale > 0) || !Number.isFinite(lengthScale)) throw new Error('Degenerate transverse node chart.');
  const n = (nx - 1) * (nt - 1), index = (i, j) => !i || !j || i === nx || j === nt ? -1 : (i - 1) * (nt - 1) + j - 1;
  const validate = grid => {
    if (!shape(grid)) throw new Error('Invalid transverse grid state.');
    for (let i = 0; i <= nx; i++) for (let j = 0; j <= nt; j++) {
      const dx = grid[i][j].x - initial[i][j].x, dy = grid[i][j].y - initial[i][j].y;
      if (index(i, j) < 0 && (dx !== 0 || dy !== 0)) throw new Error('Transverse grid boundaries must remain fixed.');
      if (Math.abs(dx * guide[i][j].y - dy * guide[i][j].x) > 64 * Number.EPSILON * lengthScale)
        throw new Error('A transverse grid node left its prescribed guide line.');
    }
  };
  const move = (grid, delta, step = 1) => {
    validate(grid);
    if (delta?.length !== n || !Array.from(delta).every(Number.isFinite) || !Number.isFinite(step)) throw new Error('Invalid transverse grid direction.');
    return grid.map((row, i) => row.map((p, j) => {
      const k = index(i, j), d = guide[i][j];
      return k < 0 ? { ...p } : { x: p.x + step * delta[k] * d.x, y: p.y + step * delta[k] * d.y };
    }));
  };
  return { nx, nt, n, initial, directions: guide, lengthScale, validate, move };
}
