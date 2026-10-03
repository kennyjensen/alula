// SPDX-License-Identifier: GPL-2.0-or-later
// Local to ONE Jacobian and its fixed geometry, density and capture state.
// A tube retains at most four section maps and two transport-speed maps.
// The caller must not mutate these inputs during the assembly.
const zero = Object.freeze({ x: 0, y: 0 });
const invalid = () => { throw new Error('Invalid Jacobian transport cache topology or query.'); };

export function createStreamtubeJacobianTransport({ chains, geometry, layout, massDerivatives } = {}) {
  if (!Number.isInteger(layout?.nx) || layout.nx < 2 || typeof layout.densityIndex !== 'function'
    || !Array.isArray(layout.tubes) || !layout.tubes.length
    || ![chains, geometry, massDerivatives].every(a => Array.isArray(a) && a.length === layout.tubes.length)) invalid();
  const { nx, tubes } = layout;
  const cache = tubes.map((count, g) => {
    if (!Number.isInteger(count) || count < 1 || !Array.isArray(chains[g]) || chains[g].length !== count
      || !Array.isArray(massDerivatives[g]) || massDerivatives[g].length !== count
      || !Array.isArray(geometry[g]) || geometry[g].length !== nx + 1
      || geometry[g].some(row => !Array.isArray(row) || row.length !== count + 1 || row.some(m => !(m instanceof Map)))) invalid();
    return chains[g].map((chain, j) => {
      if (chain?.sections?.length !== nx || typeof chain.sectionTangent !== 'function'
        || typeof chain.transportTangent !== 'function' || !(massDerivatives[g][j] instanceof Map)) invalid();
      return { i: 0, sections: new Map(), speeds: new Map() };
    });
  });
  const forCell = (i, g, j, col) => {
    if (!Number.isInteger(i) || i < 1 || i >= nx || !Number.isInteger(g) || g < 0 || g >= tubes.length
      || !Number.isInteger(j) || j < 0 || j >= tubes[g] || !Number.isInteger(col) || col < 0 || col >= 0x80000000) invalid();
    const c = cache[g][j], chain = chains[g][j];
    if (i < c.i) throw new Error('Jacobian transport cell queries must advance monotonically within each tube.');
    if (i > c.i) {
      for (const k of c.sections.keys()) if (k < i - 3) c.sections.delete(k);
      for (const k of c.speeds.keys()) if (k < i - 1) c.speeds.delete(k);
      c.i = i;
    }
    const at = k => {
      if (!Number.isInteger(k) || k < Math.max(0, i - 3) || k > i) invalid();
      if (!c.sections.has(k)) c.sections.set(k, new Map());
      const row = c.sections.get(k);
      if (!row.has(col)) row.set(col, chain.sectionTangent(k, {
        lower: [geometry[g][k][j].get(col) ?? zero, geometry[g][k + 1][j].get(col) ?? zero],
        upper: [geometry[g][k][j + 1].get(col) ?? zero, geometry[g][k + 1][j + 1].get(col) ?? zero],
        density: layout.densityIndex(k, g, j) === col ? chain.sections[k].rho : 0,
        massFlow: massDerivatives[g][j].get(col) ?? 0,
      }));
      return row.get(col);
    };
    const speed = k => {
      if (!c.speeds.has(k)) c.speeds.set(k, new Map());
      const row = c.speeds.get(k);
      if (!row.has(col)) row.set(col, chain.transportTangent(k, at));
      return row.get(col);
    };
    // Fresh pair: a caller may modify the returned array without corrupting
    // either cached scalar. Map.has preserves explicit zero and signed zero.
    return [speed(i - 1), speed(i)];
  };
  return { forCell };
}
