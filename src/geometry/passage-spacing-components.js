// SPDX-License-Identifier: GPL-2.0-or-later
// Spacing can propagate through a cut, but not through a solid body.
// Build connected components of actual passage blocks; do not split an
// exterior block at a LE/TE event shielded by its neighboring solid wall.
export function createPassageSpacingComponents(blocks) {
  const { passages, rank, inletRank, outletRank, constraints } = blocks ?? {};
  const bodies = constraints?.upper?.length;
  if (!Number.isInteger(bodies) || bodies < 1 || !Array.isArray(passages) || passages.length !== bodies + 1
    || !rank || !Number.isFinite(inletRank) || !Number.isFinite(outletRank) || !(inletRank < outletRank))
    throw new Error('Spacing components require an ordered panel passage skeleton.');
  const extents = Array.from({ length: bodies }, (_, b) => [rank[`${b}:LE`], rank[`${b}:TE`]]);
  if (extents.some(([le, te]) => ![le, te].every(Number.isFinite) || !(inletRank < le && le < te && te < outletRank)))
    throw new Error('Spacing components require ordered body endpoints.');
  const pieces = passages.flatMap((events, passage) => {
    if (!Array.isArray(events)) throw new Error('Invalid passage events.');
    const bounds = [inletRank, ...events.map(e => rank[e.id]), outletRank];
    if (bounds.some((r, i) => !Number.isFinite(r) || i && !(r > bounds[i - 1]))) throw new Error('Passage spacing bounds must increase.');
    return bounds.slice(1).map((end, k) => ({ passage, start: bounds[k], end }));
  });
  const bodyState = (b, start, end) => {
    const [le, te] = extents[b];
    if (end <= le) return 'upstream'; if (start >= te) return 'wake';
    if (start >= le && end <= te) return 'wall';
    throw new Error('A spacing block crosses an unsplit body endpoint.');
  };
  const byKey = new Map(pieces.map(p => [`${p.passage}:${p.start}:${p.end}`, p]));
  const visited = new Set(), components = [];
  for (const piece of pieces) {
    if (visited.has(piece)) continue;
    const { start, end } = piece, members = [], stack = [piece];
    while (stack.length) {
      const p = stack.pop(); if (visited.has(p)) continue;
      visited.add(p); members.push(p.passage);
      for (const direction of [-1, 1]) {
        const b = direction < 0 ? p.passage - 1 : p.passage;
        if (b < 0 || b >= bodies || bodyState(b, start, end) === 'wall') continue;
        const neighbor = byKey.get(`${p.passage + direction}:${start}:${end}`);
        if (!neighbor) throw new Error('The two banks of a shared cut need identical actual block bounds.');
        if (!visited.has(neighbor)) stack.push(neighbor);
      }
    }
    members.sort((a, b) => a - b);
    const first = members[0], last = members.at(-1), boundaries = [];
    boundaries.push(first === 0 ? { kind: 'outer', side: 'lower' } : { kind: 'wall', body: first - 1, side: 'upper' });
    for (let b = first; b < last; b++) boundaries.push({ kind: 'cut', body: b, end: bodyState(b, start, end) });
    boundaries.push(last === bodies ? { kind: 'outer', side: 'upper' } : { kind: 'wall', body: last, side: 'lower' });
    components.push({ start, end, passages: members, boundaries });
  }
  return components.sort((a, b) => a.start - b.start || a.end - b.end || a.passages[0] - b.passages[0]);
}
