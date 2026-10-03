// SPDX-License-Identifier: GPL-2.0-or-later
// A necessary ordering constraint for a panel-isopotential block skeleton.
// Each passage has its own potential branch. Never sort all body TE events
// by the mean of their upper/lower potentials: circulation changes their
// relative order on facing surfaces. This is a topology utility, not the
// recovered MSET reconciliation algorithm or a geometric grid constructor.
export function orderPotentialBlockEvents(passages) {
  if (!Array.isArray(passages) || !passages.length)
    throw new Error('Potential blocks require passage events.');
  const vertices = new Map(), edges = new Map(), incoming = new Map();
  const orderedPassages = passages.map((events, passage) => {
    if (!Array.isArray(events) || events.length < 2)
      throw new Error('Each potential passage requires at least two events.');
    const seen = new Set();
    const sorted = events.map(event => {
      if (!event || !Number.isInteger(event.body) || event.body < 0
        || !['LE', 'TE'].includes(event.edge) || !Number.isFinite(event.potential))
        throw new Error('Invalid potential block event.');
      const id = `${event.body}:${event.edge}`;
      if (seen.has(id)) throw new Error('Duplicate potential event in a passage.');
      seen.add(id);
      if (!vertices.has(id)) {
        vertices.set(id, { id, body: event.body, edge: event.edge });
        edges.set(id, new Set()); incoming.set(id, 0);
      }
      return { ...event, id };
    }).sort((a, b) => a.potential - b.potential);
    for (let i = 1; i < sorted.length; i++) {
      const a = sorted[i - 1], b = sorted[i];
      if (a.potential === b.potential)
        throw new Error(`Coincident potential block events in passage ${passage} require an explicit shared block.`);
      if (!edges.get(a.id).has(b.id)) {
        edges.get(a.id).add(b.id); incoming.set(b.id, incoming.get(b.id) + 1);
      }
    }
    for (const event of sorted) if (event.edge === 'LE') {
      const trailing = sorted.find(e => e.body === event.body && e.edge === 'TE');
      if (trailing && trailing.potential < event.potential)
        throw new Error('Surface potential must increase from stagnation to the trailing edge.');
    }
    return sorted;
  });
  // A stable linear extension of passage-local orders. Ties here concern
  // unrelated events, so no comparison of potentials on different branches
  // is meaningful or needed. A cycle needs a different block/index topology.
  const compare = (a, b) => vertices.get(a).body - vertices.get(b).body
    || Number(vertices.get(a).edge === 'TE') - Number(vertices.get(b).edge === 'TE');
  const ready = [...vertices.keys()].filter(id => incoming.get(id) === 0).sort(compare), order = [];
  while (ready.length) {
    const id = ready.shift(); order.push(vertices.get(id));
    for (const next of edges.get(id)) {
      incoming.set(next, incoming.get(next) - 1);
      if (incoming.get(next) === 0) ready.push(next);
    }
    ready.sort(compare);
  }
  if (order.length !== vertices.size)
    throw new Error('Conflicting passage potential orders require a different block/index topology.');
  return { order, passages: orderedPassages, rank: Object.fromEntries(order.map((event, i) => [event.id, i])) };
}
