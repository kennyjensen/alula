// SPDX-License-Identifier: GPL-2.0-or-later
// Scalar reconstruction of LE/TE isopotential block constraints for a
// bottom-to-top stack of lifting bodies. This is not the recovered MSET
// algorithm or a geometric skeleton: shielded boundary targets stay unknown.
import { orderPotentialBlockEvents } from './potential-block-order.js';

export function createPanelPotentialBlocks({ bodies, outer } = {}) {
  if (!Array.isArray(bodies) || !bodies.length || !Array.isArray(outer) || outer.length !== 2)
    throw new Error('Panel potential blocks require bodies and two outer boundaries.');
  const bodyData = bodies.map((body, b) => {
    if (!body || !body.trailing
      || ![body.leading, body.trailing.upper, body.trailing.lower, body.inlet, body.outletIncrement].every(Number.isFinite))
      throw new Error('Invalid body potential block data.');
    if (!(body.inlet < body.leading && body.leading < body.trailing.upper
      && body.leading < body.trailing.lower && body.outletIncrement > 0))
      throw new Error('Body potentials must increase from inlet through LE and TE to outlet.');
    const outlets = { lower: body.trailing.lower + body.outletIncrement,
      upper: body.trailing.upper + body.outletIncrement };
    if (!Object.values(outlets).every(Number.isFinite)
      || !(outlets.lower > body.trailing.lower && outlets.upper > body.trailing.upper))
      throw new Error('Body outlet potential increment is not representable.');
    return { leading: body.leading, trailing: { ...body.trailing }, inlet: body.inlet, outlets, body: b };
  });
  const outerData = outer.map(boundary => {
    if (!boundary || !Number.isFinite(boundary.inlet) || !Number.isFinite(boundary.outlet)
      || !(boundary.inlet < boundary.outlet))
      throw new Error('Invalid outer potential block interval.');
    return { inlet: boundary.inlet, outlet: boundary.outlet };
  });
  const count = bodyData.length, passages = Array.from({ length: count + 1 }, () => []);
  const checkInterior = (potential, inlet, outlet, where) => {
    if (!Number.isFinite(potential) || potential < inlet || potential > outlet)
      throw new Error('Potential block target is outside ' + where + '.');
    if (potential === inlet || potential === outlet)
      throw new Error('Coincident potential block target at ' + where + ' endpoint requires an explicit shared block.');
  };
  const events = bodyData.flatMap((body, b) => ['LE', 'TE'].map(edge => {
    const event = { id: b + ':' + edge, body: b, edge,
      targets: { lower: Array(count).fill(null), upper: Array(count).fill(null), outer: [null, null] },
      passagePotentials: Array(count + 1).fill(null) };
    for (const direction of [-1, 1]) {
      const sourceSide = direction < 0 ? 'lower' : 'upper';
      let potential = edge === 'LE' ? body.leading : body.trailing[sourceSide];
      event.targets[sourceSide][b] = potential;
      let passage = direction < 0 ? b : b + 1;
      for (;;) {
        event.passagePotentials[passage] = potential;
        passages[passage].push({ id: event.id, body: b, edge, potential });
        if (direction < 0 && passage === 0 || direction > 0 && passage === count) {
          const boundary = direction < 0 ? 0 : 1, interval = outerData[boundary];
          checkInterior(potential, interval.inlet, interval.outlet, 'outer boundary ' + boundary);
          event.targets.outer[boundary] = potential;
          break;
        }
        const receiving = direction < 0 ? passage - 1 : passage;
        const side = direction < 0 ? 'upper' : 'lower', opposite = direction < 0 ? 'lower' : 'upper';
        const receiver = bodyData[receiving], trailing = receiver.trailing[side];
        checkInterior(potential, receiver.inlet, receiver.outlets[side], 'body ' + receiving + ' ' + side);
        if (potential === receiver.leading || potential === trailing)
          throw new Error('Coincident LE/TE potential block targets require an explicit shared block.');
        event.targets[side][receiving] = potential;
        if (potential > receiver.leading && potential < trailing) break;
        // Upstream of LE, both faces share the dividing cut potential.
        // Downstream of TE, preserve the wake increment on the other branch.
        if (potential > trailing) {
          potential = receiver.trailing[opposite] + (potential - trailing);
          if (!(potential > receiver.trailing[opposite]))
            throw new Error('Potential block wake increment is not representable.');
        }
        if (!Number.isFinite(potential))
          throw new Error('Nonfinite potential block wake transfer.');
        checkInterior(potential, receiver.inlet, receiver.outlets[opposite], 'body ' + receiving + ' ' + opposite);
        event.targets[opposite][receiving] = potential;
        passage += direction;
      }
    }
    return event;
  }));
  const ordered = orderPotentialBlockEvents(passages), byId = new Map(events.map(event => [event.id, event]));
  const rank = Object.fromEntries(ordered.order.map((event, i) => [event.id, i + 1]));
  for (const event of events) event.rank = rank[event.id];
  const order = ordered.order.map(event => byId.get(event.id)), inletRank = 0, outletRank = events.length + 1;
  const known = (side, boundary, inlet, outlet) => [
    { rank: inletRank, potential: inlet, event: 'inlet' },
    ...order.filter(event => event.targets[side][boundary] !== null)
      .map(event => ({ rank: event.rank, potential: event.targets[side][boundary], event: event.id })),
    { rank: outletRank, potential: outlet, event: 'outlet' },
  ];
  return { events, order, passages: ordered.passages, rank, inletRank, outletRank,
    constraints: {
      lower: bodyData.map((body, b) => known('lower', b, body.inlet, body.outlets.lower)),
      upper: bodyData.map((body, b) => known('upper', b, body.inlet, body.outlets.upper)),
      outer: outerData.map((boundary, b) => known('outer', b, boundary.inlet, boundary.outlet)),
    },
    scope: 'Necessary scalar block constraints only; no receiving geometry, shielded station values, or MSET equivalence is established.' };
}
