// SPDX-License-Identifier: GPL-2.0-or-later

// Fixed material TE data. TECALC's projection algebra is applied to the
// actual solid C2 curve used by this Euler model. Its natural endpoint spline
// differs from the SPLIND endpoint convention of XFOIL's panel initializer.
export function streamtubeBaseGeometry(bodies, curves) {
  return bodies.map((body, b) => {
    if (body.trailingEdge?.kind !== 'finite-base') return null;
    const upper = curves[b].evaluate(0), lower = curves[b].evaluate(curves[b].length);
    const vector = { x: upper.point.x - lower.point.x, y: upper.point.y - lower.point.y };
    const tangentDerivative = { x: .5 * (-upper.derivative.x + lower.derivative.x),
      y: .5 * (-upper.derivative.y + lower.derivative.y) };
    const width = tangentDerivative.x * vector.y - tangentDerivative.y * vector.x;
    const tangentialProjection = tangentDerivative.x * vector.x + tangentDerivative.y * vector.y;
    if (!(width > 0) || !Number.isFinite(width) || !Number.isFinite(tangentialProjection))
      throw new Error('Finite-base geometry needs positive TECALC normal projection on its material surface tangents.');
    return { width, tangentialProjection, vector, tangentDerivative,
      upperDerivative: { ...upper.derivative }, lowerDerivative: { ...lower.derivative },
      center: { x: .5 * (upper.point.x + lower.point.x), y: .5 * (upper.point.y + lower.point.y) },
      convention: 'TECALC projection algebra on the Euler natural-C2 material curve; fixed in the flow solve.' };
  });
}

// MISES user guide 2.63, section 5.1.7: a fresh inviscid finite-TE wake has
// constant width and moves to impose zero pressure jump. This is not the
// shrinking prescribed dead-air component of a viscous wake.
export function initialStreamtubeDisplacement({ bodies, nx }, baseGeometry = bodies.map(() => null)) {
  return {
    surfaces: bodies.map(body => ({ upper: Array(body.trailingIndex - body.leadingIndex + 1).fill(0),
      lower: Array(body.trailingIndex - body.leadingIndex + 1).fill(0) })),
    wakes: bodies.map((body, b) => Array(nx - body.trailingIndex).fill(baseGeometry[b]?.width ?? 0)),
  };
}

// Normal-coordinate directions, compared with Giles NCALC pp.196–197.
// In its single-blade listing, the one-sided LE/TE neighbor selection
// persists across the transverse loop. The optional multielement extension
// applies this to each passage adjacent to that body. Other shared cuts
// retain their own centered direction; a cut has only one movement unknown.
// Coordinates are the physical flow grid. With a displaced wake, the shared
// unknown is its centerline, so its secant uses the mean of the two banks.
export function streamtubeMotionDirections(layout, nodes, mode = 'centered') {
  if (!['centered', 'body-stations'].includes(mode)) throw new Error('Unknown streamtube normal stencil.');
  const { nx, elements, tubes, bodies } = layout, directions = new Map();
  for (let g = 0; g <= elements; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++) {
    const node = layout.nodes[g][i][j], col = node.column;
    if (col === null || directions.has(col)) continue;
    let left = Math.max(0, i - 1), right = Math.min(nx, i + 1);
    if (mode === 'body-stations' && node.kind !== 'cut') {
      const banks = [bodies[g - 1], bodies[g]].filter(Boolean);
      const leading = banks.some(b => b.leadingIndex === i), trailing = banks.some(b => b.trailingIndex === i);
      if (leading && trailing) throw new Error('Conflicting LE/TE motion stencils in one passage.');
      if (leading) right = i;
      if (trailing) left = i;
    }
    const point = i => {
      if (node.kind !== 'cut' || (layout.independentWakeBanks && node.side)) return nodes[g][i][j];
      const a = nodes[node.body][i].at(-1), b = nodes[node.body + 1][i][0];
      return { x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) };
    };
    const a = point(left), b = point(right), dx = b.x - a.x, dy = b.y - a.y;
    const length = Math.hypot(dx, dy);
    if (!(length > 0) || !Number.isFinite(length)) throw new Error('Degenerate body grid movement direction.');
    directions.set(col, { x: -dy / length, y: dx / length });
  }
  return directions;
}

// Change only the wake-coordinate parameterization. The receiving system
// must have the same contours, stations, mass weights and displacement data.
export function transferStreamtubeGeometry(source, state, target) {
  const a = source.layout, b = target.layout;
  if (a.nx !== b.nx || a.elements !== b.elements || a.densityCount !== b.densityCount
    || a.tubes.some((n, i) => n !== b.tubes[i])
    || a.bodies.some((p, i) => p.leadingIndex !== b.bodies[i].leadingIndex || p.trailingIndex !== b.bodies[i].trailingIndex
      || (p.trailingEdge?.kind === 'finite-base' || b.bodies[i].trailingEdge?.kind === 'finite-base')
        && ['kind', 'upperIndex', 'lowerIndex'].some(k => p.trailingEdge?.[k] !== b.bodies[i].trailingEdge?.[k])))
    throw new Error('Wake geometry transfer requires identical stations and passages.');
  const x = target.initial.slice(); x.set(state.subarray(0, a.densityCount));
  for (const key of Object.keys(a.globals)) {
    const old = Array.isArray(a.globals[key]) ? a.globals[key] : [a.globals[key]];
    const next = Array.isArray(b.globals[key]) ? b.globals[key] : [b.globals[key]];
    old.forEach((col, i) => {
      if ((col === null) !== (next[i] === null)) throw new Error('Incompatible global wake-transfer unknowns.');
      if (col !== null) x[next[i]] = state[col];
    });
  }
  return target.adoptGeometry(x, source.decode(state).nodes);
}
