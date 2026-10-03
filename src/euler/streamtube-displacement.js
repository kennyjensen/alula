// SPDX-License-Identifier: GPL-2.0-or-later
// Displacement-surface geometry, Drela (1986), section 6.8, eqs. 6.57–59.
// Thicknesses here are prescribed physical lengths. A coupled BL assembly
// must eventually supply their unknowns and matching equations.
const zero = { x: 0, y: 0 };
const sub = (a, b) => ({ x: a.x - b.x, y: a.y - b.y });
const mean = (a, b) => ({ x: .5 * (a.x + b.x), y: .5 * (a.y + b.y) });
const offset = (p, n, d) => ({ x: p.x + d * n.x, y: p.y + d * n.y });
const add = (row, col, p, weight = 1) => {
  const a = row.get(col) ?? zero; row.set(col, { x: a.x + weight * p.x, y: a.y + weight * p.y });
};
const combine = (a, b, wa = .5, wb = .5) => {
  const result = new Map(); for (const [col, p] of a) add(result, col, p, wa);
  for (const [col, p] of b) add(result, col, p, wb); return result;
};
const normal = (v, sign = 1) => {
  const length = Math.hypot(v.x, v.y);
  if (!(length > 0)) throw new Error('Degenerate displacement-surface tangent.');
  const t = { x: v.x / length, y: v.y / length };
  return { value: { x: -sign * t.y, y: sign * t.x }, apply: dv => {
    const along = t.x * dv.x + t.y * dv.y;
    return { x: -sign * (dv.y - t.y * along) / length, y: sign * (dv.x - t.x * along) / length };
  } };
};

export function createStreamtubeDisplacement({ layout, curves, fractions, thicknesses }) {
  const { bodies, nx, tubes, elements, globals } = layout;
  const values = structuredClone(thicknesses);
  // Optional declared tangential bank pairing. B is the FIXED solid TE
  // vector, not the variable displaced TE difference. The wake vector is
  // D*n + (I-n*n^T)*B: retain the native normal gap without discarding a
  // finite base's material tangential separation. Independent-bank mode
  // uses this only when constructing its initial centerline-based seed.
  const materialWakeVectors = bodies.map((body, b) => {
    if (body.wakeTangentialReference !== undefined && body.wakeTangentialReference !== 'material-te')
      throw new Error('Unknown tangential wake reference.');
    return body.wakeTangentialReference === 'material-te'
      ? sub(curves[b].evaluate(0).point, curves[b].evaluate(curves[b].length).point) : zero;
  });
  if (!values || !Array.isArray(values.surfaces) || !Array.isArray(values.wakes)
    || values.surfaces.length !== elements || values.wakes.length !== elements)
    throw new Error('Supply upper/lower displacement thickness and a total wake gap for every body.');
  const valid = (a, count) => Array.isArray(a) && a.length === count && a.every(d => Number.isFinite(d) && d >= 0);
  for (let b = 0; b < elements; b++) {
    const count = bodies[b].trailingIndex - bodies[b].leadingIndex + 1, s = values.surfaces[b];
    if (!s || !valid(s.upper, count) || !valid(s.lower, count) || !valid(values.wakes[b], nx - bodies[b].trailingIndex))
      throw new Error('Invalid displacement thickness dimensions or negative/nonfinite thickness.');
    if (s.upper[0] !== s.lower[0]) throw new Error('Both surface limits must share the same leading-edge displacement.');
  }
  const parameters = [], surfaceColumns = [], wakeColumns = [];
  const parameter = item => { const column = layout.n + parameters.length; parameters.push(item); return column; };
  for (let b = 0; b < elements; b++) {
    const leading = parameter({ kind: 'surface', body: b, side: 'both', index: 0, value: values.surfaces[b].lower[0] });
    surfaceColumns[b] = Object.fromEntries(['lower', 'upper'].map(side => [side, values.surfaces[b][side].map((value, index) =>
      index === 0 ? leading : parameter({ kind: 'surface', body: b, side, index, value }))]));
    wakeColumns[b] = values.wakes[b].map((value, index) => parameter({ kind: 'wake', body: b, index, value }));
  }
  const wallNormal = (b, side, k, stagnation) => {
    const f = fractions[b][side][k], v = curves[b].branch(side, f, stagnation[b]);
    // Curves run CCW from upper TE to lower TE: outward is the right normal.
    const n = normal(v.derivative, -1);
    const d = n.apply(v.secondDerivative), scale = (1 - f) * curves[b].length;
    return { value: n.value, derivative: { x: scale * d.x, y: scale * d.y } };
  };
  // An absolute zero-thickness coordinate chart, not a prescribed wake
  // shape or gap. Both independent banks translate with the displaced TE
  // center while retaining their own free normal position unknowns.
  const wakeTranslation = (b, stagnation) => {
    const k = bodies[b].trailingIndex - bodies[b].leadingIndex;
    const normals = Object.fromEntries(['lower', 'upper'].map(side => [side, wallNormal(b, side, k, stagnation)]));
    const shift = {}, derivative = {};
    for (const key of ['x', 'y']) {
      shift[key] = .5 * (values.surfaces[b].lower[k] * normals.lower.value[key]
        + values.surfaces[b].upper[k] * normals.upper.value[key]);
      derivative[key] = .5 * (values.surfaces[b].lower[k] * normals.lower.derivative[key]
        + values.surfaces[b].upper[k] * normals.upper.derivative[key]);
    }
    return { shift, derivative, normals, k };
  };
  const apply = (nodes, stagnation, derivatives, includeDisplacement = false) => {
    const result = nodes.map(group => group.map(row => row.map(p => ({ ...p }))));
    const maps = derivatives?.map(group => group.map(row => row.map(m => new Map(m))));
    for (let b = 0; b < elements; b++) {
      const body = bodies[b];
      for (const side of ['lower', 'upper']) {
        const g = side === 'lower' ? b : b + 1, j = side === 'lower' ? tubes[b] : 0;
        for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
          const k = i - body.leadingIndex, thickness = values.surfaces[b][side][k];
          const n = wallNormal(b, side, k, stagnation);
          result[g][i][j] = offset(nodes[g][i][j], n.value, thickness);
          if (maps && globals.stagnation[b] !== null) add(maps[g][i][j], globals.stagnation[b], n.derivative, thickness);
          if (maps && includeDisplacement) add(maps[g][i][j], surfaceColumns[b][side][k], n.value);
        }
      }
      // Independent-bank mode solves the normal-gap condition as an Euler
      // row. The optional chart translation carries both banks together;
      // changing a wake gap still never directly repositions either bank.
      if (layout.independentWakeBanks) {
        if (layout.wakeDisplacementMotion === 'te-center') {
          const { shift, derivative, normals, k } = wakeTranslation(b, stagnation);
          for (let i = body.trailingIndex + 1; i <= nx; i++) for (const [g, j] of [[b, tubes[b]], [b + 1, 0]]) {
            result[g][i][j] = { x: nodes[g][i][j].x + shift.x, y: nodes[g][i][j].y + shift.y };
            if (maps && globals.stagnation[b] !== null && (derivative.x !== 0 || derivative.y !== 0))
              add(maps[g][i][j], globals.stagnation[b], derivative);
            if (maps && includeDisplacement) for (const side of ['lower', 'upper'])
              add(maps[g][i][j], surfaceColumns[b][side][k], normals[side].value, .5);
          }
        }
        continue;
      }
      // The wake centerline still has one position unknown at each cut.
      // Its two inviscid banks are distinct; their *total* normal gap is D.
      // The TE center uses both displaced wall endpoints. Downstream center
      // positions come from the underlying shared, moving inviscid cut.
      const te = body.trailingIndex;
      const center = i => i === te ? mean(result[b][i].at(-1), result[b + 1][i][0]) : nodes[b][i].at(-1);
      const centerMap = i => i === te ? combine(maps[b][i].at(-1), maps[b + 1][i][0]) : derivatives[b][i].at(-1);
      for (let i = te + 1; i <= nx; i++) {
        const left = i - 1, right = Math.min(nx, i + 1), n = normal(sub(center(right), center(left)));
        const width = values.wakes[b][i - te - 1], p = center(i);
        result[b][i][tubes[b]] = offset(p, n.value, -.5 * width);
        result[b + 1][i][0] = offset(p, n.value, .5 * width);
        const reference = materialWakeVectors[b], tangent = { x: n.value.y, y: -n.value.x };
        const along = reference.x * tangent.x + reference.y * tangent.y;
        if (along !== 0) {
          result[b][i][tubes[b]] = offset(result[b][i][tubes[b]], tangent, -.5 * along);
          result[b + 1][i][0] = offset(result[b + 1][i][0], tangent, .5 * along);
        }
        if (maps) {
          const tangentMap = combine(centerMap(right), centerMap(left), 1, -1);
          const lower = new Map(centerMap(i)), upper = new Map(centerMap(i));
          for (const [col, dv] of tangentMap) {
            const dn = n.apply(dv); add(lower, col, dn, -.5 * width); add(upper, col, dn, .5 * width);
            if (reference.x !== 0 || reference.y !== 0) {
              const dt = { x: dn.y, y: -dn.x }, da = reference.x * dt.x + reference.y * dt.y;
              const d = { x: da * tangent.x + along * dt.x, y: da * tangent.y + along * dt.y };
              add(lower, col, d, -.5); add(upper, col, d, .5);
            }
          }
          if (includeDisplacement) {
            add(lower, wakeColumns[b][i - te - 1], n.value, -.5);
            add(upper, wakeColumns[b][i - te - 1], n.value, .5);
          }
          maps[b][i][tubes[b]] = lower; maps[b + 1][i][0] = upper;
        }
      }
    }
    return { nodes: result, derivatives: maps };
  };
  const restore = (nodes, stagnation, base, tolerance) => {
    const raw = nodes.map(group => group.map(row => row.map(p => ({ ...p }))));
    for (let b = 0; b < elements; b++) {
      const body = bodies[b];
      for (let i = body.leadingIndex; i <= body.trailingIndex; i++) {
        raw[b][i][tubes[b]] = base[b][i][tubes[b]];
        raw[b + 1][i][0] = base[b + 1][i][0];
      }
      if (!layout.independentWakeBanks) for (let i = body.trailingIndex + 1; i <= nx; i++) {
        const p = mean(nodes[b][i].at(-1), nodes[b + 1][i][0]);
        raw[b][i][tubes[b]] = p; raw[b + 1][i][0] = p;
      }
      else if (layout.wakeDisplacementMotion === 'te-center') {
        const { shift } = wakeTranslation(b, stagnation);
        for (let i = body.trailingIndex + 1; i <= nx; i++) for (const [g, j] of [[b, tubes[b]], [b + 1, 0]])
          raw[g][i][j] = { x: nodes[g][i][j].x - shift.x, y: nodes[g][i][j].y - shift.y };
      }
    }
    const expected = apply(raw, stagnation).nodes;
    for (let g = 0; g <= elements; g++) for (let i = 0; i <= nx; i++) for (let j = 0; j <= tubes[g]; j++)
      if (Math.hypot(expected[g][i][j].x - nodes[g][i][j].x, expected[g][i][j].y - nodes[g][i][j].y) > tolerance)
        throw new Error('Restart wall or wake banks do not match the prescribed displacement geometry.');
    return raw;
  };
  for (const s of values.surfaces) { Object.freeze(s.upper); Object.freeze(s.lower); Object.freeze(s); }
  values.wakes.forEach(Object.freeze); Object.freeze(values.surfaces); Object.freeze(values.wakes); Object.freeze(values);
  parameters.forEach(Object.freeze); Object.freeze(parameters);
  return { values, parameters, apply, restore };
}

// Algebraic starting geometry for subsequent incompressible-flow relaxation.
// This does not establish pressure/gas admissibility or solve an Euler row.

// Seed a normal gap for independent banks. This is only an initialization
// operation, never part of residual evaluation or an accepted Newton step.
export function seedStreamtubeWakeBanks(system, state) {
  const { nodes, undisplacedNodes, stagnation } = system.decode(state);
  const raw = structuredClone(undisplacedNodes), { layout, curves, fractions } = system;
  for (let b = 0; b < layout.elements; b++) for (let i = layout.bodies[b].trailingIndex + 1; i <= layout.nx; i++) {
    // The centerline seed is a physical grid. Preserve the optional chart's
    // TE translation here; the subsequent independent adoption removes it
    // from its stored coordinate bases, so it is never applied twice.
    const centerGrid = layout.wakeDisplacementMotion === 'te-center' ? nodes : raw;
    const a = centerGrid[b][i].at(-1), c = centerGrid[b + 1][i][0], p = { x: .5 * (a.x + c.x), y: .5 * (a.y + c.y) };
    raw[b][i][layout.tubes[b]] = p; raw[b + 1][i][0] = { ...p };
  }
  return createStreamtubeDisplacement({ layout: { ...layout, independentWakeBanks: false }, curves, fractions,
    thicknesses: system.displacement }).apply(raw, stagnation).nodes;
}

export function extendStreamtubeDisplacement(system, state, { inletDisplacement = false } = {}) {
  if (typeof inletDisplacement !== 'boolean') throw new Error('Invalid inlet displacement extension control.');
  if (!system.displacement) throw new Error('Displacement boundaries are required for grid extension.');
  const { nodes, undisplacedNodes: raw, stagnation, allocation } = system.decode(state);
  const moved = structuredClone(raw), { layout, curves, fractions } = system;
  layout.bodies.forEach((body, b) => {
    if (inletDisplacement) {
      // The displaced stagnation point can overtake the last inlet station.
      // Extend its translation along the cut before extending across tubes.
      // Arc distance keeps this continuous on strongly clustered inlets.
      const le = body.leadingIndex, lead = nodes[b][le].at(-1), oldLead = raw[b][le].at(-1), arc = [0];
      for (let i = 1; i <= le; i++) {
        const a = raw[b][i - 1].at(-1), c = raw[b][i].at(-1);
        arc.push(arc.at(-1) + Math.hypot(c.x - a.x, c.y - a.y));
      }
      for (let i = 1; i < le; i++) {
        const a = raw[b][i].at(-1), f = arc[i] / arc[le];
        const p = { x: a.x + f * (lead.x - oldLead.x), y: a.y + f * (lead.y - oldLead.y) };
        moved[b][i][layout.tubes[b]] = p; moved[b + 1][i][0] = { ...p };
      }
    }
    const te = body.trailingIndex;
    const delta = Object.fromEntries(['x', 'y'].map(key => [key, .5 * ((nodes[b][te].at(-1)[key] - raw[b][te].at(-1)[key])
      + (nodes[b + 1][te][0][key] - raw[b + 1][te][0][key]))]));
    // Continue the displaced TE center without introducing a first-wake
    // kink. The outlet center and interior nodes are free Euler unknowns.
    for (let i = te + 1; i <= layout.nx; i++) {
      const a = raw[b][i].at(-1), c = raw[b + 1][i][0];
      const p = { x: .5 * (a.x + c.x) + delta.x, y: .5 * (a.y + c.y) + delta.y };
      // A shared cut can alias the same object through both passages.
      // Assign from immutable coordinates so its shift is applied once.
      moved[b][i][layout.tubes[b]] = p; moved[b + 1][i][0] = { ...p };
    }
  });
  const displaced = createStreamtubeDisplacement({ layout: { ...layout, independentWakeBanks: false }, curves, fractions, thicknesses: system.displacement }).apply(moved, stagnation).nodes;
  return displaced.map((group, g) => {
    const mass = allocation.groups[g].reduce((sum, tube) => sum + tube.massFlow, 0), eta = [0];
    for (const tube of allocation.groups[g]) eta.push(eta.at(-1) + tube.massFlow / mass);
    eta[eta.length - 1] = 1;
    return group.map((row, i) => {
      const n = row.length - 1, delta = [0, n].map(j => ({ x: row[j].x - raw[g][i][j].x, y: row[j].y - raw[g][i][j].y }));
      return row.map((p, j) => j === 0 || j === n ? { ...p } : {
        x: p.x + (1 - eta[j]) * delta[0].x + eta[j] * delta[1].x,
        y: p.y + (1 - eta[j]) * delta[0].y + eta[j] * delta[1].y,
      });
    });
  });
}

// Initial-guess geometry after a BL profile change on an already solved grid.
// Extend boundary increments in physical mass coordinates; never reapply the
// full displacement to the retained interior, or average independent wakes.
// This is not a residual operator and does not certify the resulting grid.
export function extendWarmBoundaryIncrements({ sourceNodes, targetNodes, masses }) {
  const require = (ok, message) => { if (!ok) throw new Error(message); };
  require(Array.isArray(sourceNodes) && sourceNodes.length > 0 && Array.isArray(targetNodes)
    && targetNodes.length === sourceNodes.length && Array.isArray(masses) && masses.length === sourceNodes.length,
  'Warm boundary extension requires matching passage groups and physical masses.');
  const copy = p => ({ ...p });
  return sourceNodes.map((group, g) => {
    const other = targetNodes[g], mass = masses[g];
    require(Array.isArray(group) && group.length >= 2 && Array.isArray(other) && other.length === group.length
      && Array.isArray(mass) && mass.length > 0 && mass.every(m => Number.isFinite(m) && m > 0),
    'Invalid warm passage stations or positive streamtube masses.');
    const total = mass.reduce((a, b) => a + b, 0), eta = [0];
    require(Number.isFinite(total) && total > 0, 'Invalid warm passage total mass.');
    for (const m of mass) eta.push(eta.at(-1) + m / total);
    eta[eta.length - 1] = 1;
    return group.map((row, i) => {
      const target = other[i], last = mass.length;
      require(Array.isArray(row) && row.length === last + 1 && Array.isArray(target) && target.length === row.length
        && [...row, ...target].every(p => p && Number.isFinite(p.x) && Number.isFinite(p.y)),
      'Invalid warm passage points or dimensions.');
      const delta = [0, last].map(j => ({ x: target[j].x - row[j].x, y: target[j].y - row[j].y }));
      const unchanged = delta.every(p => p.x === 0 && p.y === 0);
      return row.map((p, j) => {
        if (j === 0 || j === last) return copy(target[j]);
        if (unchanged) return copy(p);
        return { ...p, x: p.x + (1 - eta[j]) * delta[0].x + eta[j] * delta[1].x,
          y: p.y + (1 - eta[j]) * delta[0].y + eta[j] * delta[1].y };
      });
    });
  });
}
