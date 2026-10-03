// SPDX-License-Identifier: GPL-2.0-or-later
// Transfer a coupled state to prescribed new surface stations at fixed counts.
// This prepares a guess: every Euler/grid/BL/wake equation must reconverge.
import { createStreamtubeBodySystem } from '../streamtube-body.js';
import { createStreamtubeBoundaryLayers } from '../streamtube-boundary-layers.js';
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { initializeStreamtubeDensities } from '../streamtube-initial-state.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';

const mix = (a, b, t) => t === 0 ? a : t === 1 ? b : a + t * (b - a);
const point = (a, b, t) => ({ x: mix(a.x, b.x, t), y: mix(a.y, b.y, t) });
const sample = (row, coordinate, interpolate) => {
  const i = Math.min(row.length - 2, Math.floor(coordinate));
  return interpolate(row[i], row[i + 1], coordinate - i);
};
function indexAt(fractions, value) {
  let lo = 0, hi = fractions.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (fractions[mid] <= value) lo = mid; else hi = mid;
  }
  if (value === fractions[lo]) return lo;
  if (value === fractions[hi]) return hi;
  return lo + (value - fractions[lo]) / (fractions[hi] - fractions[lo]);
}

export function redistributeCoupledSurfaceStations(input, source, { surfaceFractions, initial = source.initial } = {}) {
  const { layout, conditions } = source.euler;
  if (!Array.isArray(surfaceFractions) || surfaceFractions.length !== layout.elements
    || surfaceFractions.some((body, b) => ['upper', 'lower'].some(side => {
      const f = body?.[side];
      return !Array.isArray(f) || f.length !== source.euler.fractions[b][side].length
        || f[0] !== 0 || f.at(-1) !== 1 || !f.every(Number.isFinite)
        || f.some((v, i) => i && v <= f[i - 1]);
    }))) throw new Error('Surface redistribution requires increasing fractions at unchanged counts and endpoints.');
  if (conditions.flowModel !== 'compressible' || conditions.streamwiseMode !== 'isentropic')
    throw new Error('Surface redistribution currently requires smooth isentropic coupled flow.');
  if (!source.admissible(initial)) throw new Error('Surface redistribution requires an admissible source.');
  const old = source.evaluate(initial), nextInput = structuredClone(input);
  nextInput.bodies.forEach((body, b) => { body.surfaceFractions = structuredClone(surfaceFractions[b]); });
  const options = { reynolds: source.conditions.reynolds, ncrit: source.conditions.ncrit,
    edgeMatching: source.conditions.edgeMatching, tripFractions: structuredClone(source.bl.trips),
    ...(source.bl.transitionMode === 'automatic' ? { transitionMode: 'automatic', transitionState: source.bl.snapshotActive() } : {}) };
  const identity = surfaceFractions.every((body, b) => ['upper', 'lower'].every(side =>
    body[side].every((f, k) => f === source.euler.fractions[b][side][k])));
  if (identity) {
    const initialEuler = { x: initial.slice(0, source.ne), nodes: old.outer.nodes, undisplacedNodes: old.outer.undisplacedNodes };
    const initialBL = initial.slice(source.ne);
    return { input: nextInput, options, initialEuler, initialBL,
      system: createCoupledStreamtubeBody(nextInput, { ...options, initialEuler, initialBL }), diagnostics: { identity: true } };
  }
  nextInput.gridSpacing = { coordinate: 'redistributed surface fractions at fixed station and tube counts',
    surfaceIntervalsByElement: nextInput.bodies.map((b, element) => ({ element, intervals: b.trailingIndex - b.leadingIndex })) };
  const maps = surfaceFractions.map((body, b) => Object.fromEntries(['upper', 'lower'].map(side =>
    [side, body[side].map(f => indexAt(source.euler.fractions[b][side], f))])));
  const zero = { surfaces: nextInput.bodies.map(b => Object.fromEntries(['upper', 'lower'].map(side =>
    [side, Array(b.trailingIndex - b.leadingIndex + 1).fill(0)]))),
    wakes: nextInput.bodies.map(b => Array(layout.nx - b.trailingIndex).fill(0)) };
  const euler = createStreamtubeBodySystem({ ...nextInput, displacement: zero });
  let state = euler.initial.slice(); state.set(initial.subarray(0, layout.densityCount));
  for (const key of Object.keys(layout.globals)) {
    const previous = Array.isArray(layout.globals[key]) ? layout.globals[key] : [layout.globals[key]];
    const next = Array.isArray(euler.layout.globals[key]) ? euler.layout.globals[key] : [euler.layout.globals[key]];
    previous.forEach((col, k) => { if (col !== null) state[next[k]] = initial[col]; });
  }
  for (const key of ['lengthScale', 'massScale']) if (euler.conditions[key] !== conditions[key])
    throw new Error(`Surface redistribution changed ${key}.`);
  const bl = createStreamtubeBoundaryLayers(euler, state, options), geometry = bl.geometry(state);
  if (options.transitionMode === 'automatic') bl.restoreActive(bl.surfaces.map(branch => {
    const location = old.layers.transitions.find(t => t.body === branch.body && t.side === branch.side).s;
    const index = branch.ids.findIndex(id => geometry.coordinates[id].s >= location);
    return index < 0 ? branch.ids.length - 1 : index;
  }));
  const values = new Float64Array(4 * bl.stations.length);
  const interpolate = (a, b, t) => Object.fromEntries(['aux', 'theta', 'deltaStar', 'ue'].map(key => [key, mix(a[key], b[key], t)]));
  const put = (id, v) => values.set([v.aux, v.theta / bl.scale, v.deltaStar / bl.scale, v.ue], 4 * id);
  for (const branch of bl.surfaces) {
    const parent = source.bl.surfaces.find(b => b.body === branch.body && b.side === branch.side);
    const before = parent.ids.map(id => old.layers.states[id]);
    for (const id of branch.ids) {
      const station = bl.stations[id], u = maps[branch.body][branch.side][station.k];
      let v;
      if (u < 1) {
        v = { ...before[0], aux: 0, ue: before[0].ue * geometry.coordinates[id].s / before[0].s };
        if (options.transitionMode === 'automatic' && ['leading-transition', 'transition', 'turbulent'].includes(station.regime))
          v.aux = bl.kernel.station({ ...v, s: geometry.coordinates[id].s, aux: .03 }, 'turbulent').transitionShear;
      }
      else {
        v = sample(before, u - 1, interpolate);
        if (!Number.isInteger(u)) {
          const left = Math.floor(u) - 1, right = left + 1;
          const turbulent = ['leading-transition', 'transition', 'turbulent'].includes(station.regime);
          const wasTurbulent = k => ['leading-transition', 'transition', 'turbulent'].includes(source.bl.stations[parent.ids[k]].regime);
          if (wasTurbulent(left) !== wasTurbulent(right)) v.aux = before[wasTurbulent(left) === turbulent ? left : right].aux;
        }
      }
      put(id, v);
    }
  }
  for (const wake of bl.wakes) {
    const parent = source.bl.wakes.find(w => w.body === wake.body);
    wake.ids.forEach((id, k) => put(id, old.layers.states[parent.ids[k]]));
  }
  if (options.transitionMode === 'automatic') {
    // The identity transfer returned above. On new stations, initialize N
    // consistently with the newly selected interval, including unchanged indices.
    bl.updateActive(values, state, { reinitializeAmplification: true }); options.transitionState = bl.snapshotActive();
  }
  if (!bl.admissible(values)) throw new Error('Redistributed boundary-layer guess is outside the model domain.');
  euler.setDisplacement(bl.thicknesses(values));
  const boundaries = euler.decode(state).nodes;
  const atBank = (body, side, i) => body >= 0 && body < layout.elements && layout.active(body, i)
    ? layout.bodies[body].leadingIndex + maps[body][side][i - layout.bodies[body].leadingIndex] : i;
  let maxNodeMovement = 0;
  const nodes = old.outer.nodes.map((grid, g) => {
    const masses = old.outer.allocation.groups[g].map(t => t.massFlow), total = masses.reduce((a, b) => a + b, 0), eta = [0];
    for (const m of masses) eta.push(eta.at(-1) + m / total); eta[masses.length] = 1;
    return grid.map((row, i) => {
      const lower = atBank(g - 1, 'upper', i), upper = atBank(g, 'lower', i);
      // A convex blend of increasing bank index maps preserves ordering on
      // every interior streamline. It is a local reparameterization of the
      // existing physical grid, not full-path normalized arc transfer.
      const sampled = row.map((_, j) => sample(grid, mix(lower, upper, eta[j]), (a, b, t) => point(a[j], b[j], t)));
      const low = g > 0 && layout.active(g - 1, i) ? boundaries[g][i][0] : sampled[0];
      const high = g < layout.elements && layout.active(g, i) ? boundaries[g][i].at(-1) : sampled.at(-1);
      const dl = { x: low.x - sampled[0].x, y: low.y - sampled[0].y };
      const du = { x: high.x - sampled.at(-1).x, y: high.y - sampled.at(-1).y };
      return sampled.map((p, j) => {
        const q = j === 0 ? { ...low } : j === masses.length ? { ...high }
          : { x: p.x + (1 - eta[j]) * dl.x + eta[j] * du.x, y: p.y + (1 - eta[j]) * dl.y + eta[j] * du.y };
        maxNodeMovement = Math.max(maxNodeMovement, Math.hypot(q.x - row[j].x, q.y - row[j].y));
        return q;
      });
    });
  });
  state = euler.adoptGeometry(state, nodes);
  const flow = euler.decode(state), mesh = streamtubeMeshSnapshot({ system: euler, nodes: flow.nodes });
  if (!mesh.quality.valid) throw new Error('Redistributed physical grid is not positive.', { cause: mesh.quality });
  state = initializeStreamtubeDensities(euler, state);
  const initialEuler = { x: state, nodes: flow.nodes, undisplacedNodes: flow.undisplacedNodes }, initialBL = values;
  const system = createCoupledStreamtubeBody(nextInput, { ...options, initialEuler, initialBL });
  system.evaluate(system.initial);
  if (!system.admissible(system.initial)) throw new Error('Redistributed coupled guess is inadmissible.');
  return { input: nextInput, options, initialEuler, initialBL, system,
    diagnostics: { identity: false, unknowns: system.n, maxNodeMovement, quality: mesh.quality, sourceIndexMaps: maps } };
}
