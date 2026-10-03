// SPDX-License-Identifier: GPL-2.0-or-later
// Move the downstream boundary to an existing cross-line. Retained physical
// nodes, section densities, BL states, masses and trips are copied exactly.
// The changed farfield/outlet equations require a new complete coupled solve.
import { createStreamtubeBodySystem } from '../streamtube-body.js';
import { createStreamtubeBoundaryLayers } from '../streamtube-boundary-layers.js';
import { createCoupledStreamtubeBody } from '../streamtube-coupled.js';
import { streamtubeMeshSnapshot } from '../streamtube-mesh-preview.js';

export function truncateCoupledDownstreamDomain(input, source, { endIndex, initial = source.initial } = {}) {
  const oldLayout = source.euler.layout;
  if (!oldLayout.independentWakeBanks)
    throw new Error('Exact downstream truncation requires independent physical wake banks.');
  if (!Number.isInteger(endIndex) || endIndex > oldLayout.nx
    || oldLayout.bodies.some(b => endIndex < b.trailingIndex + 2))
    throw new Error('Keep at least two wake intervals after every trailing edge when truncating the domain.');
  if (!source.admissible(initial)) throw new Error('Domain truncation requires an admissible parent.');
  const previous = source.evaluate(initial), nextInput = structuredClone(input);
  nextInput.outerLower = input.outerLower.slice(0, endIndex + 1).map(p => ({ ...p }));
  nextInput.outerUpper = input.outerUpper.slice(0, endIndex + 1).map(p => ({ ...p }));
  nextInput.cutPaths = input.cutPaths.map(path => path.slice(0, endIndex + 1).map(p => ({ ...p })));
  if (endIndex !== oldLayout.nx) nextInput.gridSpacing = { coordinate: 'retained cross-lines with a nearer downstream boundary',
    parentSegments: oldLayout.nx, segments: endIndex, retainedSurfaceStations: true };
  const displacement = source.bl.thicknesses(initial.slice(source.ne));
  displacement.wakes = displacement.wakes.map((row, b) => row.slice(0, endIndex - oldLayout.bodies[b].trailingIndex));
  const target = createStreamtubeBodySystem({ ...nextInput, displacement }), layout = target.layout;
  let x = target.initial.slice();
  x.set(initial.subarray(0, layout.densityCount));
  for (const [key, columns] of Object.entries(oldLayout.globals)) {
    const old = Array.isArray(columns) ? columns : [columns];
    const next = Array.isArray(layout.globals[key]) ? layout.globals[key] : [layout.globals[key]];
    old.forEach((col, k) => { if (col !== null) x[next[k]] = initial[col]; });
  }
  for (const key of ['lengthScale', 'massScale']) if (target.conditions[key] !== source.euler.conditions[key])
    throw new Error(`Domain truncation changed ${key}.`);
  const nodes = previous.outer.nodes.map(grid => grid.slice(0, endIndex + 1));
  x = target.adoptGeometry(x, nodes);
  const flow = target.decode(x), mesh = streamtubeMeshSnapshot({ system: target, nodes: flow.nodes });
  if (!mesh.quality.valid) throw new Error('Truncated physical grid is not positive.', { cause: mesh.quality });
  const options = { reynolds: source.conditions.reynolds, ncrit: source.conditions.ncrit,
    edgeMatching: source.conditions.edgeMatching, tripFractions: structuredClone(source.bl.trips),
    ...(source.bl.transitionMode === 'automatic' ? { transitionMode: 'automatic', transitionState: source.bl.snapshotActive() } : {}) };
  const bl = createStreamtubeBoundaryLayers(target, x, options), initialBL = new Float64Array(4 * bl.stations.length);
  const identity = s => `${s.kind}/${s.body}/${s.side ?? ''}/${s.i}`;
  const oldStations = new Map(source.bl.stations.map(s => [identity(s), s.id]));
  for (const station of bl.stations) {
    const old = oldStations.get(identity(station));
    if (old === undefined) throw new Error('Truncated BL station is absent in the parent.');
    initialBL.set(initial.subarray(source.ne + 4 * old, source.ne + 4 * old + 4), 4 * station.id);
  }
  const initialEuler = { x, nodes: flow.nodes, undisplacedNodes: flow.undisplacedNodes };
  const system = createCoupledStreamtubeBody(nextInput, { ...options, initialEuler, initialBL });
  const value = system.evaluate(system.initial);
  let nodeError = 0, blError = 0;
  nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    const q = value.outer.nodes[g][i][j]; nodeError = Math.max(nodeError, Math.hypot(p.x - q.x, p.y - q.y));
  })));
  for (const station of bl.stations) {
    const a = value.layers.states[station.id], b = previous.layers.states[oldStations.get(identity(station))];
    for (const key of ['s', 'aux', 'theta', 'deltaStar', 'ue']) blError = Math.max(blError, Math.abs(a[key] - b[key]));
  }
  if (nodeError > 2e-12 || blError > 1e-14) throw new Error('Domain truncation changed retained physical data.');
  return { input: nextInput, options, initialEuler, initialBL, system,
    diagnostics: { parentSegments: oldLayout.nx, segments: layout.nx, parentUnknowns: source.n, unknowns: system.n,
      nodeError, blError, quality: mesh.quality, outletWakeCenters: input.bodies.map((_, b) => {
        const lower = value.outer.nodes[b][endIndex].at(-1), upper = value.outer.nodes[b + 1][endIndex][0];
        return { x: .5 * (lower.x + upper.x), y: .5 * (lower.y + upper.y) };
      }) } };
}
