// SPDX-License-Identifier: GPL-2.0-or-later
// Coordinate response to a recorded post-step displacement correction. The
// raw BL thickness is used only to decode geometry: it is never evaluated,
// adopted, or accepted. Ordinary candidate/domain/final-grid gates remain
// the caller's responsibility. No independent wake-gap repair is invented.
import { extendWarmBoundaryIncrements } from './streamtube-displacement.js';

export function respondToCoupledProjectionGeometry(system, proposedPacked, projection, decodedProjected) {
  const require = (ok, message) => { if (!ok) throw new Error(message); };
  const changes = projection?.displacementChanges ?? [];
  require(Array.isArray(changes) && decodedProjected && Array.isArray(decodedProjected.nodes),
    'Projection geometry requires decoded nodes and recorded displacement changes.');
  const diagnostics = { method: 'projection-boundary-increment', active: false, displacementCorrections: changes.length,
    surfaceCorrections: 0, surfaceTrailingEdgeCorrections: 0, wakeCorrections: 0, wakeStartCorrections: 0,
    rawToProjectedBoundaryNodes: 0, changedInteriorNodes: 0, maxDisplacement: 0,
    geometryDecodes: 0, boundaryExtensions: 0, rawBLAccepted: false, equationsChanged: false };
  // This branch intentionally does not inspect the system or packed state,
  // compute thicknesses, decode geometry, or allocate a new node array.
  if (changes.length === 0) return { nodes: decodedProjected.nodes, diagnostics };
  const { ne, n, bl, euler } = system;
  require(Number.isInteger(ne) && ne > 0 && Number.isInteger(n) && n > ne
    && (Array.isArray(proposedPacked) || ArrayBuffer.isView(proposedPacked)) && proposedPacked.length === n
    && Array.isArray(bl?.stations) && n - ne === 4 * bl.stations.length,
  'Projection geometry requires a complete coupled proposed state.');
  const projectedBL = proposedPacked.slice(ne), rawBL = projectedBL.slice(), seen = new Set();
  for (const change of changes) {
    const id = change?.id, station = bl.stations[id];
    require(Number.isInteger(id) && id >= 0 && station?.id === id && !seen.has(id),
      'Projection geometry correction has an invalid or duplicate station.');
    require(['surface', 'wake'].includes(station.kind)
      && ['kind', 'body', 'side', 'i'].every(key => change[key] === station[key])
      && Number.isFinite(change.beforeDeltaStar) && change.beforeDeltaStar > 0
      && Number.isFinite(change.deltaStar) && change.deltaStar > change.beforeDeltaStar
      && projectedBL[4 * id + 2] === change.deltaStar,
    'Projection geometry correction does not match the proposed station displacement.');
    seen.add(id); rawBL[4 * id + 2] = change.beforeDeltaStar;
    if (station.kind === 'surface') {
      diagnostics.surfaceCorrections++;
      if (station.i === euler.layout.bodies[station.body]?.trailingIndex) diagnostics.surfaceTrailingEdgeCorrections++;
    } else {
      diagnostics.wakeCorrections++;
      if (station.regime === 'trailing-edge') diagnostics.wakeStartCorrections++;
    }
  }
  const projectedThicknesses = bl.thicknesses(projectedBL);
  let raw;
  try {
    euler.setDisplacement(bl.thicknesses(rawBL));
    diagnostics.geometryDecodes++;
    raw = euler.decode(proposedPacked.slice(0, ne));
  } finally {
    // Restore even when setting/decoding raw geometry fails. The caller has
    // already installed the projected displacement at this seam.
    euler.setDisplacement(projectedThicknesses);
  }
  const nodes = decodedProjected.nodes;
  require(Array.isArray(raw?.nodes) && raw.nodes.length === nodes.length
    && Array.isArray(decodedProjected.allocation?.groups) && decodedProjected.allocation.groups.length === nodes.length,
  'Projection geometry requires matching raw/projected passages and proposed physical masses.');
  for (let g = 0; g < nodes.length; g++) {
    require(Array.isArray(raw.nodes[g]) && raw.nodes[g].length === nodes[g].length,
      'Projection geometry changed passage station counts.');
    for (let i = 0; i < nodes[g].length; i++) {
      const a = raw.nodes[g][i], b = nodes[g][i];
      require(Array.isArray(a) && Array.isArray(b) && a.length === b.length && a.length >= 2,
        'Projection geometry changed passage tube counts.');
      for (let j = 0; j < a.length; j++) {
        require([a[j]?.x, a[j]?.y, b[j]?.x, b[j]?.y].every(Number.isFinite),
          'Projection geometry contains a nonfinite point.');
        const same = a[j].x === b[j].x && a[j].y === b[j].y;
        // Raw thickness restoration can move only the physical boundaries;
        // a changed free interior would imply a different Euler/chart state.
        if (j !== 0 && j !== a.length - 1) require(same, 'Raw projection geometry changed a free interior node.');
        else if (!same) diagnostics.rawToProjectedBoundaryNodes++;
      }
    }
  }
  diagnostics.independentWakeBanks = euler.layout.independentWakeBanks === true;
  // In independent-bank mode a wake delta changes a residual target, not
  // the bank coordinates. Count it, but do not fabricate a boundary shift.
  if (diagnostics.rawToProjectedBoundaryNodes === 0) return { nodes, diagnostics };
  const masses = decodedProjected.allocation.groups.map(group => group.map(tube => tube.massFlow));
  diagnostics.boundaryExtensions++;
  const extended = extendWarmBoundaryIncrements({ sourceNodes: raw.nodes, targetNodes: nodes, masses });
  for (let g = 0; g < nodes.length; g++) for (let i = 0; i < nodes[g].length; i++)
    for (let j = 1; j < nodes[g][i].length - 1; j++) {
      const a = nodes[g][i][j], b = extended[g][i][j], distance = Math.hypot(b.x - a.x, b.y - a.y);
      if (distance > 0) diagnostics.changedInteriorNodes++;
      diagnostics.maxDisplacement = Math.max(diagnostics.maxDisplacement, distance);
    }
  diagnostics.active = diagnostics.changedInteriorNodes > 0;
  return { nodes: extended, diagnostics };
}
