// SPDX-License-Identifier: GPL-2.0-or-later
// Research domain control, not an MSET algorithm or a production default.
// Add only exterior tubes. Retain physical inner cells and all BL/wake data
// in the starting guess; the changed boundary requires a complete new solve.
import assert from 'node:assert/strict';
import { createStreamtubeBodySystem } from '../../src/euler/streamtube-body.js';
import { createCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';
import { initializeStreamtubeDensities } from '../../src/euler/streamtube-initial-state.js';

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const create = f => createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });

export function extendCoupledLaterally(parentCheckpoint, { lowerWidths, upperWidths, maxNodes = 50000 } = {}) {
  assert.equal(parentCheckpoint.version, 1);
  for (const widths of [lowerWidths, upperWidths]) assert.ok(Array.isArray(widths)
    && widths.every(w => Number.isFinite(w) && w > 0), 'Exterior widths must be positive physical lengths.');
  assert.ok(lowerWidths.length + upperWidths.length > 0, 'Add at least one exterior tube.');
  assert.ok(Number.isInteger(maxNodes) && maxNodes > 0, 'Invalid node budget.');
  const f = parentCheckpoint.restart, parent = create(f), old = parent.euler.layout;
  assert.equal(old.independentWakeBanks, true, 'This control requires independent physical wake banks.');
  assert.equal(parent.euler.conditions.flowModel, 'compressible');
  assert.equal(parent.euler.conditions.streamwiseMode, 'isentropic', 'New densities use the subsonic uniform-entropy branch.');
  assert.equal(parent.euler.conditions.stagnationMotion, 'walls-only');
  assert.ok((old.nx + 1) * (old.tubes.reduce((sum, n) => sum + n + 1, 0)
    + lowerWidths.length + upperWidths.length) <= maxNodes, 'Exterior extension exceeds node budget.');
  const before = parent.evaluate(parent.initial);
  assert.deepEqual(before.families, parentCheckpoint.families);
  assert.ok(parent.admissible(parent.initial), 'Parent grid and BL must be admissible.');
  const input = structuredClone(f.input), last = old.elements;
  const sum = row => row.reduce((a, b) => a + b, 0);
  const lower = sum(lowerWidths), upper = sum(upperWidths);
  const radians = parent.euler.conditions.alpha * Math.PI / 180;
  const normal = { x: -Math.sin(radians), y: Math.cos(radians) };
  const shift = (p, distance) => ({ x: p.x + distance * normal.x, y: p.y + distance * normal.y });
  // Widths are ordered from the old outer boundary outwards. Freestream
  // rho*U is one, so each added width is also its prescribed captured mass.
  const distances = widths => { let total = 0; return widths.map(w => total += w); };
  const dl = distances(lowerWidths), du = distances(upperWidths);
  const physicalNodes = before.outer.nodes.map((grid, g) => grid.map(row => [
    ...(g === 0 ? dl.map(d => shift(row[0], -d)).reverse() : []),
    ...row.map(p => ({ ...p })),
    ...(g === last ? du.map(d => shift(row.at(-1), d)) : []),
  ]));
  if (lower > 0) input.outerLower = physicalNodes[0].map(row => ({ ...row[0] }));
  if (upper > 0) input.outerUpper = physicalNodes[last].map(row => ({ ...row.at(-1) }));
  // Rebase nonprimary capture levels to the parent's *solved* values.
  // Copying its normalized capture unknowns would multiply them by the new
  // total mass and would silently change every passage's initial tube masses.
  input.captureLevels = before.outer.captured.slice();
  input.captureLevels[0] -= lower; input.captureLevels[last + 1] += upper;
  input.weights = before.outer.allocation.groups.map((group, g) => [
    ...(g === 0 ? lowerWidths.slice().reverse() : []),
    ...group.map(tube => tube.massFlow),
    ...(g === last ? upperWidths : []),
  ]);
  input.gridSpacing = { coordinate: 'retained physical inner grid with added exterior tubes',
    parentTubes: old.tubes.slice(), lowerWidths: lowerWidths.slice(), upperWidths: upperWidths.slice(),
    retainedSurfaceStations: true, retainedInnerMasses: true };
  const target = createStreamtubeBodySystem({ ...input, displacement: parent.bl.thicknesses(parent.initial.subarray(parent.ne)) });
  const layout = target.layout;
  let x = target.initial.slice();
  assert.equal(target.conditions.lengthScale, parent.euler.conditions.lengthScale);
  for (const [name, columns] of Object.entries(old.globals)) {
    if (name === 'capture') continue;
    const a = Array.isArray(columns) ? columns : [columns], cols = layout.globals[name];
    const b = Array.isArray(cols) ? cols : [cols];
    a.forEach((col, k) => { if (col !== null) x[b[k]] = parent.initial[col]; });
  }
  x = target.adoptGeometry(x, physicalNodes);
  x = initializeStreamtubeDensities(target, x);
  const offset = g => g === 0 ? lowerWidths.length : 0;
  // Preserve old log densities exactly. Gas inversion initializes only the
  // added tubes in the returned state, even though its geometry guard checks
  // the whole grid. All retained density/mass/section velocities are audited.
  for (let i = 0; i < old.nx; i++) for (let g = 0; g <= last; g++) for (let j = 0; j < old.tubes[g]; j++)
    x[layout.densityIndex(i, g, j + offset(g))] = parent.initial[old.densityIndex(i, g, j)];
  const flow = target.decode(x), options = structuredClone(f.options);
  if (parent.bl.transitionMode === 'automatic') options.transitionState = parent.bl.snapshotActive();
  const restart = serialize({ input, options, initialEuler: { x, nodes: flow.nodes, undisplacedNodes: flow.undisplacedNodes },
    initialBL: parent.initial.slice(parent.ne) });
  const child = create(restart), after = child.evaluate(child.initial);
  assert.ok(child.admissible(child.initial), 'Extended physical grid and BL must be admissible.');
  assert.deepEqual(child.bl.stations, parent.bl.stations);
  assert.deepEqual(child.bl.snapshotActive(), parent.bl.snapshotActive());
  assert.deepEqual(child.bl.trips, parent.bl.trips);
  assert.deepEqual(child.euler.fractions, parent.euler.fractions);
  assert.deepEqual(child.conditions, parent.conditions);
  assert.deepEqual(child.initial.slice(child.ne), parent.initial.slice(parent.ne));
  let nodeError = 0, massRelativeError = 0, physicalBLError = 0, sectionVelocityError = 0;
  before.outer.nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    const q = after.outer.nodes[g][i][j + offset(g)];
    nodeError = Math.max(nodeError, Math.hypot(p.x - q.x, p.y - q.y));
  })));
  before.outer.allocation.groups.forEach((group, g) => group.forEach((t, j) => {
    massRelativeError = Math.max(massRelativeError, Math.abs(after.outer.allocation.groups[g][j + offset(g)].massFlow / t.massFlow - 1));
  }));
  for (const station of parent.bl.stations) for (const key of ['s', 'aux', 'theta', 'deltaStar', 'ue'])
    physicalBLError = Math.max(physicalBLError, Math.abs(before.layers.states[station.id][key] - after.layers.states[station.id][key]));
  before.outer.sections.forEach((row, i) => row.forEach((group, g) => group.forEach((section, j) => {
    const other = after.outer.sections[i][g][j + offset(g)];
    sectionVelocityError = Math.max(sectionVelocityError, Math.abs(section.q - other.q));
  })));
  assert.ok(nodeError < 2e-12 * parent.euler.conditions.lengthScale);
  assert.ok(massRelativeError < 1e-13 && physicalBLError < 1e-14 && sectionVelocityError < 2e-11);
  assert.deepEqual(after.outer.strengths, before.outer.strengths);
  assert.deepEqual(after.outer.stagnation, before.outer.stagnation);
  // The inlet history is streamwise distance along each dividing cut, not
  // fractions across an outer passage. Those cuts and all their indices are
  // retained. Preserve the actual previous SMOVE trigger history as well.
  const checkpoint = { version: 1, families: after.families, restart,
    continuation: serialize(parentCheckpoint.continuation) };
  const diagnostics = { parentUnknowns: parent.n, unknowns: child.n, parentTubes: old.tubes, tubes: layout.tubes,
    lowerWidths, upperWidths, normal, parentMassScale: parent.euler.conditions.massScale, massScale: child.euler.conditions.massScale,
    nodeError, massRelativeError, physicalBLError, sectionVelocityError,
    parentFamilies: before.families, families: after.families, maintenanceHistoryPreserved: true, initialSMOVERepeated: false };
  return { checkpoint, diagnostics, system: child };
}
