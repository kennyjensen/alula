// SPDX-License-Identifier: GPL-2.0-or-later
// Continuation metadata for the existing nested coupled refiner. This only
// certifies a transferred guess; it neither solves it nor invokes SMOVE.
import { createCoupledStreamtubeBody } from './streamtube-coupled.js';
import { nestedIntervalMap } from '../geometry/nested-intervals.js';
import { createStreamtubeRefinementCoordinates } from '../geometry/streamtube-refinement-coordinates.js';

// Browser-safe structural checks. Sharp transfers remain identical to
// scripts/validation/nested-refinement-checkpoint.js; tests compare both paths.
// Finite wake transfers additionally distinguish fluid and total displacement.
export function checkpointDataEqual(a, b) {
  if (Object.is(a, b)) return true;
  if (a === null || b === null || typeof a !== 'object' || typeof b !== 'object'
    || Object.getPrototypeOf(a) !== Object.getPrototypeOf(b)) return false;
  if (Array.isArray(a) && a.length !== b.length) return false;
  const keys = Object.keys(a), other = Object.keys(b);
  return keys.length === other.length && keys.every(key => Object.hasOwn(b, key) && checkpointDataEqual(a[key], b[key]));
}
const assert = {
  ok(value, message = 'Invalid nested refinement checkpoint.') { if (!value) throw new Error(message); },
  equal(a, b, message = 'Nested refinement checkpoint values differ.') { this.ok(Object.is(a, b), message); },
  deepEqual(a, b, message = 'Nested refinement checkpoint data differ.') { this.ok(checkpointDataEqual(a, b), message); },
};

const serialize = value => JSON.parse(JSON.stringify(value, (_, v) => ArrayBuffer.isView(v) ? Array.from(v) : v));
const create = f => createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
const mix = (a, b, t) => a + t * (b - a);
const interpolate = (row, counts) => {
  assert.equal(row.length, counts.length + 1);
  const result = [row[0]];
  counts.forEach((n, i) => {
    for (let k = 1; k <= n; k++) result.push(k === n ? row[i + 1] : mix(row[i], row[i + 1], k / n));
  });
  return result;
};
const stationKey = (s, i = s.i) => `${s.kind}/${s.body}/${s.side ?? ''}/${i}`;
const validateHistory = (history, bodies) => {
  assert.ok(history && Array.isArray(history.fractions), 'Missing ISES inlet history.');
  assert.equal(history.fractions.length, bodies.length);
  history.fractions.forEach((row, body) => {
    assert.ok(Array.isArray(row) && row.length === bodies[body].leadingIndex + 1
      && row[0] === 0 && row.at(-1) === 1 && row.every(Number.isFinite)
      && row.every((v, i) => i === 0 || v > row[i - 1]), 'Invalid ISES checkpoint inlet fractions.');
  });
  assert.ok(Array.isArray(history.lastRedistributedStagnation)
    && history.lastRedistributedStagnation.length === bodies.length
    && history.lastRedistributedStagnation.every(Number.isFinite), 'Invalid last redistribution markers.');
  assert.ok(['amd', 'colamd'].includes(history.preferredOrdering) && [.001, 1].includes(history.pivotTolerance)
    && ['convex', 'ises-sampled'].includes(history.iterationGeometry)
    && ['listing', 'admissible'].includes(history.stepAcceptance)
    && ['listing', 'prose'].includes(history.stagnationLimiter), 'Invalid ISES continuation controls.');
};

export function nestedRefinementCheckpoint(parentCheckpoint, refinedRestart, { streamwiseSubdivisions, normalSubdivisions,
  streamwiseInterpolation = 'linear' } = {}) {
  assert.equal(parentCheckpoint?.version, 1, 'Expected a version 1 ISES checkpoint.');
  assert.ok(['linear', 'surface-pchip'].includes(streamwiseInterpolation), 'Unknown nested streamwise interpolation.');
  // Build independent charts: evaluate/rebase must not mutate caller state.
  const original = serialize(parentCheckpoint), a = original.restart, b = serialize(refinedRestart);
  assert.equal(b?.input?.gridSpacing?.streamwiseInterpolation ?? 'linear', streamwiseInterpolation,
    'Refined checkpoint streamwise interpolation does not match its certificate.');
  for (const f of [a, b]) assert.ok(f?.input && f.options && Array.isArray(f.initialEuler?.x)
    && f.initialEuler.x.length && Array.isArray(f.initialBL) && f.initialBL.length
    && Array.isArray(f.initialEuler.undisplacedNodes ?? f.initialEuler.nodes), 'Expected complete supplied restart states; cold initialization is not permitted.');
  validateHistory(original.continuation, a.input.bodies);
  assert.ok(Array.isArray(streamwiseSubdivisions), 'Supply the actual streamwise subdivisions.');
  assert.ok(Array.isArray(normalSubdivisions), 'Supply the actual normal subdivisions.');
  const parent = create(a), child = create(b), old = parent.euler.layout, next = child.euler.layout;
  const streamwise = nestedIntervalMap(old.nx, { subdivisions: streamwiseSubdivisions });
  assert.equal(next.nx, streamwise.retained.at(-1), 'Streamwise subdivisions do not describe the refined grid.');
  assert.equal(normalSubdivisions.length, old.tubes.length);
  const normal = normalSubdivisions.map((row, g) => {
    const map = nestedIntervalMap(old.tubes[g], { subdivisions: row });
    assert.equal(next.tubes[g], map.retained.at(-1), 'Normal subdivisions do not describe the refined grid.');
    return map.retained;
  });
  assert.ok(streamwise.counts.some(n => n > 1) || normalSubdivisions.some(row => row.some(n => n > 1)), 'Expected a refinement.');
  const requireConvex = original.continuation.iterationGeometry === 'convex';
  assert.equal(parent.euler.conditions.stagnationMotion, 'walls-only', 'Expected an ISES continuation chart.');
  assert.equal(parent.euler.conditions.geometryDomain, requireConvex ? 'convex' : 'positive-simple');
  assert.ok(parent.admissible(parent.initial, { requireConvex }), 'Parent checkpoint is inadmissible.');
  assert.ok(child.admissible(child.initial, { requireConvex }), 'Refined checkpoint is inadmissible.');
  const before = parent.evaluate(parent.initial), after = child.evaluate(child.initial);
  assert.deepEqual(before.families, original.families, 'Parent checkpoint does not replay exactly.');
  const correspondence = streamwiseInterpolation === 'surface-pchip' ? createStreamtubeRefinementCoordinates({
    stations: streamwise, bodies: old.bodies, fractions: parent.euler.fractions, weights: b.input.weights,
  }) : null;

  // Refiner output may spell previously implicit defaults explicitly.
  assert.deepEqual(child.conditions, parent.conditions, 'Refinement changed effective coupled physics.');
  assert.deepEqual(child.euler.conditions, parent.euler.conditions, 'Refinement changed effective Euler physics or normalization.');
  assert.equal(child.bl.transitionMode, parent.bl.transitionMode);
  assert.deepEqual(child.bl.trips, parent.bl.trips, 'Refinement moved material trips.');
  const mutableInput = new Set(['weights', 'gridSpacing', 'bodies', 'outerLower', 'outerUpper', 'cutPaths']);
  for (const key of new Set([...Object.keys(a.input), ...Object.keys(b.input)]))
    if (!mutableInput.has(key)) assert.deepEqual(b.input[key], a.input[key], `Changed input: ${key}`);
  assert.equal(b.input.bodies.length, a.input.bodies.length);
  for (let body = 0; body < old.elements; body++) {
    const { leadingIndex: le, trailingIndex: te, surfaceFractions: _a, ...geometryA } = a.input.bodies[body];
    const { leadingIndex, trailingIndex, surfaceFractions: _b, ...geometryB } = b.input.bodies[body];
    assert.deepEqual(geometryB, geometryA, 'Refinement changed the material contour.');
    assert.equal(leadingIndex, streamwise.retained[le], 'Incorrect mapped leading edge.');
    assert.equal(trailingIndex, streamwise.retained[te], 'Incorrect mapped trailing edge.');
    for (const side of ['upper', 'lower']) {
      const expected = correspondence ? correspondence.surfaceFractions[body][side]
        : interpolate(parent.euler.fractions[body][side], streamwise.counts.slice(le, te));
      const actual = child.euler.fractions[body][side];
      assert.equal(actual.length, expected.length);
      expected.forEach((f, k) => assert.ok(Math.abs(actual[k] - f) < 1e-14, 'Changed nested surface fractions.'));
      for (let i = le; i <= te; i++) assert.equal(actual[streamwise.retained[i] - leadingIndex], parent.euler.fractions[body][side][i - le]);
    }
  }
  const pathsA = [a.input.outerLower, a.input.outerUpper, ...a.input.cutPaths];
  const pathsB = [b.input.outerLower, b.input.outerUpper, ...b.input.cutPaths];
  assert.equal(pathsB.length, pathsA.length);
  pathsA.forEach((path, p) => {
    assert.equal(pathsB[p].length, next.nx + 1);
    path.forEach((point, i) => assert.deepEqual(pathsB[p][streamwise.retained[i]], point, 'Changed retained input boundary node.'));
  });

  let weightRelativeError = 0, massRelativeError = 0, nodeError = 0, physicalBLError = 0, blDistanceError = 0;
  let surfaceDistanceError = 0, wakeDistanceError = 0, inletInterpolationError = 0;
  for (const [key, indices] of Object.entries(old.globals)) {
    const from = Array.isArray(indices) ? indices : [indices], target = next.globals[key];
    const to = Array.isArray(target) ? target : [target];
    assert.equal(to.length, from.length);
    from.forEach((column, j) => {
      if (column === null) assert.equal(to[j], null);
      else assert.equal(child.initial[to[j]], parent.initial[column], `Changed global ${key}`);
    });
  }
  assert.deepEqual(after.outer.captured, before.outer.captured, 'Changed physical capture levels.');
  assert.deepEqual(after.outer.stagnation, before.outer.stagnation, 'Changed physical stagnation parameters.');
  before.outer.allocation.groups.forEach((group, g) => group.forEach((tube, j) => {
    let mass = 0;
    for (let k = normal[g][j]; k < normal[g][j + 1]; k++) {
      weightRelativeError = Math.max(weightRelativeError, Math.abs(b.input.weights[g][k] / (a.input.weights[g][j] / normalSubdivisions[g][j]) - 1));
      mass += after.outer.allocation.groups[g][k].massFlow;
    }
    massRelativeError = Math.max(massRelativeError, Math.abs(mass / tube.massFlow - 1));
  }));
  before.outer.nodes.forEach((grid, g) => grid.forEach((row, i) => row.forEach((p, j) => {
    const q = after.outer.nodes[g][streamwise.retained[i]][normal[g][j]];
    nodeError = Math.max(nodeError, Math.hypot(q.x - p.x, q.y - p.y));
  })));
  assert.ok(weightRelativeError < 1e-13, 'Child weights do not partition parent weights equally.');
  assert.ok(massRelativeError < 1e-13, 'Child masses do not conserve parent physical masses.');
  assert.ok(nodeError < 2e-12 * parent.euler.conditions.lengthScale, 'Refinement moved retained physical fluid nodes.');
  // Linear historical fractions assume the new dividing-cut points split
  // each old straight segment in the same proportions. Certify that premise
  // independently of the refiner's name or grid-spacing metadata.
  old.bodies.forEach((body, b) => {
    for (let i = 0; i < body.leadingIndex; i++) for (let k = 1; k < streamwise.counts[i]; k++) {
      const t = k / streamwise.counts[i], ii = streamwise.retained[i] + k;
      for (const [g, j] of [[b, old.tubes[b]], [b + 1, 0]]) {
        const p = before.outer.nodes[g][i][j], q = before.outer.nodes[g][i + 1][j];
        const actual = after.outer.nodes[g][ii][normal[g][j]];
        inletInterpolationError = Math.max(inletInterpolationError, Math.hypot(actual.x - mix(p.x, q.x, t), actual.y - mix(p.y, q.y, t)));
      }
    }
  });
  assert.ok(inletInterpolationError < 2e-12 * parent.euler.conditions.lengthScale, 'Inserted inlet cut nodes do not support linear history interpolation.');

  const childStations = new Map(child.bl.stations.map(s => [stationKey(s), s]));
  let auxiliaryError = 0;
  const auxiliaryChanges = [], phaseChanges = [], finiteWakeDisplacement = [];
  for (const p of parent.bl.stations) {
    const q = childStations.get(stationKey(p, streamwise.retained[p.i]));
    assert.ok(q, 'Missing retained BL station.');
    const v = before.layers.states[p.id], w = after.layers.states[q.id];
    const finiteWake = p.kind === 'wake' && !!parent.euler.baseGeometry?.[p.body];
    for (const key of finiteWake ? ['theta', 'ue'] : ['theta', 'deltaStar', 'ue'])
      physicalBLError = Math.max(physicalBLError, Math.abs(w[key] - v[key]));
    if (finiteWake) {
      assert.ok(Number.isFinite(v.wakeGap) && Number.isFinite(w.wakeGap)
        && v.wakeGap >= 0 && w.wakeGap >= 0, 'Missing physical finite-wake gap.');
      const fluidError = Math.abs((w.deltaStar - w.wakeGap) - (v.deltaStar - v.wakeGap));
      physicalBLError = Math.max(physicalBLError, fluidError);
      finiteWakeDisplacement.push({ body: p.body, parentIndex: p.i, refinedIndex: q.i,
        parentTotal: v.deltaStar, refinedTotal: w.deltaStar, parentGap: v.wakeGap, refinedGap: w.wakeGap,
        totalChange: w.deltaStar - v.deltaStar, gapChange: w.wakeGap - v.wakeGap, fluidError });
    }
    blDistanceError = Math.max(blDistanceError, Math.abs(w.s - v.s));
    if (p.kind === 'surface') surfaceDistanceError = Math.max(surfaceDistanceError, Math.abs(w.s - v.s));
    else wakeDistanceError = Math.max(wakeDistanceError, Math.abs(w.s - v.s));
    auxiliaryError = Math.max(auxiliaryError, Math.abs(w.aux - v.aux));
    const identity = { kind: p.kind, body: p.body, ...(p.side ? { side: p.side } : {}), parentIndex: p.i, refinedIndex: q.i };
    if (w.aux !== v.aux) auxiliaryChanges.push({ ...identity, before: v.aux, after: w.aux });
    if (p.regime !== q.regime) phaseChanges.push({ ...identity, before: p.regime, after: q.regime });
  }
  assert.ok(physicalBLError < 1e-14, 'Refinement changed retained physical BL thickness or edge speed.');
  assert.ok(surfaceDistanceError < 2e-12, 'Refinement changed retained material surface distances.');
  let insertedPhysicalBLError = 0;
  if (correspondence) for (const branch of child.bl.surfaces) {
    const oldBranch = parent.bl.surfaces.find(s => s.body === branch.body && s.side === branch.side);
    const sourceValues = oldBranch.ids.map(id => before.layers.states[id]);
    const le = old.bodies[branch.body].leadingIndex, coordinates = correspondence.surfaceCoordinates[branch.body][branch.side];
    for (const id of branch.ids) {
      const station = child.bl.stations[id], u = coordinates[station.i - streamwise.retained[le]] - le;
      const value = after.layers.states[id];
      for (const key of ['theta', 'deltaStar', 'ue']) {
        let expected;
        if (u < 1) expected = key === 'ue' ? sourceValues[0].ue * value.s / sourceValues[0].s : sourceValues[0][key];
        else {
          const i = Math.min(sourceValues.length - 2, Math.floor(u - 1)), t = u - 1 - i;
          expected = t === 0 ? sourceValues[i][key] : t === 1 ? sourceValues[i + 1][key]
            : mix(sourceValues[i][key], sourceValues[i + 1][key], t);
        }
        insertedPhysicalBLError = Math.max(insertedPhysicalBLError, Math.abs(value[key] - expected));
      }
    }
  }
  assert.ok(insertedPhysicalBLError < 1e-14, 'Refined surface BL values do not follow the certified station correspondence.');
  // Wake distance starts at the solid TE, whereas interpolated wake centers
  // start at the displaced TE. Splitting that first segment can change the
  // cumulative wake distance even with every retained fluid node unchanged.

  const continuation = serialize(original.continuation);
  continuation.fractions = original.continuation.fractions.map((row, body) =>
    interpolate(row, streamwise.counts.slice(0, old.bodies[body].leadingIndex)));
  validateHistory(continuation, next.bodies);
  const checkpoint = { version: 1, families: { ...after.families }, restart: b, continuation };
  const replay = create(serialize(checkpoint).restart), replayed = replay.evaluate(replay.initial);
  assert.deepEqual(replayed.residual, after.residual, 'Serialized refined residual does not replay exactly.');
  assert.deepEqual(replayed.outer.nodes, after.outer.nodes, 'Serialized refined physical grid does not replay exactly.');
  assert.deepEqual(replayed.families, checkpoint.families);
  return { checkpoint, diagnostics: {
    ...(correspondence ? { streamwiseInterpolation, insertedPhysicalBLError,
      stationCorrespondence: 'PCHIP contour fractions and surface BL primary values checked against the shared source-coordinate map; inserted interior geometry and density are admissible guesses, not certified interpolants.' } : {}),
    nodeError, inletInterpolationError, weightRelativeError, massRelativeError, physicalBLError, blDistanceError, surfaceDistanceError, wakeDistanceError,
    auxiliaryError, auxiliaryChanges, phaseChanges,
    ...(finiteWakeDisplacement.length ? { finiteWakeDisplacement } : {}),
    parentTransitionState: parent.bl.snapshotActive(), refinedTransitionState: child.bl.snapshotActive(),
    retainedStreamwiseStations: streamwise.retained, retainedNormalStations: normal,
    parentFamilies: before.families, refinedFamilies: after.families,
    maintenanceHistoryPreserved: true, inletHistoryMapping: 'Linear in each nested parent interval; retained fractions copied exactly.',
    initialSMOVERepeated: false, exactSerializedReplay: true,
    interpretation: 'Transferred guess only. Refiner amplification/phase and density initialization are retained; full coupled reconvergence is required. Preserved redistribution markers use the refined spacing for future ordinary SMOVE triggers.'
  } };
}
