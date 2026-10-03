// SPDX-License-Identifier: GPL-2.0-or-later
// Validation only: one raw Newton proposal, still carrying its OLD phase,
// enters native MRCHDU before any generic transition conversion. No iteration
// policy, matrix, linear solve, acceptance relaxation or production hook.
import { createCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';
import { relaxXfoilProfiles } from '../../src/viscous/xfoil-profile-relaxation.js';
import { extendWarmBoundaryIncrements } from '../../src/euler/streamtube-displacement.js';
import { streamtubeMeshSnapshot } from '../../src/euler/streamtube-mesh-preview.js';
import { streamtubeWakeGap } from '../../src/euler/streamtube-wake-geometry.js';

const copy = v => JSON.parse(JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x));
const same = (a, b) => JSON.stringify(copy(a)) === JSON.stringify(copy(b));
const require = (ok, message) => { if (!ok) throw new Error(message); };
const exact = (a, b, name) => require(same(a, b), `${name} changed.`);
const departure = (a, b) => Math.max(0, ...a.flatMap((g, gi) => g.flatMap((r, i) => r.map((p, j) =>
  Math.hypot(p.x - b[gi][i][j].x, p.y - b[gi][i][j].y)))));

export function prepareRawCoupledProposal(checkpoint, rawPacked, { relaxProfiles = relaxXfoilProfiles, onRecord } = {}) {
  const sourceBefore = copy(checkpoint), rawBefore = copy(rawPacked), cp = copy(checkpoint);
  const report = { validationOnly: true, prepared: false, stage: 'input',
    method: 'raw old-phase proposal -> native complete profile -> boundary increment -> canonical checkpoint',
    operations: { constructors: 0, explicitEvaluations: 0, admissibilityCalls: 0, successfulAdmissibilityEvaluations: 0,
      profileAdapterCalls: 0, genericPhaseConversions: 0, jacobians: 0, linearSolves: 0, newtonUpdates: 0 } };
  report.adapter = { name: relaxProfiles?.name ?? null, injected: relaxProfiles !== relaxXfoilProfiles };
  const records = { rawPacked: rawBefore, sourcePhase: cp?.restart?.options?.transitionState };
  const record = (kind, value) => { records[kind] = copy(value); onRecord?.(kind, copy(value)); };
  const construct = f => {
    report.operations.constructors++;
    return createCoupledStreamtubeBody(f.input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  };
  const evaluate = (system, x) => { report.operations.explicitEvaluations++; return system.evaluate(x); };
  const domains = (system, x) => {
    report.operations.admissibilityCalls++; let failure;
    const ok = system.admissible(x, { requireConvex: true, onFailure: f => { failure = f; } });
    if (!ok) { const e = new Error(`Raw-profile preparation inadmissible: ${failure?.message ?? failure?.kind ?? 'unknown domain failure'}`);
      e.admissibility = failure; throw e; }
    report.operations.successfulAdmissibilityEvaluations++;
  };
  try {
    require(cp?.version === 1 && cp.restart?.initialEuler && cp.restart.initialBL && cp.continuation,
      'A complete supplied coupled checkpoint is required.');
    require(typeof relaxProfiles === 'function', 'Supply a native profile adapter function.');
    const f = cp.restart;
    require(f.options.transitionMode === 'automatic' && f.options.tripFractions?.every(p => p.length === 2 && p.every(t => t === 1)),
      'Raw native preparation supports automatic terminal material trips only.');
    report.stage = 'source chart, no residual';
    const system = construct(f), { bl, euler, ne, n } = system;
    const oldPhase = bl.snapshotActive();
    exact(oldPhase, f.options.transitionState, 'Supplied old phase');
    exact(system.initial, [...f.initialEuler.x, ...f.initialBL], 'Supplied packed chart');
    require(rawPacked?.length === n && Array.from(rawPacked).every(Number.isFinite), 'Raw proposal must contain every finite packed unknown.');
    require(euler.layout.independentWakeBanks, 'Raw native preparation requires independently retained wake banks.');
    require(bl.surfaces.length === 2 * bl.wakes.length && bl.wakes.length > 0 && bl.trips.length === bl.wakes.length,
      'Complete body surfaces and wakes are required.');
    const raw = Float64Array.from(rawPacked), rawEuler = raw.subarray(0, ne), rawBL = raw.subarray(ne);
    report.stage = 'raw geometry, no old-phase residual';
    euler.setDisplacement(bl.thicknesses(rawBL));
    const rawFlow = euler.decode(rawEuler), rawGeometry = bl.geometry(rawEuler);
    const states = bl.stations.map(({ id }) => ({ s: rawGeometry.coordinates[id].s, aux: rawBL[4 * id],
      theta: bl.scale * rawBL[4 * id + 1], deltaStar: bl.scale * rawBL[4 * id + 2], ue: rawBL[4 * id + 3],
      ...(rawGeometry.coordinates[id].wakeGap === undefined ? {} : { wakeGap: rawGeometry.coordinates[id].wakeGap }) }));
    record('raw', { phase: oldPhase, states, nodes: rawFlow.nodes, undisplacedNodes: rawFlow.undisplacedNodes,
      captured: rawFlow.captured, stagnation: rawFlow.stagnation,
      masses: rawFlow.allocation.groups.map(g => g.map(t => t.massFlow)), parameters: bl.kernel.parameters });
    const packed = rawBL.slice(), phases = oldPhase.slice(), covered = new Set(), bodies = [];
    const put = (ids, profiles) => {
      require(Array.isArray(profiles) && profiles.length === ids.length, 'Native profile changed station count.');
      profiles.forEach((p, k) => {
        const id = ids[k]; require(!covered.has(id), 'Native body profiles overlap.'); covered.add(id);
        require(['aux', 'theta', 'deltaStar', 'ue'].every(key => Number.isFinite(p[key])), 'Native profile contains nonfinite physical fields.');
        exact(p.s, states[id].s, 'Native station arclength');
        exact(p.wakeGap ?? 0, states[id].wakeGap ?? 0, 'Native finite-base gap');
        packed.set([p.aux, p.theta / bl.scale, p.deltaStar / bl.scale, p.ue], 4 * id);
      });
    };
    report.stage = 'native old-phase full-body march';
    for (const wake of bl.wakes) {
      const indices = ['upper', 'lower'].map(side => bl.surfaces.findIndex(s => s.body === wake.body && s.side === side));
      require(indices.every(i => i >= 0), 'Native body is missing a surface.');
      const surfaces = indices.map(i => bl.surfaces[i]);
      const input = { surfaces: surfaces.map(s => s.ids.map(id => states[id])), wake: wake.ids.map(id => states[id]),
        phases: indices.map(i => oldPhase[i]), tripS: surfaces.map(s => states[s.ids.at(-1)].s),
        normalGap: states[wake.ids[0]].wakeGap ?? 0 };
      record(`body-${wake.body}-input`, input);
      report.operations.profileAdapterCalls++;
      const output = relaxProfiles(copy(input), { ...bl.kernel.parameters });
      record(`body-${wake.body}-output`, output);
      require(Array.isArray(output.surfaces) && output.surfaces.length === 2, 'Native result is missing a surface.');
      surfaces.forEach((s, k) => {
        const returned = output.surfaces[k], phase = returned.transition;
        require(Number.isInteger(phase) && phase >= 1 && phase < s.ids.length
          && Array.isArray(returned.states) && returned.states.length === s.ids.length,
        'Native output has an unsupported transition interval or profile length.');
        require(returned.states.every((p, j) => Number.isFinite(p.aux)
          && (j < phase ? p.aux >= 0 && p.aux < bl.kernel.parameters.ncrit : p.aux > 0)),
        'Native output auxiliary state disagrees with its returned phase.');
        put(s.ids, returned.states); phases[indices[k]] = phase;
      });
      require(Array.isArray(output.wake) && output.wake.every(p => Number.isFinite(p.aux) && p.aux > 0),
        'Native output wake shear must be positive.');
      put(wake.ids, output.wake);
      bodies.push({ body: wake.body, before: input.phases, after: output.surfaces.map(s => s.transition),
        warnings: output.localConvergenceWarnings ?? [], messages: output.messages ?? [] });
    }
    require(covered.size === bl.stations.length, 'Native output omitted a BL station.');
    record('native', { initialBL: packed, phase: phases, bodies });
    report.stage = 'native phase and boundary increments';
    // No active-set remarch or auxiliary rewrite may replace the native output.
    bl.restoreActive(phases);
    const x = raw.slice(); x.set(packed, ne);
    euler.setDisplacement(bl.thicknesses(packed));
    const decoded = euler.decode(x.subarray(0, ne));
    const masses = rawFlow.allocation.groups.map(g => g.map(t => t.massFlow));
    const moved = extendWarmBoundaryIncrements({ sourceNodes: rawFlow.nodes, targetNodes: decoded.nodes, masses });
    x.set(euler.adoptGeometry(x.subarray(0, ne), moved));
    const adopted = euler.decode(x.subarray(0, ne)), geometry = bl.geometry(x.subarray(0, ne));
    for (const key of ['captured', 'stagnation']) exact(adopted[key], rawFlow[key], `Raw ${key}`);
    exact(adopted.allocation.groups.map(g => g.map(t => t.massFlow)), masses, 'Raw physical tube masses');
    for (let col = 0; col < euler.layout.densityCount; col++) exact(x[col], raw[col], 'Raw log density');
    for (const col of Object.values(euler.layout.globals).flat(Infinity).filter(Number.isInteger)) exact(x[col], raw[col], 'Raw global unknown');
    const coordinateTolerance = 1e-12 * euler.conditions.lengthScale;
    for (const wake of bl.wakes) for (const id of wake.ids.slice(1)) {
      const i = bl.stations[id].i, b = wake.body;
      require(Math.hypot(adopted.nodes[b][i].at(-1).x - rawFlow.nodes[b][i].at(-1).x,
        adopted.nodes[b][i].at(-1).y - rawFlow.nodes[b][i].at(-1).y) <= coordinateTolerance, 'Lower independent wake bank changed.');
      require(Math.hypot(adopted.nodes[b + 1][i][0].x - rawFlow.nodes[b + 1][i][0].x,
        adopted.nodes[b + 1][i][0].y - rawFlow.nodes[b + 1][i][0].y) <= coordinateTolerance, 'Upper independent wake bank changed.');
      require(Math.abs((geometry.coordinates[id].wakeGap ?? 0) - (rawGeometry.coordinates[id].wakeGap ?? 0)) <= coordinateTolerance / euler.conditions.lengthScale,
        'Native wake gap no longer matches adopted physical geometry.');
    }
    report.geometryExtension = { maximumInteriorChange: departure(decoded.nodes, adopted.nodes), coordinateTolerance };
    report.maximumPhysicalDisplacementChange = Math.max(...states.map((s, id) => Math.abs((packed[4 * id + 2] - rawBL[4 * id + 2]) * bl.scale * euler.conditions.lengthScale)));
    report.displacementUnits = 'physical solver geometry length; packed deltaStar * BL scale * solver length';
    const options = { ...f.options, transitionState: phases };
    let candidate = { ...cp, restart: { input: f.input, options,
      initialEuler: { x: x.slice(0, ne), nodes: adopted.nodes, undisplacedNodes: adopted.undisplacedNodes }, initialBL: packed } };
    record('candidate-before-evaluation', candidate);
    report.stage = 'native phase eligibility';
    const targets = bl.activeTargets(x.subarray(0, ne), x.subarray(ne));
    record('native-phase-targets', targets);
    require(targets.every(t => t.from === t.to && !t.amplificationReconciliation),
      'Native phase cannot hand off unchanged: simultaneous transition selection or packed amplification differs.');
    report.stage = 'strict prepared domains'; domains(system, x);
    const value = evaluate(system, x);
    require(value.residual.every(Number.isFinite), 'Prepared residual is nonfinite.');
    report.quality = streamtubeMeshSnapshot({ system: euler, nodes: value.outer.nodes }).quality;
    require(report.quality.valid, 'Prepared grid is not strictly convex.');
    candidate = { ...candidate, families: value.families };
    report.stage = 'serialized canonical replay';
    const serialized = copy(candidate), check = construct(serialized.restart);
    exact(check.initial, x, 'Serialized packed proposal'); exact(check.bl.snapshotActive(), phases, 'Serialized native phase');
    domains(check, check.initial); const replay = evaluate(check, check.initial);
    exact(replay.residual, value.residual, 'Serialized complete residual');
    exact(replay.outer.nodes, value.outer.nodes, 'Serialized physical nodes');
    exact(serialized.continuation, cp.continuation, 'Source maintenance history');
    report.wakeGapMismatch = bl.wakes.map(w => ({ body: w.body, rows: w.ids.slice(1, -1).map(id => {
      const i = bl.stations[id].i, b = w.body, indices = [i - 1, i, i + 1];
      const actual = streamtubeWakeGap(indices.map(k => value.outer.nodes[b][k].at(-1)), indices.map(k => value.outer.nodes[b + 1][k][0])).gap;
      const prescribed = value.layers.states[id].deltaStar * euler.conditions.lengthScale;
      return { id, i, actual, prescribed, actualMinusPrescribed: actual - prescribed };
    }) }));
    report.afterFamilies = value.families; report.nativePhase = phases;
    report.canonicalReplayExact = true; report.stage = 'prepared'; report.prepared = true;
    exact(checkpoint, sourceBefore, 'Caller source checkpoint'); exact(rawPacked, rawBefore, 'Caller raw proposal');
    report.sourceUnchanged = true; report.rawUnchanged = true;
    record('prepared', { checkpoint: serialized, report, residual: replay.residual });
    return { checkpoint: serialized, report, residual: replay.residual.slice(), records };
  } catch (error) {
    error.preparation = { report: copy(report), records: copy(records),
      sourceUnchanged: same(checkpoint, sourceBefore), rawUnchanged: same(rawPacked, rawBefore) };
    throw error;
  }
}
