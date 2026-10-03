// SPDX-License-Identifier: GPL-2.0-or-later
// Process-local geometry handoff. Gas is initialized by the ordinary startup
// after reuse; this helper neither evaluates flow nor changes its equations.
import { createStreamtubeBodySystem } from './streamtube-body.js';
import { streamtubeMeshSnapshot } from './streamtube-mesh-preview.js';
import { streamtubeEquationControls } from './streamtube-equation-selection.js';

const captures = new WeakMap();
const require = (ok, message) => { if (!ok) throw new Error(message); };
const copy = structuredClone;
const same = (a, b) => {
  if (Object.is(a, b)) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object') return false;
  const ka = Object.keys(a), kb = Object.keys(b);
  return ka.length === kb.length && ka.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};
const dataOf = ({ system, ...data }) => data;
const withoutMach = ({ mach, ...rest }) => rest;
const gasKeys = new Set(['mach', 'pInf', 'pressureScale', 'h0', 'rhoTotal']);
const otherConditions = value => Object.fromEntries(Object.entries(value).filter(([k]) => !gasKeys.has(k)));
const vector = value => (Array.isArray(value) || ArrayBuffer.isView(value))
  && value.length > 0 && Array.from(value).every(Number.isFinite);

function nodeDeparture(a, b, lengthScale) {
  require(Array.isArray(a) && Array.isArray(b) && a.length === b.length, 'Prepared Euler node groups changed.');
  let maximum = 0;
  a.forEach((group, g) => {
    require(Array.isArray(b[g]) && group.length === b[g].length, 'Prepared Euler node stations changed.');
    group.forEach((row, i) => {
      require(Array.isArray(b[g][i]) && row.length === b[g][i].length, 'Prepared Euler node tubes changed.');
      row.forEach((p, j) => {
        const q = b[g][i][j];
        require(p && q && [p.x, p.y, q.x, q.y].every(Number.isFinite), 'Nonfinite prepared Euler node.');
        const d = Math.hypot(p.x - q.x, p.y - q.y);
        const scale = Math.max(lengthScale, Math.abs(p.x), Math.abs(p.y), Math.abs(q.x), Math.abs(q.y));
        require(d <= 64 * Number.EPSILON * scale, 'Prepared Euler reuse changed physical geometry beyond roundoff.');
        maximum = Math.max(maximum, d);
      });
    });
  });
  return maximum;
}

function validate(prepared, caseData) {
  const s = prepared?.system, mach = caseData?.mach ?? .2;
  require(s?.conditions?.flowModel === 'compressible' && typeof s.decode === 'function'
    && typeof s.geometryChart === 'function' && prepared.input && !prepared.experimental
    && (!s.layout.displacedBoundaries || s.inviscidBaseWake === true), 'Supply a complete prepared inviscid Euler mesh.');
  require(Number.isFinite(mach) && mach > 0 && mach < 1 && prepared.mach === mach
    && prepared.input.mach === mach && s.conditions.mach === mach, 'Prepared Euler case, input and gas Mach disagree.');
  require(prepared.referenceChord === (caseData.referenceChord ?? 1)
    && same(prepared.momentReference, caseData.momentReference ?? { x: prepared.referenceChord / 4, y: 0 })
    && s.conditions.alpha === (caseData.alpha ?? 0), 'Prepared Euler reference conditions disagree with the case.');
  const equations = streamtubeEquationControls(caseData.eulerIsmom);
  for (const [key, value] of Object.entries(equations))
    require(same(prepared.input[key], value) && same(s.conditions[key], value), 'Prepared Euler equation selection disagrees with the case.');
  if (caseData.eulerIsmom === undefined)
    require(prepared.input.streamwiseMode === 'isentropic' && s.conditions.streamwiseMode === 'isentropic'
      && s.conditions.hybrid === undefined && s.conditions.upwind === undefined, 'Prepared Euler default equations changed.');
  require(vector(prepared.initial) && prepared.initial.length === s.layout.n
    && same(prepared.initialEuler?.x, prepared.initial), 'Prepared Euler packed seed is incomplete or stale.');
  require(prepared.mesh?.quality?.valid === true && !prepared.mesh.flow && !prepared.mesh.iteration
    && prepared.mesh.initialization?.flowSolved !== true && !prepared.diagnostics?.gasInitialization,
  'Prepared Euler reuse requires an unsolved initial mesh, not an old flow result.');
  const decoded = s.decode(prepared.initial);
  nodeDeparture(decoded.nodes, prepared.nodes, s.conditions.lengthScale);
  nodeDeparture(decoded.nodes, prepared.initialEuler.nodes, s.conditions.lengthScale);
  return decoded;
}

// Capture only at the preparation seam. The opaque identity and independent
// bindings prevent relabeling a packet or reusing a subsequently mutated chart.
export function capturePreparedStreamtubeAssembly(prepared, caseData) {
  const decoded = validate(prepared, caseData), data = copy(dataOf(prepared));
  const packet = { version: 1, caseData: copy(caseData), prepared: { ...copy(data), system: prepared.system } };
  captures.set(packet, { caseData: copy(caseData), data, system: prepared.system,
    conditions: copy(prepared.system.conditions), chart: copy(prepared.system.geometryChart()),
    fractions: copy(prepared.system.fractions), displacement: copy(prepared.system.displacement), decoded: copy(decoded) });
  return packet;
}

function boundCapture(packet) {
  const saved = captures.get(packet);
  require(saved && packet.version === 1 && packet.prepared?.system === saved.system
    && same(packet.caseData, saved.caseData) && same(dataOf(packet.prepared), saved.data),
  'Prepared Euler packet is unrecognized or its saved data changed.');
  require(same(saved.system.conditions, saved.conditions) && same(saved.system.geometryChart(), saved.chart)
    && same(saved.system.fractions, saved.fractions) && same(saved.system.displacement, saved.displacement)
    && same(saved.system.decode(saved.data.initial), saved.decoded), 'Prepared Euler source system changed after capture.');
  return saved;
}

export function restorePreparedStreamtubeAssembly(packet, caseData) {
  const saved = boundCapture(packet);
  require(same(caseData, saved.caseData), 'Prepared Euler packet belongs to a different case; retarget Mach explicitly.');
  return { ...copy(saved.data), system: saved.system };
}

export function retargetPreparedStreamtubeAssembly(packet, nextCase) {
  const saved = boundCapture(packet), mach = nextCase?.mach;
  require(Number.isFinite(mach) && mach > 0 && mach < 1
    && same(withoutMach(nextCase), withoutMach(saved.caseData)), 'Prepared Euler retargeting may change only Mach.');
  const source = saved.system, input = { ...copy(saved.data.input), mach }, system = createStreamtubeBodySystem(input);
  require(system.layout.n === source.layout.n && system.layout.nx === source.layout.nx
    && same(system.layout.tubes, source.layout.tubes) && same(system.layout.globals, source.layout.globals)
    && same(otherConditions(system.conditions), otherConditions(saved.conditions))
    && same(system.fractions, saved.fractions) && same(system.displacement, saved.displacement),
  'Prepared Euler retargeting changed geometry, normalization or equation controls.');
  const initial = system.adoptGeometry(Float64Array.from(saved.data.initial), saved.decoded.nodes);
  for (let i = 0; i < system.layout.densityCount; i++)
    require(Object.is(initial[i], saved.data.initial[i]), 'Prepared Euler retargeting changed the density seed.');
  for (let i = system.layout.globalOffset; i < system.layout.n; i++)
    require(Object.is(initial[i], saved.data.initial[i]), 'Prepared Euler retargeting changed a global unknown.');
  const decoded = system.decode(initial), maximumNodeDeparture = nodeDeparture(decoded.nodes, saved.decoded.nodes, system.conditions.lengthScale);
  require(same(decoded.captured, saved.decoded.captured) && same(decoded.allocation, saved.decoded.allocation)
    && same(decoded.stagnation, saved.decoded.stagnation), 'Prepared Euler retargeting changed capture, tube masses or stagnation.');
  const { gamma, pInf, h0, rhoTotal } = system.conditions;
  require(pInf === 1 / (gamma * mach * mach) && h0 === 1 / ((gamma - 1) * mach * mach) + .5
    && rhoTotal === (1 + .5 * (gamma - 1) * mach * mach) ** (1 / (gamma - 1)), 'Prepared Euler retargeting retained stale gas constants.');
  const diagnostics = { ...copy(saved.data.diagnostics), preparedEulerReuse: {
    geometryPreparationMach: saved.data.diagnostics?.preparedEulerReuse?.geometryPreparationMach ?? saved.conditions.mach,
    sourceMach: saved.conditions.mach, actualMach: mach, maximumNodeDeparture,
    equationsChanged: false, initialGuessOnly: true, meshRebuilt: false, gasInitialized: false } };
  const mesh = streamtubeMeshSnapshot({ system, initial, nodes: decoded.nodes, diagnostics });
  require(mesh.quality.valid, 'Retargeted prepared Euler mesh failed its convexity gate.');
  const prepared = { ...copy(saved.data), input, system, initial, nodes: decoded.nodes,
    initialEuler: { x: initial.slice(), nodes: copy(decoded.nodes) }, diagnostics, mesh, mach };
  return capturePreparedStreamtubeAssembly(prepared, nextCase);
}
