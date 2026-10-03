// SPDX-License-Identifier: GPL-2.0-or-later
// Warm-start the same physical body grid without resetting its density field.
// A change of freestream Mach changes h0 and the entropy reference: retaining
// rho and tube mass is a continuation guess, not exact entropy transport
// between different operating points. No isentropic inversion occurs here.
import { createStreamtubeBodySystem } from '../streamtube-body.js';

const same = (a, b) => {
  if (Object.is(a, b)) return true;
  if (!a || !b || typeof a !== 'object' || typeof b !== 'object'
    || Array.isArray(a) !== Array.isArray(b)) return false;
  const keys = Object.keys(a);
  return keys.length === Object.keys(b).length && keys.every(k => Object.hasOwn(b, k) && same(a[k], b[k]));
};
const require = (condition, message) => { if (!condition) throw new Error(message); };
const topology = system => {
  const l = system.layout;
  return { nx: l.nx, n: l.n, densityCount: l.densityCount, tubes: l.tubes,
    primaryBody: l.primaryBody, independentWakeBanks: l.independentWakeBanks,
    globals: l.globals,
    bodies: l.bodies.map(b => ({ points: b.points, leadingIndex: b.leadingIndex,
      trailingIndex: b.trailingIndex, roundedLeadingEdge: b.roundedLeadingEdge,
      ...(b.trailingEdge ? { trailingEdge: b.trailingEdge } : {}) })),
    fractions: system.fractions, initialStagnation: system.initialStagnation,
    originalNodes: system.originalNodes };
};
const entropy = (s, c) => Math.log(c.rhoTotal / s.rho) + Math.log(s.enthalpy / c.h0) / (c.gamma - 1);

export function initializeStreamtubeBodyFromFlow(input, source, { initial = source?.initial } = {}) {
  require(source?.conditions?.flowModel === 'compressible' && source.layout?.densityCount > 0
    && typeof source.evaluate === 'function' && typeof source.decode === 'function'
    && Array.isArray(source.layout.bodies) && Array.isArray(source.fractions)
    && Array.isArray(source.originalNodes) && Array.isArray(source.initialStagnation),
  'Flow restart requires a complete compressible source system with gas conditions and topology.');
  require(initial instanceof Float64Array && initial.length === source.layout.n && initial.every(Number.isFinite),
    'Flow restart requires the complete finite encoded source state.');
  const a = source.conditions;
  require([a.mach, a.gamma, a.h0, a.rhoTotal, a.lengthScale, a.massScale].every(v => Number.isFinite(v) && v > 0)
    && a.gamma > 1 && a.mach < 1, 'Invalid source gas or reference normalization.');
  require(a.h0 === 1 / ((a.gamma - 1) * a.mach * a.mach) + .5
    && a.rhoTotal === (1 + .5 * (a.gamma - 1) * a.mach * a.mach) ** (1 / (a.gamma - 1)),
  'Source conditions do not use the body rhoInfinity=qInfinity=1 normalization.');

  const targetInput = structuredClone(input), system = createStreamtubeBodySystem(targetInput), b = system.conditions;
  require(b.flowModel === 'compressible', 'Flow restart requires a compressible target.');
  require(same(topology(source), topology(system)),
    'Flow restart requires identical topology, material stations, contours and base grid.');
  for (const key of ['alpha', 'gamma', 'lengthScale', 'massScale', 'normalStencil', 'stagnationMotion',
    'geometryDomain', 'wakeGeometry', 'wakeOutlet', 'wakeDisplacementMotion', 'pressureCorrectionFactor', 'streamwiseMode', 'upwind', 'hybrid'])
    require(same(a[key], b[key]), `Flow restart changed ${key}; only freestream Mach may change.`);
  require(same(source.displacement, system.displacement),
    'Flow restart requires identical current displacement and wake thicknesses.');
  const sourceBase = source.decode(new Float64Array(source.layout.n)), targetBase = system.decode(system.initial);
  require(same(sourceBase.captured, targetBase.captured), 'Flow restart changed the captured-mass reference levels.');

  // Evaluate accepted source data without rebasing or changing its chart.
  const before = source.evaluate(initial), state = system.initial.slice();
  state.set(initial.subarray(0, source.layout.densityCount));
  for (const [key, oldColumns] of Object.entries(source.layout.globals)) {
    const from = Array.isArray(oldColumns) ? oldColumns : [oldColumns];
    const to = Array.isArray(system.layout.globals[key]) ? system.layout.globals[key] : [system.layout.globals[key]];
    from.forEach((col, i) => { if (col !== null) state[to[i]] = initial[col]; });
  }
  // With identical references these copied globals preserve physical capture,
  // stagnation parameters and multipoles, including their last stored bits.
  const transferred = system.adoptGeometry(state, before.nodes), flow = system.evaluate(transferred);
  require(same(before.captured, flow.captured) && same(before.stagnation, flow.stagnation)
    && same(before.strengths, flow.strengths), 'Flow restart did not preserve its physical global variables.');
  for (let g = 0; g < source.layout.tubes.length; g++) for (let j = 0; j < source.layout.tubes[g]; j++)
    require(before.allocation.groups[g][j].massFlow === flow.allocation.groups[g][j].massFlow,
      'Flow restart changed a physical streamtube mass.');

  let maximumGeometryChange = 0, maximumSpeedChange = 0, maximumEntropyChange = 0;
  for (let g = 0; g < before.nodes.length; g++) for (let i = 0; i < before.nodes[g].length; i++)
    for (let j = 0; j < before.nodes[g][i].length; j++) {
      const p = before.nodes[g][i][j], q = flow.nodes[g][i][j];
      maximumGeometryChange = Math.max(maximumGeometryChange, Math.hypot(p.x - q.x, p.y - q.y));
    }
  for (let i = 0; i < source.layout.nx; i++) for (let g = 0; g < source.layout.tubes.length; g++)
    for (let j = 0; j < source.layout.tubes[g]; j++) {
      const old = before.sections[i][g][j], next = flow.sections[i][g][j];
      require(old.rho === next.rho, 'Flow restart changed a physical section density.');
      maximumSpeedChange = Math.max(maximumSpeedChange, Math.abs(old.q - next.q));
      maximumEntropyChange = Math.max(maximumEntropyChange, Math.abs(entropy(next, b) - entropy(old, a)));
    }
  const sameMach = a.mach === b.mach;
  return { input: targetInput, system, initial: transferred,
    initialEuler: { x: transferred.slice(), nodes: structuredClone(flow.nodes),
      ...(flow.undisplacedNodes ? { undisplacedNodes: structuredClone(flow.undisplacedNodes) } : {}) },
    flow, diagnostics: { sourceMach: a.mach, targetMach: b.mach, sameMach,
      densityInitialization: 'copied physical density; no isentropic inversion',
      physicalMassPreserved: true, physicalDensityPreserved: true,
      maximumGeometryChange, maximumSpeedChange, maximumEntropyChange,
      sourceH0: a.h0, targetH0: b.h0, sourceRhoTotal: a.rhoTotal, targetRhoTotal: b.rhoTotal,
      entropyInterpretation: sameMach ? 'same-gas entropy replay; coordinate roundoff is reported'
        : 'Mach continuation retains density and mass; relative entropy changes with the target gas',
      converged: false },
    status: 'physical-flow warm start; target equations still require convergence' };
}
