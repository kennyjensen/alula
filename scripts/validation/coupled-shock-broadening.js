// SPDX-License-Identifier: GPL-2.0-or-later
// Validation only. MSES manual section 1.2.7 temporary shock broadening.
// Density change here means the APPLIED/DAMPED physical density change;
// the manual's update convention is ambiguous. No iteration policy is added.
import { createCoupledStreamtubeBody } from '../../src/euler/streamtube-coupled.js';

const require = (ok, message) => { if (!ok) throw new Error(message); };
const serial = v => JSON.stringify(v, (_, x) => ArrayBuffer.isView(x) ? Array.from(x) : x);
const same = (a, b) => serial(a) === serial(b);
const threshold = x => Number.isFinite(x) && x >= .75 && x <= 1;
const maximum = values => Array.from(values).reduce((m, v) => Math.max(m, Math.abs(v)), 0);

export function temporaryShockMcrit({ targetMcrit, densityChange }) {
  require(threshold(targetMcrit), 'Target shock Mcrit must be finite and in [.75,1].');
  require(Number.isFinite(densityChange) && densityChange >= 0, 'Applied density change must be finite and nonnegative.');
  const a = .15 / 4, d = densityChange;
  // Algebraically d^3/[.15*(d^2+a^2)], without the large-d Infinity/Infinity.
  // Overflow of r or r*r is harmless: exp(-Infinity) is its exact zero limit.
  const t = d < a ? d / a : a / d;
  const r = d < a ? (d / .15) * (t * t / (1 + t * t)) : (d / .15) / (1 + t * t);
  return .75 + (targetMcrit - .75) * Math.exp(-r * r);
}

export function maximumAppliedDensityChange(beforeEuler, afterEuler, densityCount) {
  require((Array.isArray(beforeEuler) || ArrayBuffer.isView(beforeEuler))
    && (Array.isArray(afterEuler) || ArrayBuffer.isView(afterEuler))
    && Number.isInteger(densityCount) && densityCount >= 0
    && beforeEuler.length >= densityCount && afterEuler.length >= densityCount,
  'Supply matching packed Euler density prefixes and a nonnegative density count.');
  let maximum = 0;
  for (let i = 0; i < densityCount; i++) {
    require(Number.isFinite(beforeEuler[i]) && Number.isFinite(afterEuler[i]), 'Packed log densities must be finite.');
    const change = Math.abs(Math.expm1(afterEuler[i] - beforeEuler[i]));
    require(Number.isFinite(change), 'Applied relative density change exceeds the finite representable domain.');
    maximum = Math.max(maximum, change);
  }
  return maximum;
}

export function rebaseCoupledMcrit(checkpoint, mcrit) {
  require(threshold(mcrit), 'Temporary shock Mcrit must be finite and in [.75,1].');
  require(checkpoint?.version === 1 && checkpoint.restart && checkpoint.continuation && checkpoint.families,
    'Shock rebasing requires a complete version-1 coupled checkpoint.');
  require(!Object.hasOwn(checkpoint, 'residual'), 'Use a standard coupled checkpoint without a separate residual array.');
  const original = serial(checkpoint), saved = structuredClone(checkpoint), f = saved.restart, c = saved.continuation;
  require(f.input?.streamwiseMode === 'hybrid' && (f.input.flowModel === undefined || f.input.flowModel === 'compressible')
    && f.input.upwind && threshold(f.input.upwind.mcrit) && Number.isFinite(f.input.upwind.mucon)
    && f.input.upwind.boundary?.kind === 'unfiltered-first-two', 'Shock rebasing requires explicit compressible hybrid/upwind controls.');
  require(['listing', 'admissible'].includes(c.stepAcceptance) && ['convex', 'ises-sampled'].includes(c.iterationGeometry)
    && ['listing', 'prose'].includes(c.stagnationLimiter)
    && f.input.geometryDomain === (c.iterationGeometry === 'convex' ? 'convex' : 'positive-simple'),
  'Shock rebasing requires ordinary controls and the checkpoint iteration geometry domain.');
  require(f.initialEuler?.x && f.initialEuler.nodes && f.initialEuler.undisplacedNodes && f.initialBL,
    'Shock rebasing requires complete supplied physical state and chart.');
  const diagnostics = { validationOnly: true, physicalAcceptance: false, sourceMcrit: f.input.upwind.mcrit, targetMcrit: mcrit,
    sameThreshold: mcrit === f.input.upwind.mcrit, requireConvex: c.iterationGeometry === 'convex',
    densityChangeConvention: 'Maximum applied/damped abs(expm1(new log density - old log density)); explicit port interpretation, not an undamped direction norm.',
    operations: { constructors: 0, explicitEvaluations: 0, strictAdmissibilityCalls: 0,
      globalJacobians: 0, globalLinearSolves: 0, globalNewtonUpdates: 0, phaseConversions: 0, nativeMarches: 0 }, stage: 'source replay' };
  const construct = input => {
    diagnostics.operations.constructors++;
    return createCoupledStreamtubeBody(input, { ...f.options, initialEuler: f.initialEuler, initialBL: f.initialBL });
  };
  const evaluate = system => { diagnostics.operations.explicitEvaluations++; return system.evaluate(system.initial); };
  const domain = system => {
    diagnostics.operations.strictAdmissibilityCalls++; let failure;
    const ok = system.admissible(system.initial, { requireConvex: diagnostics.requireConvex, onFailure: d => { failure = d; } });
    if (!ok) throw Object.assign(new Error(`Shock Mcrit rebase ${diagnostics.stage} is inadmissible: ${failure?.message ?? failure?.kind}`), { admissibility: failure });
  };
  try {
    const source = construct(f.input), before = evaluate(source);
    diagnostics.densityCount = source.euler.layout.densityCount; diagnostics.n = source.n; diagnostics.ne = source.ne;
    diagnostics.sourceReplay = { packedExact: same(source.initial, [...f.initialEuler.x, ...f.initialBL]),
      nodesExact: same(before.outer.nodes, f.initialEuler.nodes), familiesExact: same(before.families, saved.families),
      phaseExact: f.options.transitionMode !== 'automatic' || same(source.bl.snapshotActive(), f.options.transitionState),
      fullResidualEvaluated: true, historicalFullResidualCompared: false };
    require(Object.entries(diagnostics.sourceReplay).filter(([k]) => k.endsWith('Exact')).every(([, v]) => v),
      'Shock Mcrit source checkpoint does not replay exactly.');
    require(before.residual.every(Number.isFinite), 'Shock Mcrit source residual must be finite.');
    domain(source);
    diagnostics.sourceFamilies = { ...before.families }; diagnostics.sourceMaximumResidual = maximum(before.residual);
    let result = saved, value = before;
    if (!diagnostics.sameThreshold) {
      diagnostics.stage = 'target threshold';
      const input = { ...structuredClone(f.input), upwind: { ...structuredClone(f.input.upwind), mcrit } };
      const target = construct(input); value = evaluate(target);
      require(value.residual.every(Number.isFinite), 'Shock Mcrit target residual must be finite.');
      domain(target);
      require(target.n === source.n && target.ne === source.ne && same(target.initial, source.initial),
        'Shock Mcrit changed a physical packed unknown.');
      require(same(value.outer.nodes, before.outer.nodes) && same(value.outer.nodes, f.initialEuler.nodes)
        && same(value.outer.undisplacedNodes, before.outer.undisplacedNodes), 'Shock Mcrit changed actual or undisplaced geometry.');
      require(same(target.bl.snapshotActive(), source.bl.snapshotActive()), 'Shock Mcrit changed BL phase.');
      result = { ...saved, families: { ...value.families }, restart: { ...f, input } };
      // These are the only permitted data changes; source arrays and history
      // stay in their original serialized chart, including undisplaced nodes.
      const restored = structuredClone(result); restored.families = saved.families; restored.restart.input.upwind.mcrit = f.input.upwind.mcrit;
      require(same(restored, saved), 'Shock Mcrit changed fields beyond threshold/families.');
      diagnostics.residualMaximumChange = maximum(value.residual.map((r, i) => r - before.residual[i]));
      diagnostics.targetResidualChanges = value.residual.reduce((n, r, i) => n + (r !== before.residual[i] ? 1 : 0), 0);
    } else {
      diagnostics.residualMaximumChange = 0; diagnostics.targetResidualChanges = 0;
    }
    require(serial(checkpoint) === original, 'Caller checkpoint was changed during shock rebasing.');
    diagnostics.targetFamilies = { ...value.families }; diagnostics.targetMaximumResidual = maximum(value.residual);
    diagnostics.physicalStatePreserved = true; diagnostics.phasePreserved = true; diagnostics.maintenancePreserved = true;
    diagnostics.sourceUnchanged = true; diagnostics.stage = 'prepared threshold';
    return { checkpoint: result, residual: value.residual.slice(), diagnostics };
  } catch (error) {
    diagnostics.sourceUnchanged = serial(checkpoint) === original;
    error.shockBroadening = structuredClone(diagnostics); throw error;
  }
}
