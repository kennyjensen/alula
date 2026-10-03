// SPDX-License-Identifier: GPL-2.0-or-later
// MSES manual §2.3.3 stops on RMS changes and maxima below 10x the RMS
// limits. EPS.INC is unavailable: these are our documented, dimensionless
// defaults, not claimed MSES constants. Keep strict equation-root stopping
// available for verification and never infer convergence from damping alone.
export const COUPLED_CHANGE_TOLERANCE = 1e-4;
export const COUPLED_CHANGE_MAXIMUM = 10 * COUPLED_CHANGE_TOLERANCE;
const laminar = regime => ['laminar', 'similarity'].includes(regime);
const pressureRows = new Set(['streamwise', 'internalPressure', 'farfieldPressure', 'cutPressure', 'trailingKutta', 'leadingKutta']);

function samples() {
  const groups = {};
  return {
    add(name, v) {
      const g = groups[name] ??= { count: 0, squared: 0, maximum: 0 };
      g.count++; g.squared += v * v; g.maximum = Math.max(g.maximum, Math.abs(v));
    },
    finish() { return Object.fromEntries(Object.entries(groups).map(([name, g]) =>
      [name, { rms: Math.sqrt(g.squared / g.count), maximum: g.maximum, count: g.count }])); },
  };
}
export function smallCoupledChanges(groups) {
  return !!groups && Object.keys(groups).length > 0 && Object.values(groups).every(g =>
    Number.isInteger(g.count) && g.count > 0 && Number.isFinite(g.rms) && g.rms >= 0
    && Number.isFinite(g.maximum) && g.maximum >= 0
    && g.rms <= COUPLED_CHANGE_TOLERANCE && g.maximum <= COUPLED_CHANGE_MAXIMUM);
}

export function coupledResidualChanges(system, residual) {
  const s = samples();
  residual.forEach((v, i) => {
    const family = i < system.ne ? system.euler.layout.rows[i].kind : (i - system.ne) % 4 === 3 ? 'edgeMatching' : 'boundaryLayer';
    // Match the line-search pressure units so low Mach cannot hide imbalance.
    const weight = i < system.ne && system.conditions.edgeMatching === 'pressure' && pressureRows.has(family)
      ? system.euler.conditions.pressureScale : 1;
    s.add(family, weight * v);
  });
  return s.finish();
}

export function coupledChangeSnapshot(system, state, value) {
  return { state: state.slice(), nodes: value.outer.nodes, layers: value.layers.states,
    regimes: system.bl.stations.map(s => s.regime), phase: JSON.stringify(system.bl.snapshotActive()),
    mcrit: system.euler.conditions.upwind?.mcrit, length: system.euler.conditions.lengthScale };
}

// Measure the actual maintained grid, not its rebased normal coordinates.
export function acceptedCoupledChanges(system, before, state, value) {
  const s = samples(), { layout } = system.euler;
  for (let k = 0; k < layout.densityCount; k++) s.add('density', Math.expm1(state[k] - before.state[k]));
  value.outer.nodes.forEach((group, g) => group.forEach((row, i) => row.forEach((p, j) => {
    const q = before.nodes[g][i][j]; s.add('grid', Math.hypot(p.x - q.x, p.y - q.y) / before.length);
  })));
  for (let k = layout.globalOffset; k < system.ne; k++) s.add('global', state[k] - before.state[k]);
  system.bl.stations.forEach(({ id, regime }) => {
    const a = before.layers[id], b = value.layers.states[id];
    const da = a.deltaStar - (a.wakeGap ?? 0), db = b.deltaStar - (b.wakeGap ?? 0);
    s.add('theta', (b.theta - a.theta) / a.theta);
    s.add('displacement', (db - da) / da);
    s.add('edgeSpeed', b.ue - a.ue); // Ue / Uinf
    s.add('shape', db / b.theta - da / a.theta);
    s.add(laminar(regime) ? 'amplification' : 'shear', (b.aux - a.aux) / (laminar(regime) ? 10 : a.aux));
  });
  return s.finish();
}

// Original certified Newton correction, before any step limiter, projection,
// logarithmic trial map, or regularized fallback replaces it. Theta/delta
// scaling cancels; finite-base gaps must be removed and differentiated.
export function newtonCoupledChanges(system, state, direction, value) {
  const s = samples(), { layout } = system.euler, scale = system.bl.scale;
  for (let k = 0; k < layout.densityCount; k++) s.add('density', direction[k]);
  for (const p of layout.positions) s.add('grid', direction[p.column]);
  for (let k = layout.globalOffset; k < system.ne; k++) s.add('global', direction[k]);
  const geo = system.bl.geometry(state.subarray(0, system.ne), true);
  system.bl.stations.forEach(({ id, regime }) => {
    const k = system.ne + 4 * id, a = value.layers.states[id];
    let gapChange = 0;
    for (const [col, d] of geo.coordinates[id].wakeGapDerivatives ?? []) gapChange += d * direction[col];
    const fluid = a.deltaStar - (a.wakeGap ?? 0), dt = scale * direction[k + 1], dd = scale * direction[k + 2] - gapChange;
    s.add('theta', dt / a.theta); s.add('displacement', dd / fluid);
    s.add('edgeSpeed', direction[k + 3]);
    s.add('shape', (dd - fluid / a.theta * dt) / a.theta);
    s.add(laminar(regime) ? 'amplification' : 'shear', direction[k] / (laminar(regime) ? 10 : a.aux));
  });
  return s.finish();
}

export function coupledChangeDecision({ accepted, newton, residual, families, stable }) {
  const passed = stable && smallCoupledChanges(accepted) && smallCoupledChanges(newton) && smallCoupledChanges(residual);
  return { method: 'mses-changes', rmsTolerance: COUPLED_CHANGE_TOLERANCE, maximumTolerance: COUPLED_CHANGE_MAXIMUM,
    converged: Boolean(passed),
    stable, accepted, newton, residual, families: families && { ...families } };
}

// Public consumers must distinguish practical convergence from a strict root.
// Check the evidence as well as the flag; a stale/nonfinite residual fails shut.
export function coupledConvergenceSatisfied(result, tolerance) {
  const families = result?.families, c = result?.convergence;
  if (!families || Object.keys(families).length !== 3
    || !['euler', 'boundaryLayer', 'edgeMatching'].every(k => Number.isFinite(families[k]) && families[k] >= 0)) return false;
  if (Number.isFinite(tolerance) && tolerance > 0 && Object.values(families).every(v => v <= tolerance)) return true;
  return c?.method === 'mses-changes' && c.converged === true && c.stable === true
    && smallCoupledChanges(c.accepted) && smallCoupledChanges(c.newton) && smallCoupledChanges(c.residual)
    && Object.entries(families).every(([k, v]) => v <= COUPLED_CHANGE_MAXIMUM && v === c.families?.[k]);
}
