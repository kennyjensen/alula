// SPDX-License-Identifier: GPL-2.0-or-later
// Density/stagnation part of Giles's UPDATE, thesis pp.260–264. This is
// deliberately separate from SMOVE/DEKINK, which are not implemented here.
// The existing Jacobian uses z=log(rho): delta_z is delta_rho/rho in the
// linear system. log1p(r*delta_z) applies an additive physical-density step
// without replacing the state representation or changing any equation.
export function proposeDensityNewton(system, state, direction, { stagnationLimiter = 'listing', maximumStep = 1 } = {}) {
  if (!['listing', 'prose'].includes(stagnationLimiter)) throw new Error('Unknown source stagnation limiter.');
  if (!Number.isFinite(maximumStep) || maximumStep <= 0 || maximumStep > 1) throw new Error('Invalid maximum density-Newton step.');
  // The prose p.71 and appendix p.260 disagree. Keep the listed 1.5 DSLE
  // default, with the half-spacing prose interpreted using the same DSLE
  // measure only for an explicit comparison, not as an undocumented tuning.
  const capFactor = stagnationLimiter === 'listing' ? 1.5 : .5;
  const { layout, curves } = system, { densityCount, globals } = layout;
  if (!densityCount || state.length !== layout.n || direction.length !== layout.n
    || !Array.from(state).every(Number.isFinite) || !Array.from(direction).every(Number.isFinite))
    throw new Error('Invalid compressible density-Newton proposal.');
  let step = maximumStep, limiter = { kind: maximumStep === 1 ? 'full-step' : 'maximum-step' };
  const limit = (value, cause) => {
    if (!(value > 0) || !Number.isFinite(value)) throw new Error('Invalid density-Newton limiter.');
    if (value < step) { step = value; limiter = cause; }
  };
  for (let col = 0; col < densityCount; col++) {
    const change = direction[col];
    if (change > 1) limit(1 / change, { kind: 'density-increase', column: col });
    if (change < -.5) limit(-.5 / change, { kind: 'density-decrease', column: col });
  }
  const { nodes, captured } = system.decode(state), stagnation = [];
  // Capture unknowns move dividing streamlines in units of total mass.
  // Bound each passage's decrease, including both moving endpoints, before
  // constructing a candidate: decoding an inverted capture interval fails
  // before the geometry/residual line search can reduce its step.
  if (globals.capture?.some(col => col !== null)) {
    const massScale = captured.at(-1) - captured[0];
    const changes = captured.map((_, i) => {
      const col = globals.capture[i - 1];
      return col === null || col === undefined ? 0 : massScale * direction[col];
    });
    for (let group = 0; group < captured.length - 1; group++) {
      const mass = captured[group + 1] - captured[group];
      const change = changes[group + 1] - changes[group];
      if (change < 0) limit(.5 * mass / -change, { kind: 'passage-mass-decrease', group });
    }
  }
  globals.stagnation.forEach((col, body) => {
    if (col === null) return;
    const i = layout.bodies[body].leadingIndex + 1;
    const lower = nodes[body][i].at(-1), upper = nodes[body + 1][i][0];
    // Isolated-airfoil coordinates already share the two surface banks;
    // there is no cascade pitch to remove (Giles p.260).
    const spacing = .5 * Math.hypot(upper.x - lower.x, upper.y - lower.y);
    if (!(spacing > 0) || !Number.isFinite(spacing)) throw new Error('Invalid stagnation-bank spacing.');
    const change = curves[body].length * direction[col];
    if (!Number.isFinite(change)) throw new Error('Nonfinite stagnation correction.');
    if (Math.abs(change) > capFactor * spacing)
      limit(capFactor * spacing / Math.abs(change), { kind: 'stagnation-motion', body, column: col });
    stagnation.push({ body, change, spacing, cap: capFactor * spacing });
  });
  const x = Float64Array.from(state, (v, col) => v + (col < densityCount
    ? Math.log1p(step * direction[col]) : step * direction[col]));
  if (!x.every(Number.isFinite)) throw new Error('Nonfinite density-Newton candidate.');
  const norms = values => {
    let maximum = 0;
    for (const v of values) maximum = Math.max(maximum, Math.abs(v));
    const rms = maximum ? maximum * Math.sqrt(values.reduce((s, v) => s + (v / maximum) ** 2, 0) / values.length) : 0;
    return { rms, maximum, count: values.length };
  };
  return { x, step, stepKind: 'density-newton', limiter, stagnationLimiter,
    undampedUpdate: {
      relativeDensity: norms(Array.from(direction).slice(0, densityCount)),
      normalOverReferenceLength: norms(layout.positions.map(p => direction[p.column])),
    }, stagnation: stagnation.map(s => ({ ...s, acceptedChange: step * s.change })) };
}

// Giles UPDATE pp.261–262: relative Ue/theta/displacement limits and the
// conditional H>1 safeguard. Modern XFOIL aux slots retain their native
// admissibility and event transfers; no historical aux limiter is invented.

export function proposeCoupledDensityNewton(system, state, direction, controls = {}) {
  if (state.length !== system.n || direction.length !== system.n
    || !state.every(Number.isFinite) || !direction.every(Number.isFinite))
    throw new Error('Invalid coupled density-Newton state or direction.');
  let viscousStep = 1, viscousLimiter = { kind: 'full-step' };
  const limit = (value, cause) => {
    if (!(value > 0) || !Number.isFinite(value)) throw new Error('Invalid viscous Newton limiter.');
    if (value < viscousStep) { viscousStep = value; viscousLimiter = cause; }
  };
  for (const { id } of system.bl.stations) {
    const k = system.ne + 4 * id;
    // Follow the source order U, theta, displacement, then H. The H guard
    // applies only if the current viscous proposal crosses delta*=theta.
    for (const [slot, variable] of [[3, 'edge-speed'], [1, 'theta'], [2, 'delta-star']]) {
      if (!(state[k + slot] > 0)) throw new Error('Nonpositive initial viscous variable.');
      const relative = direction[k + slot] / state[k + slot];
      if (viscousStep * relative > 1.7) limit(1.7 / relative, { kind: 'bl-increase', station: id, variable });
      if (viscousStep * relative < -.6) limit(-.6 / relative, { kind: 'bl-decrease', station: id, variable });
    }
    const gap = state[k + 2] - state[k + 1], decrease = direction[k + 1] - direction[k + 2];
    if (!(gap > 0)) throw new Error('Invalid initial viscous shape.');
    if (gap < viscousStep * decrease) limit(.75 * gap / decrease, { kind: 'bl-shape', station: id });
  }
  // DSLE must be measured on the accepted displaced geometry, not on a
  // rejected candidate left behind by an earlier residual evaluation.
  system.euler.setDisplacement(system.bl.thicknesses(state.subarray(system.ne)));
  const maximumStep = controls.maximumStep ?? 1;
  if (!Number.isFinite(maximumStep) || maximumStep <= 0 || maximumStep > 1) throw new Error('Invalid maximum coupled Newton step.');
  const density = proposeDensityNewton(system.euler, state.subarray(0, system.ne), direction.subarray(0, system.ne),
    { ...controls, maximumStep: Math.min(maximumStep, viscousStep) });
  const x = state.map((v, col) => v + density.step * direction[col]); x.set(density.x);
  const limiter = viscousStep <= maximumStep && density.step === viscousStep && viscousStep < 1 ? viscousLimiter : density.limiter;
  return { ...density, x, limiter, stepKind: 'coupled-density-newton', viscousStep, viscousLimiter };
}
