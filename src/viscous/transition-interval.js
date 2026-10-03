// SPDX-License-Identifier: GPL-2.0-or-later
// Local transition root and its contribution to a simultaneous BL Jacobian.
// Both laminar and turbulent portions use the existing native-derived kernel.
// Differencing re-solves the local root; it never holds the transition point
// fixed while perturbing thickness, amplification, velocity or coordinates.
import { boundedResidualPartial } from '../numerics/bounded-partial.js';

const keys = ['aux', 'theta', 'deltaStar', 'ue', 's'];

// TRCHEK2 selects min(free fraction, forced fraction), then compares the
// interpolated XT against XIFORC with strict < to label the result. When a
// trip controls XT, roundoff in X1*WF1+X2*WF2 can flip that label. Classify
// the already selected root instead; preserve the original flag separately.
export function transitionRootIsForced(input, location) {
  const trip = input.tripS;
  return Number.isFinite(trip) && trip > input.upstream.s && trip <= input.downstream.s
    && Math.abs(location - trip) <= 8 * Number.EPSILON * Math.max(Math.abs(location), Math.abs(trip));
}

export function checkAutomaticTransition(kernel, input) {
  const natural = kernel.transitionCheck({ ...input, tripS: Number.MAX_VALUE });
  const forced = input.tripS > input.upstream.s && input.tripS <= input.downstream.s
    && (!natural.transition || input.tripS < natural.s);
  if (!forced) return { ...natural, forced: false, nativeForced: natural.forced };
  const raw = kernel.transitionCheck(input);
  if (!raw.transition || !transitionRootIsForced(input, raw.s))
    throw new Error('The forced transition root disagrees with the earlier-trip selection.');
  return { ...raw, forced: true, nativeForced: raw.forced };
}

export function evaluateTransitionInterval(kernel, input, { jacobian = false } = {}) {
  input = { ...input, regime: 'transition', tripS: input.tripS ?? Number.MAX_VALUE };
  const activeValue = candidate => {
    const selected = checkAutomaticTransition(kernel, candidate);
    if (!selected.transition) throw new Error('Transition is outside this active interval.');
    // A later trip must not affect a free root or its derivatives. TRCHEK's
    // finite internal stopping tolerance can otherwise leave a trip-dependent
    // free root when the two locations are close. Use the untripped root when
    // natural transition wins; use the native mixed interval when the trip wins.
    const raw = kernel.interval({ ...candidate, tripS: selected.forced ? candidate.tripS : Number.MAX_VALUE });
    return { ...raw, transition: { ...raw.transition, nativeForced: raw.transition.forced, forced: selected.forced } };
  };
  const value = activeValue(input), forced = value.transition.forced;
  if (!jacobian) return value;
  // Coordinate sensitivity is set by the local interval, not its distance
  // from stagnation. In a short forced-TE interval, steps proportional to s
  // leave substantial one-sided truncation errors in the separate endpoint
  // and trip partials, even when their physical effects should cancel.
  const coordinateStep = Math.cbrt(Number.EPSILON) * (input.downstream.s - input.upstream.s);
  if (!(coordinateStep > 0) || !Number.isFinite(coordinateStep)
    || [input.upstream.s, input.downstream.s, ...(forced ? [input.tripS] : [])]
      .some(s => s + coordinateStep === s || s - coordinateStep === s))
    throw new Error('Unresolved transition-coordinate derivative step.');
  const packed = r => [...r.residual, r.transition.s];
  const base = packed(value);
  // The derivative belongs to this active branch. At a natural/forced switch
  // or interval boundary the other side differentiates different equations.
  const evaluate = candidate => {
    const r = activeValue(candidate);
    if (r.transition.forced !== forced) throw new Error('Transition derivative crossed the natural/forced branch.');
    return packed(r);
  };
  const partial = (f, v, h, lower = -Infinity, upper = Infinity) => {
    // First use the established bounded second-order stencil. For a natural
    // transition the branch boundary also depends on the thermodynamic state;
    // test the samples explicitly instead of guessing its position.
    try { return boundedResidualPartial(f, v, { step: h, lower, upper, base }); }
    catch (cause) {
      for (const sign of [1, -1]) {
        if (v + 2 * sign * h > upper || v + 2 * sign * h < lower) continue;
        try {
          const a = f(v + sign * h), b = f(v + 2 * sign * h);
          const d = a.map((w, k) => (4 * (w - base[k]) - (b[k] - base[k])) / (2 * sign * h));
          if (d.every(Number.isFinite)) return d;
        } catch { /* This side is outside the active kernel domain. */ }
      }
      throw new Error('No resolved second-order derivative inside this transition branch.', { cause });
    }
  };
  // Amplification can place a natural root very close to the downstream
  // endpoint. Its local curvature then makes a second-order stencil too
  // inaccurate even with a cube-root-epsilon relative step. Re-solve the
  // root at four samples on the selected branch; at a branch boundary use
  // four samples on its admissible side. Keep the existing bounded stencil
  // if that wider stencil does not fit the active domain.
  const amplificationPartial = (f, v, h) => {
    for (const sign of [0, 1, -1]) {
      try {
        const offsets = sign ? [1, 2, 3, 4].map(k => sign * k) : [-2, -1, 1, 2];
        if (offsets.some(o => v + o * h === v)) continue;
        const [a, b, c, d] = offsets.map(o => f(v + o * h));
        const result = base.map((w, r) => sign
          ? (48 * (a[r] - w) - 36 * (b[r] - w) + 16 * (c[r] - w) - 3 * (d[r] - w)) / (12 * sign * h)
          : (8 * (c[r] - b[r]) - (d[r] - a[r])) / (12 * h));
        if (result.every(Number.isFinite)) return result;
      } catch { /* Samples must retain this transition branch. */ }
    }
    return partial(f, v, h);
  };
  const derivatives = {};
  for (const side of ['upstream', 'downstream']) {
    const columns = keys.map((key, k) => {
      // Downstream aux is turbulent shear, not laminar amplification. It
      // changes neither TRCHEK's root nor the interpolated transition state.
      // TRDIF therefore has BL2(K,1)=0 and BT2(K,1)=VS2(K,1). Reuse that
      // analytic column instead of differencing tiny Ctau with a .01 floor.
      // Its selected Hk-floor policy is already included in `value`.
      if (side === 'downstream' && key === 'aux') return [...value.downstream.map(row => row[0]), 0];
      const v = input[side][key];
      const h = key === 's' ? coordinateStep
        : Math.cbrt(Number.EPSILON) * Math.max(Math.abs(v), k === 0 ? .01 : k < 3 ? 1e-7 : 1e-6);
      const lower = key === 's' ? (side === 'upstream' ? 0 : forced ? input.tripS : input.upstream.s) : -Infinity;
      const upper = key === 's' && side === 'upstream' ? (forced ? input.tripS : input.downstream.s) : Infinity;
      const f = w => evaluate({ ...input, [side]: { ...input[side], [key]: w } });
      return side === 'upstream' && key === 'aux' ? amplificationPartial(f, v, h)
        : partial(f, v, h, lower, upper);
    });
    derivatives[side] = Array.from({ length: 4 }, (_, r) => columns.map(c => c[r]));
  }
  // Once natural transition precedes a trip, moving that later trip has no
  // effect within this branch. Avoid differencing a sentinel/no-trip value.
  derivatives.trip = forced ? partial(tripS => evaluate({ ...input, tripS }), input.tripS,
    coordinateStep, input.upstream.s, input.downstream.s) : [0, 0, 0, 0];
  if (value.hkFloorLinearizationUsed) {
    // Keep the resolved-root derivatives above. Only the optional native
    // Hk-floor policy adds XFOIL's quasi-Newton continuation through the
    // clipped shape branch; all other native analytic approximations cancel.
    const correction = kernel.interval({ ...input, tripS: forced ? input.tripS : Number.MAX_VALUE },
      { hkFloorCorrection: true }).hkFloorJacobianCorrection;
    if (!correction) throw new Error('Missing native Hk-floor linearization correction.');
    for (const side of ['upstream', 'downstream']) for (let r = 0; r < 3; r++) for (let k = 0; k < keys.length; k++)
      if (!(side === 'downstream' && k === 0)) derivatives[side][r][k] += correction[side][r][k];
    for (let r = 0; r < 3; r++) derivatives.trip[r] += correction.trip[r];
  }
  return { ...value, partials: {
    upstream: derivatives.upstream.slice(0, 3), downstream: derivatives.downstream.slice(0, 3), trip: derivatives.trip.slice(0, 3),
    location: { upstream: derivatives.upstream[3], downstream: derivatives.downstream[3], trip: derivatives.trip[3] },
  } };
}
